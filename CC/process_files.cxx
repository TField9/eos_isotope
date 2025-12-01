#pragma once
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraph.h"
#include <TH3D.h>
#include "TH2D.h"
#include "TF1.h"
#include "TString.h"
#include <map>
#include <unordered_map>
#include <TString.h>
#include <vector>

Double_t simplified_expgausexp(Double_t* xx, Double_t* par) {
    Double_t x = xx[0];
    Double_t A = par[0];     // 振幅
    Double_t x0 = par[1];    // 峰值位置
    Double_t sigma = par[2]; // 主sigma
    Double_t alphaL = par[3];// 左尾参数
    Double_t alphaR = par[4];// 右尾参数
    
    // 保护参数
    sigma = TMath::Max(sigma, 0.01);
    alphaL = TMath::Max(alphaL, 0.1);
    alphaR = TMath::Max(alphaR, 0.1);
    
    // 定义函数值
    Double_t value = 0.0;
    
    if (x < x0) {
        Double_t t = (x0 - x) / sigma;
        if (t > alphaL) {
            Double_t a = 0.5 * alphaL * alphaL;
            value = TMath::Exp(-a - alphaL * (t - alphaL));
        } else {
            value = TMath::Exp(-0.5 * t * t);
        }
    } else {
        Double_t t = (x - x0) / sigma;
        if (t > alphaR) {
            Double_t a = 0.5 * alphaR * alphaR;
            value = TMath::Exp(-a - alphaR * (t - alphaR));
        } else {
            value = TMath::Exp(-0.5 * t * t);
        }
    }
    
    return A * value;
}

void Fit(TH3D *h1, TH2D *&gSigma3D, TH2D *&gMean3D, const std::string& baseName,TFile *f)
{
    // 能量分箱定义
    const int nkbin = 73;
    const double Be_bins[] = {
        0.08, 0.13, 0.17, 0.21, 0.27, 0.33, 0.41, 0.49, 0.59,
        0.70, 0.82, 0.96, 1.11, 1.28, 1.47, 1.68, 1.91, 2.16,
        2.44, 2.73, 3.06, 3.41, 3.79, 4.20, 4.65, 5.14, 5.64,
        6.18, 6.78, 7.42, 8.12, 8.86, 9.66, 10.51, 11.45, 12.45,
        13.50, 14.65, 15.84, 17.14, 18.54, 20.04, 21.64, 23.34, 25.19,
        27.13, 29.23, 31.48, 33.93, 36.53, 39.33, 42.33, 45.58, 49.08,
        53.08, 57.08, 61.58, 66.58, 72.57, 79.07, 86.57, 95.07, 104.57,
        115.57, 128.57, 144.57, 164.07, 188.57, 219.57, 261.57, 329.07, 439.07,
        649.07, 1649.07
    };

    // 方向余弦分箱定义
    const int n_costheta_bin = 4;
    Double_t costheta_bins[n_costheta_bin+1];
    const Double_t ct_start = 0.94;
    const Double_t ct_end = 1.0;
    const Double_t ct_step = (ct_end - ct_start) / n_costheta_bin;
    
    for(int i = 0; i <= n_costheta_bin; i++) {
        costheta_bins[i] = ct_start + i * ct_step;
    }

    // 创建输出直方图
    gSigma3D = new TH2D(Form("gSig_%s", baseName.c_str()),
                        Form("gSig_%s", baseName.c_str()),
                        nkbin, Be_bins, n_costheta_bin, costheta_bins);
    
    gMean3D = new TH2D(Form("gMean_%s", baseName.c_str()),
                       Form("gMean_%s", baseName.c_str()),
                       nkbin, Be_bins, n_costheta_bin, costheta_bins);

    // 循环处理每个能量-角度单元
    const int nBinsX = h1->GetXaxis()->GetNbins();
    const int nBinsY = h1->GetYaxis()->GetNbins();
    
    for(int i = 1; i <= nBinsX; i++) {
        for(int j = 1; j <= nBinsY; j++) {
            TH1D *hProjection = h1->ProjectionZ(
                Form("hProj_%d_%d_%s", i, j, baseName.c_str()), i, i, j, j
            );

            const double energyValue = h1->GetXaxis()->GetBinCenter(i);
            const double costhetaValue = h1->GetYaxis()->GetBinCenter(j);
            
            // 跳过无效能量区域
            if(energyValue < 0.21 || energyValue > 20) {
                delete hProjection;  // 关键修复：防止内存泄漏
                gSigma3D->SetBinContent(i, j, 0);
                gMean3D->SetBinContent(i, j, 0);
                continue;
            }
            f->cd();
            hProjection->SetName(Form("hProj_%d_%d_%s", i, j, baseName.c_str()));
            hProjection->Rebin(2);
            hProjection->Write();
            // 高斯拟合参数设置
            TF1 *gaussFit = new TF1(
                "gaussFit", "gaus", 
                hProjection->GetXaxis()->GetXmin(),
                hProjection->GetXaxis()->GetXmax()
            );
            gaussFit->SetParameters(
                hProjection->GetMaximum(), 
                hProjection->GetMean(),
                hProjection->GetRMS()
            );

            if(hProjection->GetSumOfWeights() == 0){
                // cout<< baseName.c_str()<<"Warning: Empty histogram for bin (" 
                //     << i << ", " << j << "). Skipping fit." << endl;
                delete hProjection;  // 清理内存
                delete gaussFit;     // 清理拟合函数
                gSigma3D->SetBinContent(i, j, 0);
                gMean3D->SetBinContent(i, j, 0);
                continue;
            }
            // 迭代拟合直至收敛
            const double tolerance = 0.00002;
            const int maxIterations = 100;
            bool converged = false;
            double prevSigma = gaussFit->GetParameter(2);
            double prevMean = gaussFit->GetParameter(1);
            
            for(int iter = 0; iter < maxIterations && !converged; iter++) {
                hProjection->Fit(gaussFit, "Q0R");  // 静默模式
                const double currSigma = gaussFit->GetParameter(2);
                const double currMean = gaussFit->GetParameter(1);
                
                if(fabs(currSigma - prevSigma) < tolerance && 
                   fabs(currMean - prevMean) < tolerance) {
                    converged = true;
                    prevSigma = currSigma;
                    prevMean = currMean;
                    break;
                }
                
                prevSigma = currSigma;
                prevMean = currMean;
                gaussFit->SetRange(
                    currMean -  1.5* currSigma, 
                    currMean +  1.5* currSigma
                );
            }
            // if(converged) {
            //     std::cout<<"converged in 100"<<std::endl;
            // } else {
            //     std::cout<<"not converged in 100"<<std::endl;
            // }

            // 存储拟合结果
            gSigma3D->SetBinContent(i, j, prevSigma);
            gMean3D->SetBinContent(i, j, prevMean);
            gaussFit->SetName(
                Form("gaussFit_%d_%d_%s", i, j, baseName.c_str())
            );
            gaussFit->Write();

            // 清理内存
            delete hProjection;
            delete gaussFit;
        }
    }
}





void process_files() 
{
 
    string iso_type[] = {"he","li","be","b"};
    std::unordered_map<string, std::vector<int>> iso_map = {
        {"he", {3,4}},
        {"li", {6,7}},
        {"be", {7,9,10}},
        {"b",  {10,11}}
    };
    cout<<"Starting processing isotope FDA files..."<<endl;
    // string iso_type[] = {"he"};
    // std::unordered_map<std::string, std::vector<int>> iso_map = {
    //     {"he", {3,4}},
    // };

    TString iso_fda_dir="/eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda";

    // 创建输出文件


    for(string iso : iso_type) {
        // 获取当前同位素对应的质量数向量
        cout<<"Processing isotope: "<<iso<<endl;
        std::vector<int>& beam_isos = iso_map[iso];
            
        TFile* f4 = new TFile(Form("/eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda_norm/%s_norm2D.root", iso.c_str()), "recreate");
        //TFile* f5 = new TFile(Form("/eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda_norm/%s_proj.root", iso.c_str()), "recreate");


        // 遍历质量数，同时获取索引
        for(int idx = 0; idx < beam_isos.size(); idx++) {
            int beam_iso = beam_isos[idx];  // 当前质量数
            
            TFile* current_file = TFile::Open(iso_fda_dir + Form("/%s./%s%d_fda.root", 
                iso.c_str(), iso.c_str(), beam_iso));
            if(!current_file || current_file->IsZombie()) {
                std::cerr << "Error opening: " << iso_fda_dir + Form("/%s./%s%d_fda.root", 
                    iso.c_str(), iso.c_str(), beam_iso) << std::endl;
                continue;
            }

            // 遍历所有组合
            for(int i = 0; i < 2; i++) {
                for(int k = 0; k < 2; k++) {
                    for(int m = 1; m <= 3; m++) {
                        // 使用索引 idx 而不是 beam_iso
                        const TString hist_name = Form("hist04%d%d%d%d", 
                            idx, i, k, m);  // 使用 idx 作为第一个数字
                            
                        const std::string obj_id = Form(
                            "fda_%d_%d_%d_%d", 
                            idx, i, k, m   // 使用 idx 作为第一个数字
                        );

                        TH3D* h3 = dynamic_cast<TH3D*>(current_file->Get(hist_name));
                        if(!h3) {
                            std::cerr << "Missing histogram: " 
                                    << hist_name << std::endl;
                            continue;
                        }

                        // 执行拟合
                        TH2D *gSig = NULL;
                        TH2D *gMean = NULL;
                        Fit(h3, gSig, gMean, obj_id, f4);

                        // 保存结果
                        f4->cd();
                        if(gSig) gSig->Write();
                        if(gMean) gMean->Write();
                        
                        // 清理内存
                        if(gSig) delete gSig;
                        if(gMean) delete gMean;
                    }
                }
            }
            
            // 关闭当前文件
            current_file->Close();
            delete current_file;
        }
        f4->Close();
        delete f4;

    }


    
    cout<<"Finished processing isotope FDA files."<<endl;
}