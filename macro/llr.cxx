#include <cmath>
#include <iostream>
#include <unordered_map>
#include "TFile.h"
#include "TKey.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "sdst.h"
#include <string>
#include <vector>

Int_t nkbin = 73;
const double Be_bins[] = {
    0.08,   0.13,   0.17,   0.21,   0.27,   0.33,   0.41,   0.49,   0.59,
    0.70,   0.82,   0.96,   1.11,   1.28,   1.47,   1.68,   1.91,   2.16,
    2.44,   2.73,   3.06,   3.41,   3.79,   4.20,   4.65,   5.14,   5.64,
    6.18,   6.78,   7.42,   8.12,   8.86,   9.66,   10.51,  11.45,  12.45,
    13.50,  14.65,  15.84,  17.14,  18.54,  20.04,  21.64,  23.34,  25.19,
    27.13,  29.23,  31.48,  33.93,  36.53,  39.33,  42.33,  45.58,  49.08,
    53.08,  57.08,  61.58,  66.58,  72.57,  79.07,  86.57,  95.07,  104.57,
    115.57, 128.57, 144.57, 164.07, 188.57, 219.57, 261.57, 329.07, 439.07,
    649.07, 1649.07};



double k_s(TH1D*h1,TH1D*h2){
    double max_dis=-1;
    h1->Scale(1.0/h1->Integral()/h1->GetBinWidth(1));
    h2->Scale(1.0/h2->Integral()/h2->GetBinWidth(1));
    double cdf1=0;
    double cdf2=0;
    for(int i=1;i<=h1->GetNbinsX();i++){
        cdf1+=h1->GetBinContent(i)*h1->GetBinWidth(i);
        cdf2+=h2->GetBinContent(i)*h2->GetBinWidth(i);
        double dis=fabs(cdf1-cdf2);
        if(dis>max_dis) max_dis=dis;
    }
    return max_dis;

}


double overlap(TH1D*h1,TH1D*h2){
    double ov=0;
    h1->Scale(1.0/h1->Integral()/h1->GetBinWidth(1));
    h2->Scale(1.0/h2->Integral()/h2->GetBinWidth(1));
    for(int i=1;i<=h1->GetNbinsX();i++){
        double min_bin=fmin(h1->GetBinContent(i),h2->GetBinContent(i));
        ov+=min_bin*h1->GetBinWidth(i);
    }
    return ov;
}
void llr() {
    std::string iso_A[] = {"He", "Li", "Be", "B"};
    std::unordered_map<std::string, std::string> A_a = {
        {"He","he"}, {"Li","li"}, {"Be","be"}, {"B","b"}
    };
    std::unordered_map<std::string, std::vector<int>> iso_map = {
        {"He", {3,4}}, {"Li", {6,7}}, {"Be", {7,9,10}}, {"B", {10,11}}
    };

    int rebin = 20;
    int mass_rebin=2;  
    
    // 定义探测器范围和对应的索引
    struct DetectorRange {
        double min;
        double max;
        int index;
        const char* name;
    };
    
    DetectorRange ranges[] = {
        {0.5, 1.5, 0, "TOF"},
        {1.5, 3.0, 1, "NAF"}, 
        {3.0, 10.0, 2, "AGL"}
    };

    for (std::string A : iso_A) {
        std::string input_dir = Form("/eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_llr/%s./", A.c_str());
        
        // 打开文件
        std::vector<TFile*> flist;
        for (int mass : iso_map[A]) {
            std::string filename = input_dir + Form("%s%d_llr_temp.root", A_a[A].c_str(), mass);
            flist.push_back(TFile::Open(filename.c_str()));
        }

        // 获取直方图
        std::vector<std::vector<TH2D*>> h2list(flist.size(), std::vector<TH2D*>(3));
        std::vector<std::vector<TH2D*>> h2mass_list(flist.size(), std::vector<TH2D*>(3));
        
        for (int i = 0; i < flist.size(); i++) {
            for (int j = 0; j < 3; j++) {
                h2list[i][j] = (TH2D*)flist[i]->Get(Form("hist09900%d", j+1));
                h2mass_list[i][j] = (TH2D*)flist[i]->Get(Form("ty_hist0%d0", j+1));
                
                if (!h2list[i][j] || !h2mass_list[i][j]) {
                    std::cerr << "Error loading histograms for " << A << iso_map[A][i] << std::endl;
                }
            }
        }

        int nbins_x = h2list[0][0]->GetNbinsX();
        TCanvas *c1 = new TCanvas("c1", "c1", 1200, 600);
        c1->Print(Form("/eos/ams/user/s/selu/mdst/tianye/pdf/%s_llr_significance.pdf[", A.c_str()));
        gStyle->SetOptStat(0);
        // 为每个bin创建图形的函数
        auto createPlot = [&](int binx, int range_idx) {
            double bin_center = h2list[0][1]->GetXaxis()->GetBinCenter(binx);
            DetectorRange range = ranges[range_idx];
            
            if (bin_center <= range.min || bin_center >= range.max) return;
            
            c1->Clear();
            c1->Divide(2, 1);
            gStyle->SetErrorX(0);
            // 准备llr和mass直方图
            std::vector<TH1D*> h1_llr(flist.size());
            std::vector<TH1D*> h1_mass(flist.size());
            double max_llr = 0, max_mass = 0;
            
            for (int i = 0; i < flist.size(); i++) {
                // LLR直方图
                h1_llr[i] = h2list[i][range.index]->ProjectionY(
                    Form("h1_%s%d_llr_bin%.2f", A_a[A].c_str(), iso_map[A][i], bin_center), binx, binx);
                h1_llr[i]->Rebin(rebin);
                h1_llr[i]->Scale(1.0 / h1_llr[i]->Integral() / h1_llr[i]->GetBinWidth(1));
                max_llr = std::max(max_llr, h1_llr[i]->GetMaximum());
                
                // Mass直方图
                h1_mass[i] = h2mass_list[i][range.index]->ProjectionY(
                    Form("h1_%s%d_mass_bin%.2f", A_a[A].c_str(), iso_map[A][i], bin_center), binx, binx);
                h1_mass[i]->Rebin( mass_rebin);
                h1_mass[i]->Scale(1.0 / h1_mass[i]->Integral() / h1_mass[i]->GetBinWidth(1));
                max_mass = std::max(max_mass, h1_mass[i]->GetMaximum());
            }
            
            // 绘制LLR
            c1->cd(1);
            for (int i = 0; i < flist.size(); i++) {
                h1_llr[i]->SetLineColor(i+1);
                h1_llr[i]->SetTitle(Form("%s llr at %s=%.2f GeV/n", A.c_str(), range.name, bin_center));
                h1_llr[i]->GetXaxis()->SetTitle("llr");
                h1_llr[i]->GetYaxis()->SetTitle("Counts");
                h1_llr[i]->SetMaximum(1.2 * max_llr);
                h1_llr[i]->Draw(i == 0 ? "HIST" : "same HIST");
            }
            TLegend *legend = new TLegend(0.3, 0.7, 0.9, 0.9);
            if(A!="Be"){
                double k_val=k_s(h1_llr[0],h1_llr[1]);
                double ov_val=overlap(h1_llr[0],h1_llr[1]);
                legend->AddEntry(h1_llr[1], Form("%s%d, K-S=%.3f, Overlap=%.3f", A_a[A].c_str(), iso_map[A][1],k_val, ov_val), "l");
            }
            else{
                double k_val=k_s(h1_llr[1],h1_llr[2]);
                double ov_val=overlap(h1_llr[1],h1_llr[2]);
                legend->AddEntry(h1_llr[1], Form("%s%d, K-S=%.3f, Overlap=%.3f", A_a[A].c_str(), iso_map[A][1],k_val, ov_val), "l");
            }
            legend->SetTextSize(0.03);
            legend->Draw("same");
            
            // 绘制Mass
            c1->cd(2);
            for (int i = 0; i < flist.size(); i++) {
                h1_mass[i]->SetLineColor(i+1);
                h1_mass[i]->SetTitle(Form("%s mass at %s=%.2f GeV/n", A.c_str(), range.name, bin_center));
                h1_mass[i]->GetXaxis()->SetTitle("mass");
                h1_mass[i]->GetYaxis()->SetTitle("Counts");
                h1_mass[i]->SetMaximum(1.2 * max_mass);
                h1_mass[i]->Draw(i == 0 ? "HIST" : "same HIST");
            }

            TLegend *legend1 = new TLegend(0.3, 0.7, 0.9, 0.9);
            if(A!="Be"){
                double k_val=k_s(h1_mass[0],h1_mass[1]);
                double ov_val=overlap(h1_mass[0],h1_mass[1]);
                legend1->AddEntry(h1_mass[1], Form("%s%d, K-S=%.3f, Overlap=%.3f", A_a[A].c_str(), iso_map[A][1],k_val, ov_val), "l");
            }
            else{
                double k_val=k_s(h1_mass[1],h1_mass[2]);
                double ov_val=overlap(h1_mass[1],h1_mass[2]);
                legend1->AddEntry(h1_mass[1], Form("%s%d, K-S=%.3f, Overlap=%.3f", A_a[A].c_str(), iso_map[A][1],k_val, ov_val), "l");
            }
            legend1->SetTextSize(0.03);
            legend1->Draw("same");
            
            c1->Print(Form("/eos/ams/user/s/selu/mdst/tianye/pdf/%s_llr_significance.pdf", A.c_str()));
            
            // 清理内存
            for (auto h : h1_llr) delete h;
            for (auto h : h1_mass) delete h;
        };
        
        // 遍历所有bin和所有探测器范围
        for (int binx = 1; binx <= nbins_x; binx++) {
            for (int range_idx = 0; range_idx < 3; range_idx++) {
                createPlot(binx, range_idx);
            }
        }
        
        c1->Print(Form("/eos/ams/user/s/selu/mdst/tianye/pdf/%s_llr_significance.pdf]", A.c_str()));
        
        // 清理内存
        delete c1;
        for (auto file : flist) {
            if (file) file->Close();
        }
    }
}