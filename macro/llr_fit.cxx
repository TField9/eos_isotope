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

// RooFit 头文件
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooArgList.h"
#include "RooFit.h"

void addEpsilon(TH1* hist, double epsilon = 1e-12) {
    for (int i = 1; i <= hist->GetNbinsX(); i++) {
        double content = hist->GetBinContent(i);
        double error = hist->GetBinError(i);
        
        if (content <= 0) {
            hist->SetBinContent(i, epsilon);
        }
        
        // 确保误差不为0
        if (error <= 0) {
            hist->SetBinError(i, 1e10);
        }
    }
}

void refill(TH1D* hist){

    int nbin=hist->GetNbinsX();
    for(int i = 1; i <= nbin; i++) {
        double original_value = hist->GetBinContent(i);
        
        // 使用高斯随机数
        // 均值 = 原始值
        // 标准差 = sqrt(原始值) [假设方差与均值成正比，类似于泊松]
        double sigma = sqrt(fabs(hist->GetBinError(i)));  // 使用fabs避免负值问题
        
        double new_value = gRandom->Gaus(original_value, sigma);
        
        // 确保非负（直方图bin内容不能为负）
        if (new_value < 0) new_value = 0;
        
        // hist->SetBinContent(i, new_value);
        // hist->SetBinError(i, sigma);  // 误差设为标准差
        hist->SetBinContent(i, hist->GetBinContent(i));

    }


}

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

double k_s(TH1D* h1, TH1D* h2) {
    double max_dis = -1;
    h1->Scale(1.0 / h1->GetSumOfWeights() / h1->GetBinWidth(1));
    h2->Scale(1.0 / h2->GetSumOfWeights() / h2->GetBinWidth(1));
    double cdf1 = 0;
    double cdf2 = 0;
    for (int i = 1; i <= h1->GetNbinsX(); i++) {
        cdf1 += h1->GetBinContent(i) * h1->GetBinWidth(i);
        cdf2 += h2->GetBinContent(i) * h2->GetBinWidth(i);
        double dis = fabs(cdf1 - cdf2);
        if (dis > max_dis) max_dis = dis;
    }
    return max_dis;
}

double overlap(TH1D* h1, TH1D* h2) {
    double ov = 0;
    h1->Scale(1.0 / h1->GetSumOfWeights() / h1->GetBinWidth(1));
    h2->Scale(1.0 / h2->GetSumOfWeights() / h2->GetBinWidth(1));
    for (int i = 1; i <= h1->GetNbinsX(); i++) {
        double min_bin = fmin(h1->GetBinContent(i), h2->GetBinContent(i));
        ov += min_bin * h1->GetBinWidth(i);
    }
    return ov;
}
double calculateSignificance(std::vector<TH1D*> hlist, string iso, double eff,double &bincenter,double &sig_integral,double &bkg_integral) {

    std::unordered_map<std::string, std::vector<int>> iso_hist_idx = {
    {"He", {1,0}}, {"Li", {0,1}}, {"Be", {1,2}}, {"B", {0,1}}
    };
    double s_sum=0;
    double b_sum=0;
    double significance=0;
    TH1D* h_signal = hlist[iso_hist_idx[iso][1]];
    TH1D* h_background = hlist[iso_hist_idx[iso][0]];
    double total_signal = h_signal->Integral("width");
    //cout<<"Integral signal: "<<h_signal->Integral("width")<<", Integral background: "<<h_background->Integral("width")<<endl;
    if(iso!="He"){
        
        for(int i=1;i<=h_signal->GetNbinsX();i++){
            double s = h_signal->GetBinContent(i)*h_signal->GetBinWidth(i);
            double b = h_background->GetBinContent(i)*h_background->GetBinWidth(i);
            s_sum += s;
            b_sum += b;
            if(s_sum/total_signal>=eff){
                bincenter = h_signal->GetBinCenter(i);
                break;
            }
        }
        if(b_sum>0){
            sig_integral=s_sum;
            bkg_integral=b_sum;
            significance = s_sum / sqrt(b_sum);
        }
    }
    else{
        for(int i=h_signal->GetNbinsX();i>=1;i--){
            double s = h_signal->GetBinContent(i)*h_signal->GetBinWidth(i);
            double b = h_background->GetBinContent(i)*h_background->GetBinWidth(i);
            s_sum += s;
            b_sum += b;
            if(s_sum/total_signal>=eff){
                bincenter = h_signal->GetBinCenter(i);
                break;
            }
        }
        if(b_sum>0){
             sig_integral=s_sum;
            bkg_integral=b_sum;
            // significance = s_sum / sqrt(s_sum + b_sum);
            significance = s_sum / sqrt(b_sum);
        }
    }
    return significance;
}





// RooFit 拟合函数
void performRooFit(TString iso,TH1D* data_hist, std::vector<TH1D*> template_hists, 
                   std::vector<double> initial_ratios, 
                   double& chi2_ndf, std::vector<double>& fitted_ratios,
                   std::vector<double>& fitted_ratios_err,
                   const std::string& element_name, double kinetic_energy, 
                   const std::string& detector_name, int binx) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    // 创建RooFit变量
    RooRealVar x("x", "LLR", data_hist->GetXaxis()->GetXmin(), data_hist->GetXaxis()->GetXmax());
    
    // 创建数据直方图
    addEpsilon(data_hist);
    RooDataHist data("data", "Data", RooArgList(x), data_hist);
    
    
    // 创建模板PDFs
    std::vector<RooRealVar*> fractions;
    std::vector<RooHistPdf*> pdfs;
    RooArgList pdfList;
    RooArgList fracList;

    
    for (int i = 0; i < template_hists.size(); i++) {
        // 创建模板PDF
        addEpsilon(template_hists[i]);
        RooDataHist* template_data = new RooDataHist(Form("template_data_%d", i), 
                                                    Form("Template %d", i), 
                                                    RooArgList(x), 
                                                    template_hists[i]);
        RooHistPdf* pdf = new RooHistPdf(Form("pdf_%d", i), Form("PDF %d", i), x, *template_data);
        pdfs.push_back(pdf);
        pdfList.add(*pdf);
        
        // 创建分数参数（最后一个分数由1-其他分数之和决定）
        // if (i < template_hists.size() - 1) {
        //     RooRealVar* frac = new RooRealVar(Form("f%d", i), Form("Fraction %d", i), 
        //                                     0.5, 0.0, 1.0);
        //     fractions.push_back(frac);
        //     fracList.add(*frac);
        // }
                // 创建分数参数（最后一个分数由1-其他分数之和决定）
        if (i < template_hists.size() - 1) {
            // choose initial fraction based on provided initial_ratios if available,
            // otherwise fall back to defaults per isotope (iso)
            std::string iso_str = std::string(iso.Data());
            std::unordered_map<std::string, std::vector<double>> iso_defaults = {
            {"He", {0.15, 0.85}},
            {"Li", {0.5, 0.5}},
            {"Be", {0.6, 0.3, 0.1}},
            {"B",  {0.7, 0.3}}
            };

            double init = 0.5;

            auto it = iso_defaults.find(iso_str);
            if (it != iso_defaults.end() && i < (int)it->second.size()) {
                init = it->second[i];
            }
            
            // clamp to [0,1]
            if (init < 0.0) init = 0.0;
            if (init > 1.0) init = 1.0;

            RooRealVar* frac = new RooRealVar(Form("f%d", i), Form("Fraction %d", i),
                              init, 0.0, 1.0);
            fractions.push_back(frac);
            fracList.add(*frac);
        }
    }
    
    // 创建叠加PDF
    RooAddPdf model("model", "Total PDF", pdfList, fracList);
    
    // // 执行拟合
     RooChi2Var chi2("chi2", "Chi2", model, data);
     RooMinimizer minimizer(chi2);
    // minimizer.setStrategy(1);
    // minimizer.setPrintLevel(-1);
    // minimizer.minimize("Minuit", "migrad");
    // minimizer.hesse();
    // // 获取拟合结果
    // double chi2_value = chi2.getVal(); 

    // fitted_ratios.clear();

    // 执行拟合 - 改为NLL拟合
    // 改为扩展似然，考虑重加权
// 创建常数变量

    // RooNLLVar nll("nll", "Negative Log Likelihood", model, data);
    // RooMinimizer minimizer(nll);
// 替换这两行：
    minimizer.setEps(1e-8);  // 设置容差
    minimizer.setMaxIterations(100000);  // 设置最大迭代次数
    minimizer.setStrategy(0);
    minimizer.setPrintLevel(-1);
    minimizer.minimize("Minuit", "scan");

    // 阶段2：标准梯度下降
    minimizer.setStrategy(1);
    minimizer.minimize("Minuit", "migrad");

    // 阶段3：改进误差估计
    minimizer.minimize("Minuit", "improve");

    // 阶段4：高精度最终拟合
    minimizer.setStrategy(2);
        minimizer.setEps(1e-8);  // 设置容差
    minimizer.setMaxIterations(100000);  // 设置最大迭代次数
    minimizer.hesse();
    // 获取拟合结果
    double chi2_value = chi2.getVal(); 

    fitted_ratios.clear();






    
    // 获取拟合的比例
    double sum = 0;
    double sum_var = 0.0; // 用于累加方差

    for (size_t i = 0; i < fractions.size(); i++) {
        double val = fractions[i]->getVal();
        double err = fractions[i]->getError();
        
        fitted_ratios.push_back(val);
        fitted_ratios_err.push_back(err);
        
        sum += val;
        sum_var += err * err; // 方差累加
    }

    // 最后一个比例
    double last_ratio = 1.0 - sum;
    fitted_ratios.push_back(last_ratio);

    // 计算最后一个比例的误差（误差传递）
    double last_ratio_err = sqrt(sum_var);
    fitted_ratios_err.push_back(last_ratio_err);

    int nbins_nonzero = 0;
    for (int bin = 1; bin <= data_hist->GetNbinsX(); ++bin) {
        if (data_hist->GetBinContent(bin) > 1e-10) {
            nbins_nonzero++;
        }
    }
    chi2_ndf = chi2_value / (nbins_nonzero - fitted_ratios.size()-1);
    //chi2_ndf = 0;
    
    
   
}

void llr_fit() {
    std::string iso_A[] = {"He", "Li", "Be", "B"};
    //std::string iso_A[] = {"He"};
    //std::string iso_A[] = {"Be"};
    std::unordered_map<std::string, std::string> A_a = {
        {"He","he"}, {"Li","li"}, {"Be","be"}, {"B","b"}
    };
    std::unordered_map<std::string, std::vector<int>> iso_map = {
        {"He", {3,4}}, {"Li", {6,7}}, {"Be", {7,9,10}}, {"B", {10,11}}
    };

    std::unordered_map<std::string, std::vector<double>> iso_ratio = {
        {"He", {0.15,0.85}}, {"Li", {0.5,0.5}}, {"Be", {0.6,0.3,0.1}}, {"B", {0.7,0.3}}
    };


    
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


    // 创建输出文件保存拟合结果
    TFile* output_file = new TFile("/eos/ams/user/s/selu/mdst/tianye/isotope_fit_results.root", "RECREATE");

    for (std::string A : iso_A) {



    std::vector<double> eff_list = {0.9,0.8,0.7,0.6,0.5,0.4};
    std::vector<double> eff_bincenter_list(eff_list.size(), 0.0);
    std::vector<std::vector<double>> significance_llr_list(eff_list.size(), std::vector<double>(0));
    std::vector<std::vector<double>> significance_mass_list(eff_list.size(), std::vector<double>(0));
    std::vector<std::vector<double>> significance_integral_llr_sig_list(eff_list.size(), std::vector<double>(0));
    std::vector<std::vector<double>> significance_integral_llr_bkg_list(eff_list.size(), std::vector<double>(0));
    std::vector<std::vector<double>> significance_integral_mass_sig_list(eff_list.size(), std::vector<double>(0));
    std::vector<std::vector<double>> significance_integral_mass_bkg_list(eff_list.size(), std::vector<double>(0));

    std::vector<double> energy_list(0.0);
    std::vector<double> signal_ratio_llr_list(0.0);
    std::vector<double> signal_error_ratio_llr_list(0.0);
    std::vector<double> signal_ratio_mass_list(0.0);
    std::vector<double> signal_error_ratio_mass_list(0.0);


        std::string input_dir = Form("/eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_llr_reweight/%s./", A.c_str());
        
        // 打开文件
        std::vector<TFile*> flist_temp;
        std::vector<TFile*> flist_target;
        for (int mass : iso_map[A]) {
            std::string filename_temp = input_dir + Form("%s%d_llr_temp.root", A_a[A].c_str(), mass);
            std::string filename_target = input_dir + Form("%s%d_llr_target.root", A_a[A].c_str(), mass);
            flist_temp.push_back(TFile::Open(filename_temp.c_str()));
            flist_target.push_back(TFile::Open(filename_target.c_str()));
        }

        // 获取直方图
        std::vector<std::vector<TH2D*>> h2list_temp(flist_temp.size(), std::vector<TH2D*>(3));
        std::vector<std::vector<TH2D*>> h2list_target(flist_target.size(), std::vector<TH2D*>(3));
        std::vector<std::vector<TH2D*>> h2mass_list_temp(flist_temp.size(), std::vector<TH2D*>(3));
        std::vector<std::vector<TH2D*>> h2mass_list_target(flist_target.size(), std::vector<TH2D*>(3));
        
        for (int i = 0; i < flist_temp.size(); i++) {
            for (int j = 0; j < 3; j++) {
                h2list_temp[i][j] = (TH2D*)flist_temp[i]->Get(Form("hist09900%d", j+1));
                h2mass_list_temp[i][j] = (TH2D*)flist_temp[i]->Get(Form("ty_hist0%d0", j+1));
                h2list_target[i][j] = (TH2D*)flist_target[i]->Get(Form("hist09900%d", j+1));
                h2mass_list_target[i][j] = (TH2D*)flist_target[i]->Get(Form("ty_hist0%d0", j+1));
                
                if (!h2list_temp[i][j] || !h2mass_list_temp[i][j]) {
                    std::cerr << "Error loading histograms for " << A << iso_map[A][i] << std::endl;
                }
                if (!h2list_target[i][j] || !h2mass_list_target[i][j]) {
                    std::cerr << "Error loading histograms for " << A << iso_map[A][i] << std::endl;
                }
            }
        }

        int nbins_x = h2list_temp[0][0]->GetNbinsX();
        TCanvas* c1 = new TCanvas("c1", "c1", 1200, 600);
        c1->Print(Form("/eos/ams/user/s/selu/mdst/tianye/pdf/%s_llr.pdf[", A.c_str()));
        gStyle->SetOptStat(0);
        
        // 为每个bin创建图形的函数
        auto createPlot = [&](int binx, int range_idx) {
            double bin_center = h2list_temp[0][0]->GetXaxis()->GetBinCenter(binx);
            DetectorRange range = ranges[range_idx];
            
            if (bin_center <= range.min || bin_center >= range.max) return;
            energy_list.push_back(bin_center);
            c1->Clear();
            c1->Divide(2, 1);
            gStyle->SetErrorX(0);
            
            // 准备llr和mass直方图
            std::vector<TH1D*> h1_llr_temp(flist_temp.size());
            std::vector<TH1D*> h1_llr_temp_norebin(flist_temp.size());
            std::vector<TH1D*> h1_llr_target(flist_target.size());
            std::vector<TH1D*> h1_mass_temp(flist_temp.size());
            std::vector<TH1D*> h1_mass_temp_norebin(flist_temp.size());
            std::vector<TH1D*> h1_mass_target(flist_target.size());
            
            double max_llr = 0, max_mass = 0;
            int rebin, mass_rebin;
            for (int i = 0; i < flist_temp.size(); i++) {
                std::unordered_map<std::string, std::unordered_map<int,std::pair<int,int>>> rebin_map = {

                {"He", { {0,{200,80}}, {1,{200,125}}, {2,{200,100}} }},

                {"Li", { {0,{200,50}}, {1,{200,50}},  {2,{200,50}} }}, //perfect

                {"Be", { {0,{250,50}}, {1,{200,80}},  {2,{200,50}} }},

                {"B",  { {0,{200,50}}, {1,{200,50}},  {2,{200,40}} }}

                };
                rebin = rebin_map[A][range.index].first;
                mass_rebin = rebin_map[A][range.index].second;


                // LLR直方图 - temp
                h1_llr_temp[i] = h2list_temp[i][range.index]->ProjectionY(
                    Form("h1_%s%d_temp_llr_bin%.2f", A_a[A].c_str(), iso_map[A][i], bin_center), binx, binx);
                h1_llr_temp_norebin[i]= (TH1D*)h1_llr_temp[i]->Clone();
                //refill(h1_llr_temp[i]);
                //refill(h1_llr_temp_norebin[i]);
                h1_llr_temp[i]->Rebin(rebin);
                h1_llr_temp[i]->Scale(1.0 / h1_llr_temp[i]->GetSumOfWeights() / h1_llr_temp[i]->GetBinWidth(1));
                //cout<<"llr temp integral "<<A<<iso_map[A][i]<<" : "<<h1_llr_temp_norebin[i]->GetSumOfWeights()<<endl;
                h1_llr_temp_norebin[i]->Scale(1.0 / h1_llr_temp_norebin[i]->GetSumOfWeights() / h1_llr_temp_norebin[i]->GetBinWidth(1));
                //cout<<"llr temp integral "<<A<<iso_map[A][i]<<" : "<<h1_llr_temp_norebin[i]->Integral("width")<<endl;
                // LLR直方图 - target
                h1_llr_target[i] = h2list_target[i][range.index]->ProjectionY(
                    Form("h1_%s%d_target_llr_bin%.2f", A_a[A].c_str(), iso_map[A][i], bin_center), binx, binx);
                h1_llr_target[i]->Rebin(rebin);
                h1_llr_target[i]->Scale(1.0 / h1_llr_target[i]->GetSumOfWeights() / h1_llr_target[i]->GetBinWidth(1));
                
                max_llr = std::max(max_llr, std::max(h1_llr_temp[i]->GetMaximum(), h1_llr_target[i]->GetMaximum()));
                
                // Mass直方图
                h1_mass_temp[i] = h2mass_list_temp[i][range.index]->ProjectionY(
                    Form("h1_%s%d_temp_mass_bin%.2f", A_a[A].c_str(), iso_map[A][i], bin_center), binx, binx);
                h1_mass_temp_norebin[i]= (TH1D*)h1_mass_temp[i]->Clone();
                //refill(h1_mass_temp[i]);
                //refill(h1_mass_temp_norebin[i]);    
                h1_mass_temp[i]->Rebin(mass_rebin);
                //h1_mass_temp[i]->GetXaxis()->SetRangeUser(0,0.4);
                h1_mass_temp[i]->Scale(1.0 / h1_mass_temp[i]->GetSumOfWeights() / h1_mass_temp[i]->GetBinWidth(1));
                h1_mass_temp_norebin[i]->Scale(1.0 / h1_mass_temp_norebin[i]->GetSumOfWeights() / h1_mass_temp_norebin[i]->GetBinWidth(1));
  
                
                h1_mass_target[i] = h2mass_list_target[i][range.index]->ProjectionY(
                    Form("h1_%s%d_target_mass_bin%.2f", A_a[A].c_str(), iso_map[A][i], bin_center), binx, binx);
                h1_mass_target[i]->Rebin(mass_rebin);
                //h1_mass_target[i]->GetXaxis()->SetRangeUser(0,0.4);
                h1_mass_target[i]->Scale(1.0 / h1_mass_target[i]->GetSumOfWeights() / h1_mass_target[i]->GetBinWidth(1));
                
                max_mass = std::max(max_mass, std::max(h1_mass_temp[i]->GetMaximum(), h1_mass_target[i]->GetMaximum()));
            }

            TH1D* h1_llr_target_obs = (TH1D*)h1_llr_target[0]->Clone();  // 
            TH1D* h1_mass_target_obs = (TH1D*)h1_mass_target[0]->Clone();  //
            h1_llr_target_obs->Reset();
            h1_mass_target_obs->Reset();
            for(int i=0; i<flist_temp.size(); i++){
                TH1D* h1_llr_clone=(TH1D*)h1_llr_target[i]->Clone();
                TH1D* h1_mass_clone=(TH1D*)h1_mass_target[i]->Clone();
                h1_llr_clone->Scale(iso_ratio[A][i]);
                h1_mass_clone->Scale(iso_ratio[A][i]);
                h1_llr_target_obs->Add(h1_llr_clone);
                h1_mass_target_obs->Add(h1_mass_clone);
            }


            // 使用RooFit进行拟合
            double chi2_llr, chi2_mass;
            std::vector<double> fitted_ratios_llr, fitted_ratios_mass;
            std::vector<double> fitted_ratios_err_llr, fitted_ratios_err_mass;

            // 对LLR进行拟合
            performRooFit(A, h1_llr_target_obs, h1_llr_temp, iso_ratio[A], 
                        chi2_llr, fitted_ratios_llr, fitted_ratios_err_llr, A, bin_center, range.name, binx);
            // 对Mass进行拟合
            performRooFit(A, h1_mass_target_obs, h1_mass_temp, iso_ratio[A],
                        chi2_mass, fitted_ratios_mass, fitted_ratios_err_mass, A, bin_center, range.name, binx);

            // 绘制LLR
            c1->cd(1);

            // 创建堆叠直方图
            THStack *hs_llr = new THStack("hs_llr", Form("%s llr at %s=%.2f GeV/n", A.c_str(), range.name, bin_center));

            // 为每个同位素创建按拟合比例缩放的直方图
            std::vector<TH1D*> h1_llr_fitted(flist_temp.size());
            std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan}; // 颜色列表

            for (int i = 0; i < flist_temp.size(); i++) {
                h1_llr_fitted[i] = (TH1D*)h1_llr_temp[i]->Clone(Form("h1_%s%d_fitted_llr", A_a[A].c_str(), iso_map[A][i]));
                
                // 按照拟合比例缩放
                double scale_factor = fitted_ratios_llr[i] ;
                h1_llr_fitted[i]->Scale(scale_factor/h1_llr_fitted[i]->GetSumOfWeights()/h1_llr_fitted[i]->GetBinWidth(1));
                h1_llr_temp_norebin[i]->Scale(scale_factor/h1_llr_temp_norebin[i]->GetSumOfWeights()/h1_llr_temp_norebin[i]->GetBinWidth(1));
                
                // 设置颜色和填充样式
                h1_llr_fitted[i]->SetFillColorAlpha(colors[i], 0.7); // 0.5是透明度，范围0-1
                h1_llr_fitted[i]->SetLineColor(colors[i]);
                h1_llr_fitted[i]->SetFillStyle(1001);
                
                // 添加到堆叠直方图
                hs_llr->Add(h1_llr_fitted[i]);
            }
            
            for(int i=0;i<eff_list.size();i++){
                double llr_sig_integral=0;
                double llr_bkg_integral=0;
                significance_llr_list[i].push_back(calculateSignificance(h1_llr_temp_norebin,A,eff_list[i],eff_bincenter_list[i],llr_sig_integral,llr_bkg_integral));
                significance_integral_llr_sig_list[i].push_back(llr_sig_integral);
                significance_integral_llr_bkg_list[i].push_back(llr_bkg_integral);
            }
            if(A=="He"){
                signal_ratio_llr_list.push_back(fitted_ratios_llr[0]);
                signal_error_ratio_llr_list.push_back(fitted_ratios_err_llr[0]);
            }
            else if(A=="Be"){
                signal_ratio_llr_list.push_back(fitted_ratios_llr[2]);
                signal_error_ratio_llr_list.push_back(fitted_ratios_err_llr[2]);
            }
            else{
                signal_ratio_llr_list.push_back(fitted_ratios_llr[1]);
                signal_error_ratio_llr_list.push_back(fitted_ratios_err_llr[1]);
            }
            
            // 绘制堆叠直方图
            hs_llr->Draw("HIST");
            
            TH1D* hist_llr = (TH1D*)h1_llr_target_obs->Clone();
            hist_llr->Scale(1.0 / h1_llr_target_obs->GetSumOfWeights() / h1_llr_target_obs->GetBinWidth(1));
            int first_bin_llr = 1;
            int last_bin_llr = hist_llr->GetNbinsX();

            // 从前往后找第一个>0.1的bin
            while (first_bin_llr <= hist_llr->GetNbinsX() && hist_llr->GetBinContent(first_bin_llr) <= 0.1) {
                first_bin_llr++;
            }

            // 从后往前找最后一个>0.1的bin
            while (last_bin_llr >= 1 && hist_llr->GetBinContent(last_bin_llr) <= 0.1) {
                last_bin_llr--;
            }

            if (first_bin_llr <= last_bin_llr) {
                double x_min = hist_llr->GetXaxis()->GetBinLowEdge(first_bin_llr);
                double x_max = hist_llr->GetXaxis()->GetBinUpEdge(last_bin_llr);
                hs_llr->GetXaxis()->SetRangeUser(x_min, x_max);
            } else {
                hs_llr->GetXaxis()->SetRangeUser(0, 0.6);
            }
            // 或者使用当前Pad
            gPad->Modified();
            gPad->Update();
            hs_llr->GetXaxis()->SetTitle("llr");
            hs_llr->GetYaxis()->SetTitle("Counts");

            // 绘制观测数据点
            h1_llr_target_obs->SetLineColor(kBlack);
            h1_llr_target_obs->SetMarkerStyle(20);
            h1_llr_target_obs->SetMarkerSize(0.8);
            h1_llr_target_obs->Scale(1.0 / h1_llr_target_obs->GetSumOfWeights() / h1_llr_target_obs->GetBinWidth(1));
            h1_llr_target_obs->Draw("E SAME"); // 误差棒，与堆叠图叠加

            // 设置Y轴范围
            double y_max = std::max(hs_llr->GetMaximum(), h1_llr_target_obs->GetMaximum()) * 2;
            hs_llr->SetMaximum(y_max);

            // 创建图例
            TLegend* legend_llr = new TLegend(0.2, 0.6, 0.9, 0.9);
            legend_llr->SetHeader("LLR Fit Results", "C");

            // 添加观测数据
            legend_llr->AddEntry(h1_llr_target_obs, "Observed Data", "lep");

            // 添加各个拟合成分
            for (int i = 0; i < flist_temp.size(); i++) {
                TString label = Form("%s%d: %.3f #pm %.3f", A_a[A].c_str(), iso_map[A][i], fitted_ratios_llr[i], fitted_ratios_err_llr[i]);
                legend_llr->AddEntry(h1_llr_fitted[i], label, "f");
            }

            // 添加拟合质量信息
            legend_llr->AddEntry("", Form("#chi^{2}/ndf: %.2f", chi2_llr), "");
            for (int i = 0; i < eff_list.size(); i++) {
                legend_llr->AddEntry("", Form("Significance @ %.1f%% eff. (bin center=%.2f): %.3f", eff_list[i]*100, eff_bincenter_list[i], significance_llr_list[i].back()), "");
                legend_llr->AddEntry("", Form("Integral Signal: %.3f, Integral Background: %.3f", significance_integral_llr_sig_list[i].back(), significance_integral_llr_bkg_list[i].back()), "");
            }
            legend_llr->SetTextSize(0.02);
            legend_llr->Draw("same");

            // 对Mass也做同样的处理
            c1->cd(2);

            // 创建Mass的堆叠直方图
            THStack *hs_mass = new THStack("hs_mass", Form("%s mass at %s=%.2f GeV/n", A.c_str(), range.name, bin_center));

            // 为每个同位素创建按拟合比例缩放的直方图
            std::vector<TH1D*> h1_mass_fitted(flist_temp.size());

            for (int i = 0; i < flist_temp.size(); i++) {
                h1_mass_fitted[i] = (TH1D*)h1_mass_temp[i]->Clone(Form("h1_%s%d_fitted_mass", A_a[A].c_str(), iso_map[A][i]));
                
                // 按照拟合比例缩放
                double scale_factor = fitted_ratios_mass[i];
                h1_mass_fitted[i]->Scale(scale_factor/h1_mass_fitted[i]->GetSumOfWeights()/h1_mass_fitted[i]->GetBinWidth(1));
                h1_mass_temp_norebin[i]->Scale(scale_factor/h1_mass_temp_norebin[i]->GetSumOfWeights()/h1_mass_temp_norebin[i]->GetBinWidth(1));
                
                // 设置颜色和填充样式
                h1_mass_fitted[i]->SetFillColorAlpha(colors[i], 0.7); // 0.7是透明度，范围0-1
                h1_mass_fitted[i]->SetLineColor(colors[i]);
                h1_mass_fitted[i]->SetFillStyle(1001);
                
                // 添加到堆叠直方图
                hs_mass->Add(h1_mass_fitted[i]);
            }
            if(A=="He"){
                signal_ratio_mass_list.push_back(fitted_ratios_mass[0]);
                signal_error_ratio_mass_list.push_back(fitted_ratios_err_mass[0]);
            }
            else if(A=="Be"){
                signal_ratio_mass_list.push_back(fitted_ratios_mass[2]);
                signal_error_ratio_mass_list.push_back(fitted_ratios_err_mass[2]);
            }
            else{
                signal_ratio_mass_list.push_back(fitted_ratios_mass[1]);
                signal_error_ratio_mass_list.push_back(fitted_ratios_err_mass[1]);
            }

            for(int i=0;i<eff_list.size();i++){
                double mass_integral_sig=0;
                double mass_integral_bkg=0;
                significance_mass_list[i].push_back(calculateSignificance(h1_mass_temp_norebin,A,eff_list[i],eff_bincenter_list[i],mass_integral_sig,mass_integral_bkg));
                significance_integral_mass_sig_list[i].push_back(mass_integral_sig);
                significance_integral_mass_bkg_list[i].push_back(mass_integral_bkg);
            }

            // 绘制堆叠直方图
            
            hs_mass->Draw("HIST");
            TH1D* hist = (TH1D*)h1_mass_target_obs->Clone();
            hist->Scale(1.0 / hist->GetSumOfWeights() / hist->GetBinWidth(1));
            int first_bin = 1;
            int last_bin = hist->GetNbinsX();

            // 从前往后找第一个>0.1的bin
            while (first_bin <= hist->GetNbinsX() && hist->GetBinContent(first_bin) <= 0.1) {
                first_bin++;
            }

            // 从后往前找最后一个>0.1的bin
            while (last_bin >= 1 && hist->GetBinContent(last_bin) <= 0.1) {
                last_bin--;
            }

            if (first_bin <= last_bin) {
                double x_min = hist->GetXaxis()->GetBinLowEdge(first_bin);
                double x_max = hist->GetXaxis()->GetBinUpEdge(last_bin);
                hs_mass->GetXaxis()->SetRangeUser(x_min, x_max);
            } else {
                hs_mass->GetXaxis()->SetRangeUser(0, 0.6);
            }
            // 或者使用当前Pad
            gPad->Modified();
            gPad->Update();
            hs_mass->GetXaxis()->SetTitle("mass");
            hs_mass->GetYaxis()->SetTitle("Counts");

            // 绘制观测数据点
            h1_mass_target_obs->SetLineColor(kBlack);
            h1_mass_target_obs->SetMarkerStyle(20);
            h1_mass_target_obs->SetMarkerSize(0.8);
            h1_mass_target_obs->Scale(1.0 / h1_mass_target_obs->GetSumOfWeights() / h1_mass_target_obs->GetBinWidth(1));
            h1_mass_target_obs->Draw("E SAME");

            // 设置Y轴范围
            y_max = std::max(hs_mass->GetMaximum(), h1_mass_target_obs->GetMaximum()) * 2;
            hs_mass->SetMaximum(y_max);

            // 创建Mass的图例
            TLegend* legend_mass = new TLegend(0.3, 0.6, 0.9, 0.9);
            legend_mass->SetHeader("Mass Fit Results", "C");
            legend_mass->AddEntry(h1_mass_target_obs, "Observed Data", "lep");

            for (int i = 0; i < flist_temp.size(); i++) {
                TString label = Form("%s%d: %.3f #pm %.3f", A_a[A].c_str(), iso_map[A][i], fitted_ratios_mass[i], fitted_ratios_err_mass[i]);
                legend_mass->AddEntry(h1_mass_fitted[i], label, "f");
            }

            legend_mass->AddEntry("", Form("#chi^{2}/ndf: %.2f", chi2_mass), "");
            // for (int i = 0; i < eff_list.size(); i++) {
            //     legend_mass->AddEntry("", Form("Significance @ %.1f%% eff. (bin center=%.2f): %.3f", eff_list[i]*100, eff_bincenter_list[i], significance_mass_list[i].back()), "");
            // }

            for (int i = 0; i < eff_list.size(); i++) {
                legend_mass->AddEntry("", Form("Significance @ %.1f%% eff. (bin center=%.2f): %.3f. Nbinx: %d", eff_list[i]*100, eff_bincenter_list[i], significance_mass_list[i].back(),h1_mass_temp_norebin[0]->FindBin(eff_bincenter_list[i])), "");
                legend_mass->AddEntry("", Form("Integral Signal: %.3f, Integral Background: %.3f", significance_integral_mass_sig_list[i].back(), significance_integral_mass_bkg_list[i].back()), "");
            }
            legend_mass->SetTextSize(0.02);
            legend_mass->Draw("same");



            c1->Update();
            c1->Print(Form("/eos/ams/user/s/selu/mdst/tianye/pdf/%s_llr.pdf", A.c_str()));

            // 清理内存
            for (auto h : h1_llr_fitted) delete h;
            for (auto h : h1_mass_fitted) delete h;
            delete hs_llr;
            delete hs_mass;
            delete legend_llr;
            delete legend_mass;
            //delete info_box;
        };
        
        // 遍历所有bin和所有探测器范围
        for (int binx = 1; binx <= nbins_x; binx++) {
            for (int range_idx = 0; range_idx < 3; range_idx++) {
                createPlot(binx, range_idx);
            }
        }

        // 添加新的图：使用水平直线拟合LLR和Mass拟合比例的比值随能量的变化
if (energy_list.size() > 0) {
    TCanvas* c_ratio = new TCanvas("c_ratio", "LLR vs Mass Ratio Comparison with Constant Fit", 800, 600);
    
    // 计算比值和误差
    std::vector<double> ratio_values;
    std::vector<double> ratio_errors;
    std::vector<double> energy_values; // 用于拟合的能量值
    
    for (size_t i = 0; i < energy_list.size(); i++) {
        if (signal_ratio_mass_list[i] > 0) {
            // 计算比值
            double ratio = signal_ratio_llr_list[i] / signal_ratio_mass_list[i];
            ratio_values.push_back(ratio);
            
            // 计算比值的误差（误差传递公式）
            double relative_error_llr = signal_error_ratio_llr_list[i] / signal_ratio_llr_list[i];
            double relative_error_mass = signal_error_ratio_mass_list[i] / signal_ratio_mass_list[i];
            double relative_error_ratio = sqrt(relative_error_llr * relative_error_llr + 
                                             relative_error_mass * relative_error_mass);
            double ratio_error = ratio * relative_error_ratio;
            ratio_errors.push_back(ratio_error);
            
            energy_values.push_back(energy_list[i]);
        }
    }
    
    if (ratio_values.size() > 0) {
        TGraphErrors* gr_ratio = new TGraphErrors(ratio_values.size(), 
                                                 energy_values.data(), ratio_values.data(),
                                                 nullptr, ratio_errors.data());
        
        gr_ratio->SetTitle(Form("%s: LLR/Mass Signal Ratio vs Energy;Energy (GeV/n);LLR Ratio / Mass Ratio", A.c_str()));
        gr_ratio->SetMarkerStyle(20);
        gr_ratio->SetMarkerColor(kBlue);
        gr_ratio->SetLineColor(kBlue);
        gr_ratio->SetMarkerSize(1.2);
        
        // 设置坐标轴范围
        double y_min = 0.5, y_max = 1.5;
        for (size_t i = 0; i < ratio_values.size(); i++) {
            double val = ratio_values[i];
            double err = ratio_errors[i];
            y_min = std::min(y_min, val - 2*err);
            y_max = std::max(y_max, val + 2*err);
        }
        gr_ratio->GetYaxis()->SetRangeUser(y_min, y_max);
        
        gr_ratio->Draw("AP");
        
        // 进行水平直线拟合（常数函数）
        TF1* constant_fit = new TF1("constant_fit", "[0]", 
                                   energy_values.front(), energy_values.back());
        constant_fit->SetLineColor(kRed);
        constant_fit->SetLineWidth(2);
        constant_fit->SetParName(0, "Constant");
        
        // 设置初始参数值为数据的平均值
        double mean_ratio = TMath::Mean(ratio_values.begin(), ratio_values.end());
        constant_fit->SetParameter(0, mean_ratio);
        
        // 执行拟合
        TFitResultPtr fit_result = gr_ratio->Fit(constant_fit, "S"); // "S" 保存拟合结果
        
        // 添加参考线y=1
        TLine* line_one = new TLine(energy_values[0], 1.0, energy_values.back(), 1.0);
        line_one->SetLineColor(kGreen);
        line_one->SetLineStyle(2);
        line_one->SetLineWidth(2);
        line_one->Draw("same");
        
        // 添加图例
        TLegend* leg_ratio = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg_ratio->AddEntry(gr_ratio, "LLR/Mass Ratio", "lep");
        leg_ratio->AddEntry(constant_fit, "Constant Fit", "l");
        leg_ratio->AddEntry(line_one, "Reference (y=1)", "l");
        leg_ratio->Draw();
        
        // 添加拟合结果信息
        TLatex* fit_info = new TLatex();
        fit_info->SetNDC();
        fit_info->SetTextSize(0.03);
        
        double chi2 = fit_result->Chi2();
        int ndf = fit_result->Ndf();
        double constant = constant_fit->GetParameter(0);
        double constant_err = constant_fit->GetParError(0);
        
        fit_info->DrawLatex(0.15, 0.85, Form("Fit: y = constant"));
        fit_info->DrawLatex(0.15, 0.80, Form("Constant = %.3f #pm %.3f", constant, constant_err));
        fit_info->DrawLatex(0.15, 0.75, Form("#chi^{2}/ndf = %.2f/%d = %.2f", chi2, ndf, chi2/ndf));
        
        // 计算与1的差异显著性
        double diff_from_one = fabs(constant - 1.0);
        double significance = (constant_err > 0) ? diff_from_one / constant_err : 0;
        fit_info->DrawLatex(0.15, 0.70, Form("Difference from 1: %.3f (#sigma = %.2f)", 
                                           constant - 1.0, significance));
        
        // 计算加权平均值（作为比较）
        double weighted_mean = 0;
        double sum_weights = 0;
        for (size_t i = 0; i < ratio_values.size(); i++) {
            double weight = 1.0 / (ratio_errors[i] * ratio_errors[i]);
            weighted_mean += ratio_values[i] * weight;
            sum_weights += weight;
        }
        if (sum_weights > 0) weighted_mean /= sum_weights;
        fit_info->DrawLatex(0.15, 0.65, Form("Weighted mean: %.3f", weighted_mean));
        
        c_ratio->Update();
        c_ratio->Print(Form("/eos/ams/user/s/selu/mdst/tianye/pdf/%s_llr.pdf", A.c_str()));
        
        // 保存到输出文件
        output_file->cd();
        gr_ratio->Write(Form("%s_llr_mass_ratio", A.c_str()));
        constant_fit->Write(Form("%s_constant_fit", A.c_str()));
        
        delete gr_ratio;
        delete constant_fit;
        delete line_one;
        delete leg_ratio;
        delete fit_info;
    }

    //draw llr ，mass original signal fraction

    {

        TGraphErrors* llr_signal_fraction = new TGraphErrors(ratio_values.size(), 
                                                 energy_values.data(), signal_ratio_llr_list.data(),
                                                 nullptr, signal_error_ratio_mass_list.data());
        llr_signal_fraction->SetTitle(Form("%s: LLR Signal Fraction vs Energy;Energy (GeV/n);LLR Signal Fraction", A.c_str()));
        llr_signal_fraction->SetMarkerStyle(20);
        llr_signal_fraction->SetMarkerColor(kBlue);
        llr_signal_fraction->SetLineColor(kBlue);
        llr_signal_fraction->SetMarkerSize(1.2);
        TGraphErrors* mass_signal_fraction = new TGraphErrors(ratio_values.size(), 
                                                 energy_values.data(), signal_ratio_mass_list.data(),
                                                 nullptr, signal_error_ratio_mass_list.data());
        mass_signal_fraction->SetTitle(Form("%s: Mass Signal Fraction vs Energy;Energy (GeV/n);Mass Signal Fraction", A.c_str()));
        mass_signal_fraction->SetMarkerStyle(20);
        mass_signal_fraction->SetMarkerColor(kRed);
        mass_signal_fraction->SetLineColor(kRed);
        mass_signal_fraction->SetMarkerSize(1.2);

        double y_fraction=0;
        if(A=="He"){
            y_fraction=0.15;
        }
        else if(A=="Li"){
            y_fraction=0.5;
        }
        else if(A=="Be"){
            y_fraction=0.1;
        }
        else{
            y_fraction=0.3;
        }
        double y_ratio=0.5;
        mass_signal_fraction->GetYaxis()->SetRangeUser((1-y_ratio)*y_fraction, (1+y_ratio)*y_fraction);
        llr_signal_fraction->GetYaxis()->SetRangeUser((1-y_ratio)*y_fraction, (1+y_ratio)*y_fraction);

        TF1* constant_fit_fraction_llr = new TF1("constant_fit_fraction_llr", "[0]", 
                                   energy_values.front(), energy_values.back());
        constant_fit_fraction_llr->SetLineColor(kBlue);
        constant_fit_fraction_llr->SetLineWidth(2);
        constant_fit_fraction_llr->SetParName(0, "Constant");
        constant_fit_fraction_llr->SetParameter(0, y_fraction);
        llr_signal_fraction->Fit(constant_fit_fraction_llr, "S");

        TF1* constant_fit_fraction_mass = new TF1("constant_fit_fraction_mass", "[0]", 
                                   energy_values.front(), energy_values.back());
        constant_fit_fraction_mass->SetLineColor(kRed);
        constant_fit_fraction_mass->SetLineWidth(2);
        constant_fit_fraction_mass->SetParName(0, "Constant");
        constant_fit_fraction_mass->SetParameter(0, y_fraction);
        mass_signal_fraction->Fit(constant_fit_fraction_mass, "S");

        TLine * line_fraction = new TLine(energy_values[0], y_fraction, energy_values.back(), y_fraction);
        line_fraction->SetLineColor(kGreen);
        line_fraction->SetLineStyle(2);
        line_fraction->SetLineWidth(2);

        c_ratio->Clear();
        llr_signal_fraction->Draw("AP");
        mass_signal_fraction->Draw("P SAME");
        TLegend* leg_signal = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg_signal->AddEntry(llr_signal_fraction, "LLR Signal Fraction", "lep");
        leg_signal->AddEntry(mass_signal_fraction, "Mass Signal Fraction", "lep");
        leg_signal->Draw();
        line_fraction->Draw("same");
        TLegend* leg_fraction = new TLegend(0.15, 0.75, 0.5, 0.9);
        leg_fraction->AddEntry(constant_fit_fraction_llr, Form("LLR Fit: %.3f #pm %.3f", constant_fit_fraction_llr->GetParameter(0), constant_fit_fraction_llr->GetParError(0)), "l");
        leg_fraction->AddEntry(constant_fit_fraction_mass, Form("Mass Fit: %.3f #pm %.3f", constant_fit_fraction_mass->GetParameter(0), constant_fit_fraction_mass->GetParError(0)), "l");
        leg_fraction->AddEntry(line_fraction, Form("Expected: %.3f", y_fraction), "l");
        leg_fraction->Draw();
        c_ratio->Update();
        c_ratio->Print(Form("/eos/ams/user/s/selu/mdst/tianye/pdf/%s_llr.pdf", A.c_str()));
        








    }



    //draw significance vs energy

    {
        

        for(int at_eff=0; at_eff<eff_list.size(); at_eff++){
        TGraphErrors* significance_llr_graph = new TGraphErrors(energy_values.size(), 
                                                 energy_values.data(), significance_llr_list[at_eff].data(),
                                                 nullptr, nullptr);

        TGraphErrors* significance_mass_graph = new TGraphErrors(energy_values.size(), 
                                                 energy_values.data(), significance_mass_list[at_eff].data(),
                                                 nullptr, nullptr);
        significance_llr_graph->SetTitle(Form("%s: Significance vs Energy(@eff %.0f %%);Energy (GeV/n);Significance ", A.c_str(), (eff_list[at_eff]*100)));
        significance_llr_graph->SetMarkerStyle(20);
        significance_llr_graph->SetMarkerColor(kBlue);
        significance_llr_graph->SetLineColor(kBlue);
        significance_llr_graph->SetMarkerSize(1.2);
        significance_mass_graph->SetTitle(Form("%s: Significance vs Energy(@eff %.0f %%);Energy (GeV/n);Significance ", A.c_str(), (eff_list[at_eff]*100)));
        significance_mass_graph->SetMarkerStyle(20);
        significance_mass_graph->SetMarkerColor(kRed);
        significance_mass_graph->SetLineColor(kRed);
        significance_mass_graph->SetMarkerSize(1.2);

        {
            // Determine global min/max from both graphs
            double y_min = 1e9;
            double y_max = -1e9;
            int n1 = significance_llr_graph->GetN();
            for (int i = 0; i < n1; ++i) {
                double x, y;
                significance_llr_graph->GetPoint(i, x, y);
                if (y < y_min) y_min = y;
                if (y > y_max) y_max = y;
            }
            int n2 = significance_mass_graph->GetN();
            for (int i = 0; i < n2; ++i) {
                double x, y;
                significance_mass_graph->GetPoint(i, x, y);
                if (y < y_min) y_min = y;
                if (y > y_max) y_max = y;
            }

            if (y_min > y_max) { // no points
                y_min = 0.0;
                y_max = 1.0;
            }

            double y_low = y_min * 0.5;
            double y_high = y_max * 1.5;
            if (y_low == y_high) { // protect against degenerate range
                y_low = y_min - 0.5 * (fabs(y_min) + 1e-6);
                y_high = y_max + 0.5 * (fabs(y_max) + 1e-6);
            }

            significance_llr_graph->GetYaxis()->SetRangeUser(y_low, y_high);
            significance_mass_graph->GetYaxis()->SetRangeUser(y_low, y_high);

            // Draw and annotate
            significance_llr_graph->Draw("AP");
            significance_mass_graph->Draw("P SAME");
            TLegend* leg_sign = new TLegend(0.7, 0.7, 0.9, 0.9);
            leg_sign->AddEntry(significance_llr_graph, "LLR Significance", "lep");
            leg_sign->AddEntry(significance_mass_graph, "Mass Significance", "lep");
            leg_sign->Draw();

            c_ratio->Update();
            c_ratio->Print(Form("/eos/ams/user/s/selu/mdst/tianye/pdf/%s_llr.pdf", A.c_str()));

            // save graphs
            output_file->cd();
            significance_llr_graph->Write(Form("%s_significance_llr", A.c_str()));
            significance_mass_graph->Write(Form("%s_significance_mass", A.c_str()));

            delete significance_llr_graph;
            delete significance_mass_graph;
            delete leg_sign;
        }
    }

    
    
    
    
    
    
    
    
    
    }
    
    delete c_ratio;
}






        
        c1->Print(Form("/eos/ams/user/s/selu/mdst/tianye/pdf/%s_llr.pdf]", A.c_str()));
        
        // 清理内存
        delete c1;
        
        for (auto file : flist_temp) {
            if (file) file->Close();
        }
        for (auto file : flist_target) {
            if (file) file->Close();
        }
    }
    
    // 关闭输出文件
    output_file->Close();
    delete output_file;
}