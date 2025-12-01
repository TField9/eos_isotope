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

#include <cstdint>

int count_hit(uint64_t value){
#if defined(__GNUG__) || defined(__clang__)
    return __builtin_popcountll(value);
#else
    int cnt = 0;
    while(value){
        cnt += value & 1ULL;
        value >>= 1;
    }
    return cnt;
#endif
}

const double Be_bins[] = {
    0.27, 0.33, 0.41, 0.49, 0.59,
    0.70, 0.82, 0.96, 1.11, 1.28, 1.47, 1.68, 1.91, 2.16,
    2.44, 2.73, 3.06, 3.41, 3.79, 4.20, 4.65, 5.14, 5.64,
    6.18, 6.78, 7.42, 8.12, 8.86, 9.66, 10.51
};

void normalize(){

    TFile *normalize_para = new TFile("/eos/ams/user/s/selu/mdst/tianye/root/normalize_parameter/normalize_para.root","recreate");
    TChain *c_be7 = new TChain("tree");
    TChain *c_be9 = new TChain("tree");
    TChain *c_be10 = new TChain("tree");

    for(int i=0;i<100;i++){
        c_be7->Add(Form("/eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_stdcut/sdst_be7_%d_%d_stdcut.root", i*10+1, i*10+10));
        c_be9->Add(Form("/eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_stdcut/sdst_be9_%d_%d_stdcut.root", i*10+1, i*10+10));
        c_be10->Add(Form("/eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_stdcut/sdst_be10_%d_%d_stdcut.root", i*10+1, i*10+10));

    }

    DST *dst_be7 = new DST();
    DST *dst_be9 = new DST();
    DST *dst_be10 = new DST();

    dst_be7->SetAddress(c_be7);
    dst_be9->SetAddress(c_be9);
    dst_be10->SetAddress(c_be10);
    cout<<"Total entries: "<<c_be7->GetEntries()<<endl;
    TH2D *h_be7[3][7];
    TH2D *h_be9[3][7];
    TH2D *h_be10[3][7];

    int Be_bin_size = sizeof(Be_bins)/sizeof(Be_bins[0]) -1;

    for(int i=1;i<8;i++){
        for(int j=0;j<3;j++){
            float bin_edge;
            if(j==0) bin_edge = 0.15;
            else bin_edge = 0.08;

            h_be7[j][i-1]= new TH2D(Form("be7_%d_%d",j,i),Form("be7_%d_%d",j,i),Be_bin_size,Be_bins,4000,-bin_edge,bin_edge);
            h_be9[j][i-1]= new TH2D(Form("be9_%d_%d",j,i),Form("be9_%d_%d",j,i),Be_bin_size,Be_bins,4000,-bin_edge,bin_edge);
            h_be10[j][i-1]= new TH2D(Form("be10_%d_%d",j,i),Form("be10_%d_%d",j,i),Be_bin_size,Be_bins,4000,-bin_edge,bin_edge);
        }
    }
    const Int_t be7_nEntries = c_be7->GetEntries();
    for (Int_t i = 0; i < be7_nEntries; i++) {
        c_be7->GetEntry(i);
        if(dst_be7->fEk.Ek[0]<1.5&&dst_be7->fEk.Ek[0]>0.25){
            for(int j=1;j<8;j++){
                h_be7[0][j-1]->Fill(dst_be7->fEk.Ek[0], dst_be7->fTrack.cooy[j]-dst_be7->fTrack.yglob_JN[j]);
            }
        
        }
        if(dst_be7->fRichcut_flag.richcut_flag[0]==1&&dst_be7->fRichcut_flag.richcut_flag[1]==0&&dst_be7->fEk.Ek[1]<3.5&&
           dst_be7->fEk.Ek[1]>1.3){
            for(int j=1;j<8;j++){
                h_be7[1][j-1]->Fill(dst_be7->fEk.Ek[1], dst_be7->fTrack.cooy[j]-dst_be7->fTrack.yglob_JN[j]);
            }
        }
        if(dst_be7->fRichcut_flag.richcut_flag[0]==1&&dst_be7->fRichcut_flag.richcut_flag[1]==1&&
           dst_be7->fEk.Ek[1]<10&& dst_be7->fEk.Ek[1]>3.1){
            for(int j=1;j<8;j++){
                h_be7[2][j-1]->Fill(dst_be7->fEk.Ek[1], dst_be7->fTrack.cooy[j]-dst_be7->fTrack.yglob_JN[j]);
            }
        }

    }

    // process be9 entries
    cout<<"processing be9"<<endl;
    const Int_t be9_nEntries = c_be9->GetEntries();
    for (Int_t idx = 0; idx < be9_nEntries; ++idx) {
        c_be9->GetEntry(idx);
        if (dst_be9->fEk.Ek[0] > 0.25 && dst_be9->fEk.Ek[0] < 1.5) {
            for (int j = 1; j < 8; ++j) {
                h_be9[0][j-1]->Fill(dst_be9->fEk.Ek[0], dst_be9->fTrack.cooy[j] - dst_be9->fTrack.yglob_JN[j]);
            }
        }
        if (dst_be9->fRichcut_flag.richcut_flag[0] == 1 && dst_be9->fRichcut_flag.richcut_flag[1] == 0 &&
            dst_be9->fEk.Ek[1] > 1.3 && dst_be9->fEk.Ek[1] < 3.5) {
            for (int j = 1; j < 8; ++j) {
                h_be9[1][j-1]->Fill(dst_be9->fEk.Ek[1], dst_be9->fTrack.cooy[j] - dst_be9->fTrack.yglob_JN[j]);
            }
        }
        if (dst_be9->fRichcut_flag.richcut_flag[0] == 1 && dst_be9->fRichcut_flag.richcut_flag[1] == 1 &&
            dst_be9->fEk.Ek[1] > 3.1 && dst_be9->fEk.Ek[1] < 10.0) {
            for (int j = 1; j < 8; ++j) {
                h_be9[2][j-1]->Fill(dst_be9->fEk.Ek[1], dst_be9->fTrack.cooy[j] - dst_be9->fTrack.yglob_JN[j]);
            }
        }
    }

    // process be10 entries
    cout<<"processing be10"<<endl;
    const Int_t be10_nEntries = c_be10->GetEntries();
    for (Int_t idx = 0; idx < be10_nEntries; ++idx) {
        c_be10->GetEntry(idx);
        if (dst_be10->fEk.Ek[0] > 0.25 && dst_be10->fEk.Ek[0] < 1.5) {
            for (int j = 1; j < 8; ++j) {
                h_be10[0][j-1]->Fill(dst_be10->fEk.Ek[0], dst_be10->fTrack.cooy[j] - dst_be10->fTrack.yglob_JN[j]);
            }
        }
        if (dst_be10->fRichcut_flag.richcut_flag[0] == 1 && dst_be10->fRichcut_flag.richcut_flag[1] == 0 &&
            dst_be10->fEk.Ek[1] > 1.3 && dst_be10->fEk.Ek[1] < 3.5) {
            for (int j = 1; j < 8; ++j) {
                h_be10[1][j-1]->Fill(dst_be10->fEk.Ek[1], dst_be10->fTrack.cooy[j] - dst_be10->fTrack.yglob_JN[j]);
            }
        }
        if (dst_be10->fRichcut_flag.richcut_flag[0] == 1 && dst_be10->fRichcut_flag.richcut_flag[1] == 1 &&
            dst_be10->fEk.Ek[1] > 3.1 && dst_be10->fEk.Ek[1] < 10.0) {
            for (int j = 1; j < 8; ++j) {
                h_be10[2][j-1]->Fill(dst_be10->fEk.Ek[1], dst_be10->fTrack.cooy[j] - dst_be10->fTrack.yglob_JN[j]);
            }
        }
    }


    TH1D *h_mean_be7[3][7];
    TH1D *h_mean_be9[3][7];
    TH1D *h_mean_be10[3][7];

    TH1D *h_sigma_be7[3][7];
    TH1D *h_sigma_be9[3][7];
    TH1D *h_sigma_be10[3][7];

    for(int i=1;i<8;i++){
        for(int j=0;j<3;j++){
            h_mean_be7[j][i-1]= new TH1D(Form("mean_be7_%d_%d",j,i),Form("mean_be7_detec%d_layer%d",j,i),Be_bin_size,Be_bins);
            h_mean_be9[j][i-1]= new TH1D(Form("mean_be9_%d_%d",j,i),Form("mean_be9_detec%d_layer%d",j,i),Be_bin_size,Be_bins);
            h_mean_be10[j][i-1]= new TH1D(Form("mean_be10_%d_%d",j,i),Form("mean_be10_detec%d_layer%d",j,i),Be_bin_size,Be_bins);

            h_sigma_be7[j][i-1]= new TH1D(Form("sigma_be7_%d_%d",j,i),Form("sigma_be7_detec%d_layer%d",j,i),Be_bin_size,Be_bins);
            h_sigma_be9[j][i-1]= new TH1D(Form("sigma_be9_%d_%d",j,i),Form("sigma_be9_detec%d_layer%d",j,i),Be_bin_size,Be_bins);
            h_sigma_be10[j][i-1]= new TH1D(Form("sigma_be10_%d_%d",j,i),Form("sigma_be10_detec%d_layer%d",j,i),Be_bin_size,Be_bins);
        }
    }

    // Fit and extract mean and sigma for be7, be9, be10
    auto fit_and_fill = [&](TH2D* h2, TH1D* hmean, TH1D* hsigma, const char* prefix) {
        for (int bin = 1; bin <= Be_bin_size; ++bin) {
            // projection for this energy bin
            TH1D *py = h2->ProjectionY(Form("%s_proj_%d", prefix, bin), bin, bin);
            double mean = 0.0, sigma = 0.0;
            double mean_err = 0.0, sigma_err = 0.0;

            if (py->GetEntries() > 10) {
                // broad fit (quiet)
                py->Fit("gaus", "RLQ");
                TF1 *f = py->GetFunction("gaus");
                if (f) {
                    mean = f->GetParameter(1);
                    sigma = f->GetParameter(2);
                    if (std::isfinite(sigma) && sigma > 0) {
                        double xmin = mean - 1.0 * sigma;
                        double xmax = mean + 1.0 * sigma;
                        // refined fit around ±2 sigma
                        py->Fit("gaus", "RQL", "", xmin, xmax);
                        f = py->GetFunction("gaus");
                        if (f) {
                            mean = f->GetParameter(1);
                            sigma = f->GetParameter(2);
                            mean_err = f->GetParError(1);
                            sigma_err = f->GetParError(2);
                        }
                    }
                }
            } else {
                // fallback to moments if too few entries
                mean = -9.0;
                sigma = -9.0;
            }

            // fill the mean and sigma histograms (bin index matches energy bin)
            hmean->SetBinContent(bin, mean);
            // try to set a reasonable error for the mean: use mean error from projection if available
            hmean->SetBinError(bin, (mean_err > 0.0) ? mean_err : 0.0);
            hsigma->SetBinContent(bin, sigma);
            // no simple built-in error for sigma from projection; set to 0 or estimate if needed
            hsigma->SetBinError(bin, (sigma_err > 0.0) ? sigma_err : 0.0);

            delete py;
        }
    };

    // Loop over layers (i) and categories (j) to process each 2D histogram
    for (int i = 1; i < 8; ++i) {
        for (int j = 0; j < 3; ++j) {
            fit_and_fill(h_be7[j][i-1], h_mean_be7[j][i-1], h_sigma_be7[j][i-1],
                         Form("be7_%d_%d", j, i));
            fit_and_fill(h_be9[j][i-1], h_mean_be9[j][i-1], h_sigma_be9[j][i-1],
                         Form("be9_%d_%d", j, i));
            fit_and_fill(h_be10[j][i-1], h_mean_be10[j][i-1], h_sigma_be10[j][i-1],
                         Form("be10_%d_%d", j, i));
        }
    }


    
    TCanvas *c1= new TCanvas("c1","c1",800,600);
    c1->Print("/eos/ams/user/s/selu/mdst/tianye/transformer/pdf/normalize.pdf[");
    for(int i=0;i<7;i++){
        for(int j=0;j<3;j++){
            c1->Clear();
            c1->Divide(2,1);
            c1->cd(1);
            h_mean_be10[j][i]->SetLineColor(kRed);
            h_mean_be10[j][i]->SetMarkerColor(kRed);
            h_mean_be10[j][i]->SetMarkerStyle(20);
            h_mean_be10[j][i]->GetYaxis()->SetRangeUser(-0.1,0.1);
            //h_mean_be10[j][i]->SetTitle(Form("Mean Y residuals Be10 Layer %d Category %d", i+1, j));
            h_mean_be10[j][i]->GetXaxis()->SetTitle("Kinetic Energy (GeV)");
            h_mean_be10[j][i]->GetYaxis()->SetTitle("Mean Y Residuals (cm)");
            h_mean_be10[j][i]->Draw("E");
            c1->cd(2);
            h_sigma_be10[j][i]->SetLineColor(kBlue);
            h_sigma_be10[j][i]->SetMarkerColor(kBlue);
            h_sigma_be10[j][i]->SetMarkerStyle(20);
            h_sigma_be10[j][i]->GetYaxis()->SetRangeUser(0,0.015);
            //h_sigma_be10[j][i]->SetTitle(Form("Sigma Y residuals Be10 Layer %d Category %d", i+1, j));
            h_sigma_be10[j][i]->GetXaxis()->SetTitle("Kinetic Energy (GeV)");
            h_sigma_be10[j][i]->GetYaxis()->SetTitle("Sigma Y Residuals (cm)");
            h_sigma_be10[j][i]->Draw("E");
            c1->Print("/eos/ams/user/s/selu/mdst/tianye/transformer/pdf/normalize.pdf");

        }
    }

    c1->Print("/eos/ams/user/s/selu/mdst/tianye/transformer/pdf/normalize.pdf]");

    TCanvas *c2= new TCanvas("c2","c2",800,600);
    c2->Print("/eos/ams/user/s/selu/mdst/tianye/transformer/pdf/normalize_be7.pdf[");
    //gStyle->SetOptStat(0);
        for(int bin=1;bin<=Be_bin_size;bin++){ 
                for(int layer=6;layer<8;layer++){
                    TH1D *proj=h_be10[0][layer-1]->ProjectionY(Form("proj_mean_be10_layer%d_bin%d",layer,bin),bin,bin);
                    if (!proj) continue;

                    double mean = -9.0, sigma = -9.0, mean_err = 0.0, sigma_err = 0.0;

                    if (proj->GetEntries() > 10) {
                        // first broad fit (quiet)
                        proj->Fit("gaus", "RQL");
                        TF1 *f = proj->GetFunction("gaus");
                        if (f) {
                            mean = f->GetParameter(1);
                            sigma = f->GetParameter(2);
                            if (std::isfinite(sigma) && sigma > 0.0) {
                                // refined fit within ±1 sigma
                                double xmin = mean - 1.0 * sigma;
                                double xmax = mean + 1.0 * sigma;
                                proj->Fit("gaus", "RLQ", "", xmin, xmax);
                                f = proj->GetFunction("gaus");
                                if (f) {
                                    mean = f->GetParameter(1);
                                    sigma = f->GetParameter(2);
                                    mean_err = f->GetParError(1);
                                    sigma_err = f->GetParError(2);
                                }
                            }
                        }
                    } else {
                        // too few entries: mark invalid
                        mean = -9.0;
                        sigma = -9.0;
                    }

                    // optional: draw and save the projection with fit
                    c2->Clear();
                    proj->SetMarkerStyle(20);
                    proj->SetMarkerColor(kBlack);
                    proj->Draw("E");
                    TF1 *finalFit = proj->GetFunction("gaus");
                    if (finalFit) finalFit->SetLineColor(kRed);
                    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
                    legend->AddEntry(proj, "Y Residuals", "lep");
                    if (finalFit) legend->AddEntry(finalFit, "Gaussian Fit", "l");

                    if (mean == -9.0 && sigma == -9.0) {
                        legend->AddEntry((TObject*)0, "Insufficient entries for reliable fit", "");
                    } else {
                        legend->AddEntry((TObject*)0, Form("Mean = %.5g #pm %.5g cm", mean, mean_err), "");
                        legend->AddEntry((TObject*)0, Form("Sigma = %.5g #pm %.5g cm", sigma, sigma_err), "");
                    }

                    legend->SetTextSize(0.04);
                    legend->Draw();
                    c2->Print("/eos/ams/user/s/selu/mdst/tianye/transformer/pdf/normalize_be7.pdf");

                    // clean up
                    delete proj;
                }
            
            
        }

        c2->Print("/eos/ams/user/s/selu/mdst/tianye/transformer/pdf/normalize_be7.pdf]");

    normalize_para->cd();
    for(int i=1;i<8;i++){
        for(int j=0;j<3;j++){
            h_mean_be7[j][i-1]->Write();
            h_mean_be9[j][i-1]->Write();
            h_mean_be10[j][i-1]->Write();
            h_sigma_be7[j][i-1]->Write();
            h_sigma_be9[j][i-1]->Write();
            h_sigma_be10[j][i-1]->Write();
            h_be7[j][i-1]->Write();
            h_be9[j][i-1]->Write();
            h_be10[j][i-1]->Write();
        }
    }
    normalize_para->Close();






}