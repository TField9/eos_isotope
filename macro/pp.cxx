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

void pp() {
    
    TFile *f1 = new TFile("/eos/ams/user/s/selu/mdst/tianye/root/normalize_parameter/para_store.root","read");
    TFile *f2= new TFile("/eos/ams/user/s/selu/mdst/tianye/root/normalize_parameter/normalize_fit_function.root","recreate");
    TH1D *h_mean_be7[3][7];
    TH1D *h_mean_be9[3][7];
    TH1D *h_mean_be10[3][7];

    TH1D *h_sigma_be7[3][7];
    TH1D *h_sigma_be9[3][7];
    TH1D *h_sigma_be10[3][7];

    for(int i=1;i<8;i++){
        for(int j=0;j<3;j++){
            h_mean_be7[j][i-1]=(TH1D*)f1->Get(Form("mean_be7_%d_%d",j,i));
            h_mean_be9[j][i-1]=(TH1D*)f1->Get(Form("mean_be9_%d_%d",j,i));
            h_mean_be10[j][i-1]=(TH1D*)f1->Get(Form("mean_be10_%d_%d",j,i));

            h_sigma_be7[j][i-1]=(TH1D*)f1->Get(Form("sigma_be7_%d_%d",j,i));
            h_sigma_be9[j][i-1]=(TH1D*)f1->Get(Form("sigma_be9_%d_%d",j,i));
            h_sigma_be10[j][i-1]=(TH1D*)f1->Get(Form("sigma_be10_%d_%d",j,i));
        }
    }

    
    // fit each histogram with an exponential and draw 21 canvases (j=0..2, i=1..7)
    std::vector<TCanvas*> canvases;
    for (int i = 1; i < 8; ++i) {
        for (int j = 0; j < 3; ++j) {
            TH1D *hm7  = h_mean_be7[j][i-1];
            TH1D *hm9  = h_mean_be9[j][i-1];
            TH1D *hm10 = h_mean_be10[j][i-1];

            TH1D *hs7  = h_sigma_be7[j][i-1];
            TH1D *hs9  = h_sigma_be9[j][i-1];
            TH1D *hs10 = h_sigma_be10[j][i-1];

        }
    }

    TF1 *f_mean_be7[3][7];
    TF1 *f_mean_be9[3][7];
    TF1 *f_mean_be10[3][7];

    TF1 *f_sigma_be7[3][7];
    TF1 *f_sigma_be9[3][7];
    TF1 *f_sigma_be10[3][7];

    {


        for (int i = 1; i < 8; ++i) {
            for (int j = 0; j < 3; ++j) {
                // helpers to access histograms from existing variables in scope
                TH1D *hist_mean7  = h_mean_be7[j][i-1];
                TH1D *hist_mean9  = h_mean_be9[j][i-1];
                TH1D *hist_mean10 = h_mean_be10[j][i-1];

                TH1D *hist_sigma7  = h_sigma_be7[j][i-1];
                TH1D *hist_sigma9  = h_sigma_be9[j][i-1];
                TH1D *hist_sigma10 = h_sigma_be10[j][i-1];

                auto fitHist = [&](TH1D* h, TF1*& dest, const char* nameBase) {
                    if (!h) { dest = nullptr; return; }
                    double xmin = h->GetXaxis()->GetXmin();
                    double xmax = h->GetXaxis()->GetXmax();
                    TString name = Form("%s_%d_%d", nameBase, j, i);
                    TF1 *f = new TF1(name, "[0]*exp([1]*x)+[2]", xmin, xmax);
                    f->SetNpx(1000);
                    h->GetYaxis()->SetRangeUser(-0.06, 0.06);
                    f->SetParameter(0, h->GetMaximum());
                    f->SetParameter(1, -0.1);
                    h->Fit(f, "RQ"); // R: use range, Q: quiet
                    dest = f;
                };

                fitHist(hist_mean7,  f_mean_be7[j][i-1],  "f_mean_be7");
                fitHist(hist_mean9,  f_mean_be9[j][i-1],  "f_mean_be9");
                fitHist(hist_mean10, f_mean_be10[j][i-1], "f_mean_be10");

                fitHist(hist_sigma7,  f_sigma_be7[j][i-1],  "f_sigma_be7");
                fitHist(hist_sigma9,  f_sigma_be9[j][i-1],  "f_sigma_be9");
                fitHist(hist_sigma10, f_sigma_be10[j][i-1], "f_sigma_be10");
            }
        }
        

        TCanvas *c1= new TCanvas("c1","c1",1600,1000);
        c1->Divide(3,2);
        c1->Print("/eos/ams/user/s/selu/mdst/tianye/transformer/pdf/fit_normalize.pdf[");
        for(int i=1;i<8;i++){
            for(int j=0;j<3;j++){
                c1->Clear();
                c1->Divide(3,2);
                gStyle->SetOptStat(0);

                // top row: mean histograms + fits
                c1->cd(1);
                if (h_mean_be7[j][i-1]) {
                    h_mean_be7[j][i-1]->SetMarkerColor(kRed);
                    h_mean_be7[j][i-1]->SetMarkerStyle(20);
                    h_mean_be7[j][i-1]->Draw("E");
                }
                if (f_mean_be7[j][i-1]) {
                    f_mean_be7[j][i-1]->SetLineColor(kBlack);
                    f_mean_be7[j][i-1]->SetLineWidth(2);
                    f_mean_be7[j][i-1]->Draw("same");
                }

                c1->cd(2);
                if (h_mean_be9[j][i-1]) {
                    h_mean_be9[j][i-1]->SetMarkerColor(kBlue);
                    h_mean_be9[j][i-1]->SetMarkerStyle(20);
                    h_mean_be9[j][i-1]->Draw("E");
                }
                if (f_mean_be9[j][i-1]) {
                    f_mean_be9[j][i-1]->SetLineColor(kBlack);
                    f_mean_be9[j][i-1]->SetLineWidth(2);
                    f_mean_be9[j][i-1]->Draw("same");
                }

                c1->cd(3);
                if (h_mean_be10[j][i-1]) {
                    h_mean_be10[j][i-1]->SetMarkerColor(kGreen+2);
                    h_mean_be10[j][i-1]->SetMarkerStyle(20);
                    h_mean_be10[j][i-1]->Draw("E");
                }
                if (f_mean_be10[j][i-1]) {
                    f_mean_be10[j][i-1]->SetLineColor(kBlack);
                    f_mean_be10[j][i-1]->SetLineWidth(2);
                    f_mean_be10[j][i-1]->Draw("same");
                }

                // bottom row: sigma histograms + fits
                c1->cd(4);
                if (h_sigma_be7[j][i-1]) {
                    h_sigma_be7[j][i-1]->SetMarkerColor(kMagenta);
                    h_sigma_be7[j][i-1]->SetMarkerStyle(21);
                    h_sigma_be7[j][i-1]->Draw("E");
                }
                if (f_sigma_be7[j][i-1]) {
                    f_sigma_be7[j][i-1]->SetLineColor(kBlack);
                    f_sigma_be7[j][i-1]->SetLineWidth(2);
                    f_sigma_be7[j][i-1]->Draw("same");
                }

                c1->cd(5);
                if (h_sigma_be9[j][i-1]) {
                    h_sigma_be9[j][i-1]->SetMarkerColor(kCyan);
                    h_sigma_be9[j][i-1]->SetMarkerStyle(21);
                    h_sigma_be9[j][i-1]->Draw("E");
                }
                if (f_sigma_be9[j][i-1]) {
                    f_sigma_be9[j][i-1]->SetLineColor(kBlack);
                    f_sigma_be9[j][i-1]->SetLineWidth(2);
                    f_sigma_be9[j][i-1]->Draw("same");
                }

                c1->cd(6);
                if (h_sigma_be10[j][i-1]) {
                    h_sigma_be10[j][i-1]->SetMarkerColor(kOrange+7);
                    h_sigma_be10[j][i-1]->SetMarkerStyle(21);
                    h_sigma_be10[j][i-1]->Draw("E");
                }
                if (f_sigma_be10[j][i-1]) {
                    f_sigma_be10[j][i-1]->SetLineColor(kBlack);
                    f_sigma_be10[j][i-1]->SetLineWidth(2);
                    f_sigma_be10[j][i-1]->Draw("same");
                }

                c1->Print("/eos/ams/user/s/selu/mdst/tianye/transformer/pdf/fit_normalize.pdf");
            }
        }
        c1->Print("/eos/ams/user/s/selu/mdst/tianye/transformer/pdf/fit_normalize.pdf]");
    }
    // Save fit functions to file
    f2->cd();
    for(int i=1;i<8;i++){
        for(int j=0;j<3;j++){
            if (f_mean_be7[j][i-1])  f_mean_be7[j][i-1]->Write();
            if (f_mean_be9[j][i-1])   f_mean_be9[j][i-1]->Write();
            if (f_mean_be10[j][i-1]) f_mean_be10[j][i-1]->Write(); 
            if (f_sigma_be7[j][i-1])  f_sigma_be7[j][i-1]->Write();
            if (f_sigma_be9[j][i-1])   f_sigma_be9[j][i-1]->Write();
            if (f_sigma_be10[j][i-1]) f_sigma_be10[j][i-1]->Write(); 
        }
    }
    f2->Close();
}