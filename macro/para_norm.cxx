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

void para_norm() {
    
    TFile *f1 = new TFile("/eos/ams/user/s/selu/mdst/tianye/root/normalize_parameter/normalize_para.root","read");

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

            // if nothing to draw for this (j,i), skip
            if (!hm7 && !hm9 && !hm10 && !hs7 && !hs9 && !hs10) continue;

            TCanvas *c = new TCanvas(Form("c_mean_sigma_j%d_i%d", j, i),
                                     Form("j=%d i=%d", j, i), 800, 900);

            // top pad for means
            TPad *padTop = new TPad(Form("padTop_j%d_i%d", j, i), "", 0.0, 0.45, 1.0, 1.0);
            padTop->SetBottomMargin(0.02);
            padTop->Draw();
            padTop->cd();

            TLegend *legTop = new TLegend(0.65, 0.65, 0.88, 0.88);
            legTop->SetBorderSize(0);

            bool firstMeanDrawn = false;
            auto draw_mean = [&](TH1D *h, Color_t col, const char *label) {
                if (!h) return;
                h->SetLineColor(col);
                h->SetMarkerColor(col);
                h->SetMarkerStyle(20);
                if (!firstMeanDrawn) {
                    h->SetTitle(Form("Means (j=%d i=%d);x;Entries", j, i));
                    h->Draw("E");
                    firstMeanDrawn = true;
                } else {
                    h->Draw("E SAME");
                }
                // fit exponential [0]*exp([1]*x) in histogram x range
                double xmin = h->GetXaxis()->GetXmin();
                double xmax = h->GetXaxis()->GetXmax();
                TF1 *f = new TF1(Form("f_mean_%d_%d_%s", j, i, label), "[0]*exp([1]*x)", xmin, xmax);
                f->SetLineColor(col);
                // initial guesses
                f->SetParameter(0, h->GetMaximum());
                f->SetParameter(1, -0.1);
                h->Fit(f, "RQ"); // R: use range, Q: quiet
                f->Draw("SAME");
                legTop->AddEntry(h, label, "lep");
            };

            draw_mean(hm7, kRed,   "mean_be7");
            draw_mean(hm9,  kBlue,  "mean_be9");
            draw_mean(hm10, kGreen, "mean_be10");

            legTop->Draw();

            // bottom pad for sigmas
            c->cd();
            TPad *padBot = new TPad(Form("padBot_j%d_i%d", j, i), "", 0.0, 0.0, 1.0, 0.45);
            padBot->SetTopMargin(0.02);
            padBot->SetBottomMargin(0.15);
            padBot->Draw();
            padBot->cd();

            TLegend *legBot = new TLegend(0.65, 0.55, 0.88, 0.88);
            legBot->SetBorderSize(0);

            bool firstSigmaDrawn = false;
            auto draw_sigma = [&](TH1D *h, Color_t col, const char *label) {
                if (!h) return;
                h->SetLineColor(col);
                h->SetMarkerColor(col);
                h->SetMarkerStyle(21);
                if (!firstSigmaDrawn) {
                    h->SetTitle(Form("Sigmas (j=%d i=%d);x;Entries", j, i));
                    h->Draw("E");
                    firstSigmaDrawn = true;
                } else {
                    h->Draw("E SAME");
                }
                double xmin = h->GetXaxis()->GetXmin();
                double xmax = h->GetXaxis()->GetXmax();
                TF1 *f = new TF1(Form("f_sigma_%d_%d_%s", j, i, label), "[0]*exp([1]*x)", xmin, xmax);
                f->SetLineColor(col);
                f->SetParameter(0, h->GetMaximum());
                f->SetParameter(1, -0.1);
                h->Fit(f, "RQ");
                f->Draw("SAME");
                legBot->AddEntry(h, label, "lep");
            };

            draw_sigma(hs7,  kRed,   "sigma_be7");
            draw_sigma(hs9,  kBlue,  "sigma_be9");
            draw_sigma(hs10, kGreen, "sigma_be10");

            legBot->Draw();

            // finalize
            c->Update();
            canvases.push_back(c);
        }
    }

    // save all canvases into a single multi-page PDF
    const char *outPdf = "mean_sigma_all.pdf";
    if (canvases.empty()) {
        std::cout << "No canvases to save to PDF\n";
    } else {
        for (size_t idx = 0; idx < canvases.size(); ++idx) {
            if (idx == 0) {
                canvases[idx]->Print(Form("%s(", outPdf)); // open PDF
            } else if (idx == canvases.size() - 1) {
                canvases[idx]->Print(Form("%s)", outPdf)); // close PDF
            } else {
                canvases[idx]->Print(outPdf); // append page
            }
        }
    }

    // optional: clean up canvases if desired
    // for (auto c : canvases) { delete c; }



}