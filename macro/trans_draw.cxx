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

const double Be_bins[] = {
    0.27, 0.33, 0.41, 0.49, 0.59,
    0.70, 0.82, 0.96, 1.11, 1.28, 1.47, 1.68, 1.91, 2.16,
    2.44, 2.73, 3.06, 3.41, 3.79, 4.20, 4.65, 5.14, 5.64,
    6.18, 6.78, 7.42, 8.12, 8.86, 9.66, 10.51
};

void trans_draw(){
    TFile *f_be7 = TFile::Open("/eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_plot/sdst_be7_plot.root","read");
    TFile *f_be9 = TFile::Open("/eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_plot/sdst_be9_plot.root","read");
    TFile *f_be10 = TFile::Open("/eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_plot/sdst_be10_plot.root","read");
    TFile *f_list[] = {f_be7, f_be9, f_be10};
    TH2D* h_be_prob[3][3]; // [isotope][detector]
    TH2D* h_be_mass[3][3];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            h_be_prob[i][j] = (TH2D*)f_list[i]->Get(Form("h_be_%d_prob", j));
            h_be_mass[i][j] = (TH2D*)f_list[i]->Get(Form("h_be_%d_mass", j));
        }
    }

    TH2D *h_ref = (TH2D*)h_be_prob[0][0]->Clone("h_ref");
    h_ref->Reset();
    TH1D *h_proj_prob[3][sizeof(Be_bins)/sizeof(Be_bins[0])-1];
    TH1D *h_proj_mass[3][sizeof(Be_bins)/sizeof(Be_bins[0])-1];
    //loop
    for(int i=1;i<sizeof(Be_bins)/sizeof(Be_bins[0]);i++){
        if(Be_bins[i]<1.5){
            for(int j=0;j<3;j++){
                h_proj_prob[j][i-1] = h_be_prob[j][0]->ProjectionY(Form("h_be_%d_prob_bin%d", j, i), i, i);
                h_proj_mass[j][i-1] = h_be_mass[j][0]->ProjectionY(Form("h_be_%d_mass_bin%d", j, i), i, i);
            }
        }
        else if(Be_bins[i]<3.5){
            for(int j=0;j<3;j++){
                h_proj_prob[j][i-1] = h_be_prob[j][1]->ProjectionY(Form("h_be_%d_prob_bin%d", j, i), i, i);
                h_proj_mass[j][i-1] = h_be_mass[j][1]->ProjectionY(Form("h_be_%d_mass_bin%d", j, i), i, i);
            }
        }
        else{
            for(int j=0;j<3;j++){
                h_proj_prob[j][i-1] = h_be_prob[j][2]->ProjectionY(Form("h_be_%d_prob_bin%d", j, i), i, i);
                h_proj_mass[j][i-1] = h_be_mass[j][2]->ProjectionY(Form("h_be_%d_mass_bin%d", j, i), i, i);
            }
        }
    }

    TCanvas *c = new TCanvas("c", "Be Projections", 1200, 800);
    c->Print("/eos/ams/user/s/selu/mdst/tianye/transformer/pdf/Be_Projections.pdf");
    c->Divide(2,1);

    for(int i=0; i<28; i++){
        c->Clear();
        c->Divide(2,1);
        gStyle->SetOptStat(0);
        c->cd(1);
        h_proj_mass[0][i]->SetLineColor(kRed);
        h_proj_mass[1][i]->SetLineColor(kBlue);
        h_proj_mass[2][i]->SetLineColor(kGreen+2);
        h_proj_mass[0][i]->SetTitle(Form("Mass, %.2f < Ek/n < %.2f", Be_bins[i], Be_bins[i+1]));
        for(int j=0; j<3; j++) {
            h_proj_mass[j][i]->Rebin(1);
        }
        double max_mass = std::max({h_proj_mass[0][i]->GetMaximum(), h_proj_mass[1][i]->GetMaximum(), h_proj_mass[2][i]->GetMaximum()});
        for(int j=0; j<3; j++) {
            h_proj_mass[j][i]->SetMaximum(1.2 * max_mass);
        }
        h_proj_mass[0][i]->GetXaxis()->SetRangeUser(0.05, 0.25);

        h_proj_mass[0][i]->Draw("hist");
        h_proj_mass[1][i]->Draw("hist same");
        h_proj_mass[2][i]->Draw("hist same");

        auto legend_mass = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend_mass->AddEntry(h_proj_mass[0][i], "Be7", "l");
        legend_mass->AddEntry(h_proj_mass[1][i], "Be9", "l");
        legend_mass->AddEntry(h_proj_mass[2][i], "Be10", "l");
        legend_mass->Draw();

        c->cd(2);
        h_proj_prob[0][i]->SetLineColor(kRed);
        h_proj_prob[1][i]->SetLineColor(kBlue);
        h_proj_prob[2][i]->SetLineColor(kGreen+2);
        h_proj_prob[0][i]->SetTitle(Form("transformer be10 Prob, %.2f < Ek/n < %.2f", Be_bins[i], Be_bins[i+1]));
         for(int j=0; j<3; j++) {
            h_proj_prob[j][i]->Rebin(2);
        }
        double max_prob = std::max({h_proj_prob[0][i]->GetMaximum(), h_proj_prob[1][i]->GetMaximum(), h_proj_prob[2][i]->GetMaximum()});
        for(int j=0; j<3; j++) {
            h_proj_prob[j][i]->SetMaximum(1.2 * max_prob);
        }
        
        h_proj_prob[0][i]->Draw("hist");
        h_proj_prob[1][i]->Draw("hist same");
        h_proj_prob[2][i]->Draw("hist same");
        auto legend_prob = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend_prob->AddEntry(h_proj_prob[0][i], "Be7", "l");
        legend_prob->AddEntry(h_proj_prob[1][i], "Be9", "l");
        legend_prob->AddEntry(h_proj_prob[2][i], "Be10", "l");
        legend_prob->Draw();

        c->Print("/eos/ams/user/s/selu/mdst/tianye/transformer/pdf/Be_Projections.pdf(");
    }
    c->Print("/eos/ams/user/s/selu/mdst/tianye/transformer/pdf/Be_Projections.pdf)");

}