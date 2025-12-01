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

void fda_dis_check(){

    TFile* f = TFile::Open("/eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_fda_norm/li_norm2D.root");
    TCanvas* c1 = new TCanvas("c1","c1",800,600);
    c1->Print("/eos/ams/user/s/selu/mdst/tianye/pdf/fda_dis_check.pdf[");
    for(int i=13;i<20;i++){
        TH1D* h = (TH1D*)f->Get(Form("hProj_%d_4_fda_1_0_1_2",i));
        TF1* f1 = (TF1*)f->Get(Form("gaussFit_%d_4_fda_1_0_1_2",i));
        //TF1* f1 = new TF1("f1","gaus",-5,5);
        //f1->SetParameters(h->GetMaximum(),h->GetMean(),h->GetRMS());
        //h->Fit(f1,"RQ");
        c1->Clear();
        cout<<h->GetEntries()<<endl;
        h->Draw();
        f1->Draw("same");
        c1->Print("/eos/ams/user/s/selu/mdst/tianye/pdf/fda_dis_check.pdf");
    }
    c1->Print("/eos/ams/user/s/selu/mdst/tianye/pdf/fda_dis_check.pdf]");

}