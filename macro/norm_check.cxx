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

void norm_check(){
    TFile *f = new TFile("/eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_residual_norm/sdst_be7_61_70_residual_norm.root","read");
    DST *dst = new DST();
    
    TTree *t = (TTree*)f->Get("tree");
    dst->SetAddress(t);
    Float_t y_norm[3][7];
    t->SetBranchAddress("y_residual", y_norm);
    Long64_t nentries = t->GetEntries();
    TH1D* h1= new TH1D("h1","y_residual_norm",100,-10,10);
    float minv = 1e30f;
    float maxv = -1e30f;
    for(Long64_t i=0;i<nentries;i++){
        t->GetEntry(i);
        if(dst->fEk.Ek[0]<0.25 || dst->fEk.Ek[0]>1.5) continue;
        Float_t val = y_norm[0][6];
        h1->Fill(val);
        if(val < minv) minv = val;
        if(val > maxv) maxv = val;
    }
    if(minv==1e30f){
        std::cout << "No entries passed the selection" << std::endl;
    } else {
        std::cout << "y_norm[0][6] min = " << minv << "  max = " << maxv << std::endl;
    }
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    h1->Draw();
    c1->Print("/eos/ams/user/s/selu/mdst/tianye/transformer/pdf/norm_check.pdf");



}