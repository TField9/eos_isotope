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
#include "sdst.h"


int check_bit(int bit){
    int sum=0;
    while(bit){
        sum += bit&0x1;
        bit=bit>>1;
    }
    return sum;


    
}

void he_check(){

    TChain *ch = new TChain("ftree");
    DST *dst = new DST();
    //dst->SetAddress(ch);
    for(int i=1;i<=100;i=i+10){
        TFile *f = TFile::Open(Form("/eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_llr_fillTree/He./he4_%d_%d_llr.root",i,i+9));
        ch->Add(f->GetName());
        cout<<"Adding file: "<<f->GetName()<<endl;
    }
    dst->SetAddress(ch);
    TH2D * m_llr2D = new TH2D("m_llr2D","m_llr2D",100,0,1,100,0,1);
    TH1D *h_chi2y = new TH1D("h_chi2y","h_chi2y",100,0,10);
    TH1D *fda_up = new TH1D("fda_up","fda_up",100,-0.2,0.2);
    TH1D *fda_down = new TH1D("fda_down","fda_down",100,-0.2,0.2);
    TH1D *y_hits = new TH1D("y_hits","y_hits",7,1,8);
    for(int j=0;j<ch->GetEntries();j++){

        //TCut cut= "dst->fEk.Ek[1]>2 && dst->fEk.Ek[1]<2.5 &&dst->fEstim.llr2[1]>0.2";
        if(j%10000==0) cout<<"Processing entry: "<<j<<"/"<<ch->GetEntries()<<endl;
        ch->GetEntry(j);
        if(dst->fEk.Ek[1]>2 && dst->fEk.Ek[1]<2.5&&dst->fEk.rev_mass[1]<0.1){
            //cout<<dst->fEstim.llr2[1]<<endl;
            m_llr2D->Fill(dst->fEk.rev_mass[1],dst->fEstim.llr2[1]);
            h_chi2y->Fill(dst->fTrack.csy2);
            fda_up->Fill(dst->fEstim.fDa2[1][1][0]);
            fda_down->Fill(dst->fEstim.fDa2[1][1][1]);
            y_hits->Fill(check_bit(dst->fTrack.bith));
        }
    }

    TH1D *hProj_mass=m_llr2D->ProjectionX("hProj_mass",1,100);
    TH1D *hProj_llr=m_llr2D->ProjectionY("hProj_llr",1,100);
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    c1->Print("/eos/ams/user/s/selu/mdst/tianye/pdf/he3_llr2D.pdf[");
    c1->Divide(2,2);
    c1->cd(1);
    m_llr2D->Draw("COLZ");
    c1->cd(2);
    hProj_mass->Draw();
    c1->cd(3);
    hProj_llr->Draw();
    c1->cd(4);
    h_chi2y->Draw();
    c1->Print("/eos/ams/user/s/selu/mdst/tianye/pdf/he3_llr2D.pdf");
    c1->Clear();
    c1->Divide(2,2);
    c1->cd(1);
    fda_up->Draw();
    c1->cd(2);
    fda_down->Draw();
    c1->cd(3);
    y_hits->Draw();
    c1->Print("/eos/ams/user/s/selu/mdst/tianye/pdf/he3_llr2D.pdf");
    c1->Print("/eos/ams/user/s/selu/mdst/tianye/pdf/he3_llr2D.pdf]");






}