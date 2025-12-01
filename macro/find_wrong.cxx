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

void find_wrong(){

for(int i=1;i<3901;i+=10){
    TFile *f1= TFile::Open(Form("/eos/ams/user/s/selu/mdst/tianye/isotope_paper/root/isotope_llr/Be./be10_%d_%d_llr.root",i,i+9));
    TH1D *h1=(TH1D*)f1->Get("hist099001");
    int nbiny=h1->GetNbinsY();
    if (nbiny==1000)
    cout<<i<<endl;
}

}