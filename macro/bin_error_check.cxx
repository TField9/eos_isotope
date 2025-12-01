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
#include <cmath>
#include <vector>
#include <functional>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"


void bin_error_check(){

    TH1D* no_weight= new TH1D("no_weight","no_weight",2,-1,1);
    TH1D* with_weight= new TH1D("with_weight","with_weight",2,-1,1);

    no_weight->Fill(-0.5);
    no_weight->Fill(-0.5);

    with_weight->Fill(-0.5,0.2);
     cout<<"With weight hist error: "<<with_weight->GetBinError(1)<<endl;
    with_weight->Fill(-0.5,0.1);

    cout<<"No weight hist error: "<<no_weight->GetBinError(1)<<endl;
    cout<<"With weight hist error: "<<with_weight->GetBinError(1)<<endl;







}