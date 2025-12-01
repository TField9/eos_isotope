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


void prob_draw(TString input, TString output) {
    TChain *c_be7 = new TChain("tree");
    TChain *c_be9 = new TChain("tree");
    TChain *c_be10 = new TChain("tree");

    const double Be_bins[] = {
        0.27, 0.33, 0.41, 0.49, 0.59,
        0.70, 0.82, 0.96, 1.11, 1.28, 1.47, 1.68, 1.91, 2.16,
        2.44, 2.73, 3.06, 3.41, 3.79, 4.20, 4.65, 5.14, 5.64,
        6.18, 6.78, 7.42, 8.12, 8.86, 9.66, 10.51
    };

    TFile *f_in = TFile::Open(input);
    TFile *f_out = new TFile(output, "RECREATE");
    TTree *tree = (TTtree*)f_in->Get("tree");

    // 修复2：正确声明TString数组
    const char* detec[] = {"tof", "naf", "agl"};

    // 创建直方图
    f_out->cd();
    TH2D* h_be_prob[3];
    TH2D* h_be_mass[3];

    for (int i = 0; i < 3; i++) {
        h_be_prob[i] = new TH2D(Form("h_be_%d_prob", i), Form("h_be_%d_prob", i), 29, Be_bins, 100, 0, 1);
        h_be_mass[i] = new TH2D(Form("h_be_%d_mass", i), Form("h_be_%d_mass", i), 29, Be_bins, 100, 0, 0.3);
    }

    const Int_t nEntries_be = tree->GetEntries();

    DST sdst_be;
    sdst_be.SetAddress(tree);


    // 修复3：使用正确的变量名（sdst_be7.fEk 而不是 fEk）
    for (int i = 0; i < nEntries_be; i++) {
        tree->GetEntry(i);
        if (i % 10000 == 0) {
            std::cout << "Processed Be event: " << i + 1 << "/" << nEntries_be << "\r";
            std::cout.flush();
        }

        // require track bit
        if (!((sdst_be.fTrack.bith >> 4) & 0x1)) continue;

        // detector 0 (tof)
        if (sdst_be.fEk.Ek[0] > 0.25 && sdst_be.fEk.Ek[0] < 1.5) {
            h_be_prob[0]->Fill(sdst_be.fEk.Ek[0], sdst_be.fTransformer_output.trans_pdf[0][2]);
            h_be_mass[0]->Fill(sdst_be.fEk.Ek[0], sdst_be.fEk.rev_mass[0]);
        }

        // detector 1 (naf) when richcut_flag[1] == 0
        if (sdst_be.fRichcut_flag.richcut_flag[1] == 0) {
            h_be_prob[1]->Fill(sdst_be.fEk.Ek[1], sdst_be.fTransformer_output.trans_pdf[1][2]);
            h_be_mass[1]->Fill(sdst_be.fEk.Ek[1], sdst_be.fEk.rev_mass[1]);
        }

        // detector 2 (agl) when richcut_flag[1] == 1
        if (sdst_be.fRichcut_flag.richcut_flag[1] == 1) {
            h_be_prob[2]->Fill(sdst_be.fEk.Ek[1], sdst_be.fTransformer_output.trans_pdf[2][2]);
            h_be_mass[2]->Fill(sdst_be.fEk.Ek[1], sdst_be.fEk.rev_mass[1]);
        }
    }

    
    // 保存所有直方图
    f_out->cd();
    for (int i = 0; i < 3; i++) {
        h_be_prob[i]->Write();
        h_be_mass[i]->Write();
    }
    
    f_out->Write();
    f_out->Close();
    std::cout << "All histograms saved to " << f_out->GetName() << std::endl;
}
int main(int argc, char** argv) {
    prob_draw(argv[1], argv[2]);
    cout<<"done"<<endl;
    return 0;


}