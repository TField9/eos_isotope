
enum   EStat1 { BadTrig = 0x0001, BadTRD   = 0x0002, BadRTI   = 0x0004,
		InSAA   = 0x0008, FirstSec = 0x0010, LastSec  = 0x0020,
		IsInL1  = 0x0040, IsInL9   = 0x0080, IsInEcal = 0x0100,
		HasL1   = 0x0200, HasL2    = 0x0400, HasL9    = 0x0800,
		TrdInFS = 0x1000, TrdInTr  = 0x2000, FirstRTI = 0x4000,
		NoPTrk  = 0x8000,
		TrdInL1 = 0x10000, CutTk2nd = 0x20000 };
enum   EStat2 {};

struct Status { UInt_t status;     // AMSEventR::fStatus
                UInt_t ustat; };   // Estat1 | (Estat2 << 16)

struct Header { UInt_t run, event, ient, utime, phpat, jmpat;
                Int_t ntrdhit, ntrdseg, ntrhit, ntofcls, 
		      nanticls, nrichhit, necalhit;
                Int_t error; };
struct RTI    { Float_t theta, phi, zenith, cfi, lf, dl1, dl9, cf; };

struct Part   { Float_t mon, emom, mass, emass, beta, ebeta; };
struct Track  { Int_t   bith, bitx;  // Y,X bit pattern
                Int_t   bitr;        // Recovered bit pattern
                Float_t rgt[4], csqx[4], csqy[4]; // 0:In 1:L1 2:L9 3:FS
                Float_t qin, ql1, ql9;
                Float_t p0x, p0y, theta, phi; // At Z=0
                Float_t coox[9], cooy[9], cooz[9]; // GBL, no MS, has Energy loss
                Float_t cootheta[9], coophi[9]; //// GBL, no MS, has Energy loss
                Float_t xglob_JN[9], yglob_JN[9], zglob_JN[9]; // inner sGlobal coordinates (JN version)
                Float_t resx[2], resy[2];     // 0:L1 1:L9 excl.residual
                Float_t hity[4];              // L1 PG,MD, L9 PG,MD
                Int_t   qsta[2];              // GetQStatus() for 0:L1 and 1:L9
                Float_t qrms;                 // GetInnerQ_all.RMS
                Float_t rgt2[3], csx2, csy2;  // Chikanian,Up,Low
                Int_t   itp, addl;            // Max span, AddLost
                Int_t   rec;                  // GetRecType()
              };
struct TrHit  { Int_t   nhit[6];       // 0:L1xy 1:L1y 2:L9xy 3:L9y 4:Ixy 5:Iy
                Float_t qmax[6];       // 0:L1xy 1:L1y 2:L9xy 3:L9y 4:Ixy 5:Iy
                Float_t qli [7];       // TrTrack::Qlayer (Inner)
                Float_t dx[2], dy[4];  // QmaxHit-Track 0:L1 1:L9 +2:OnlyY
                Int_t   qsi [7];       // GetQStatus()
                Int_t   nh10[7];       // IsTrackPickingUpNoise(10, 0, 0, nh10)
                Float_t rata[7];       // GetTrackerRawSignalRatio(j+2);
                Float_t tkfd[7];       // GetTkFeetDist(j+2)
                Float_t qs10[7][20];   // 20-Strip Q around hit
                Float_t xglob[7], yglob[7], zglob[7]; // inner sGlobal coordinates
              };
struct TrCls  { Int_t   tkml[9];
                Float_t xcog[9], ycog[9], mdist[9];
                Int_t   flag[9];
                Float_t eta [9], sig [9];
              };

struct Beta   { Int_t pattern; Float_t beta; };
struct BetaH  { Int_t pattern, pbit;                     // Default BetaH
                Float_t beta, chi2t, chi2c, q, ql[4];
                Int_t clsn[4], z[3]; Float_t p[3];
                Int_t type; Float_t qup, qlow;
              };
struct BetaHs { Int_t pattern, pbit;                     // TrTrack independent
                Float_t beta, chi2t, chi2c, q, ql[4]; };
struct TofHit { Int_t   nhit[4];
                Float_t qsum[4], qmax[4];
                Float_t tmin[4], tmax[4]; };

struct Trd    { Float_t coo[3], theta, phi, q;
                Float_t dtrk[3];       // 0:dx 1:dy 2:angle(deg) with TrTrack
                Float_t ql1m, dl1m[2]; // Charge, dx, dy of the closest L1 hit
                Int_t nsegt;
              };
struct TrdK   { Int_t nhits; Float_t llr[3], q, llre[3];
                Int_t nhitt; Float_t llrt[3];
                Int_t nhitd; Float_t llrd[3];
              };
struct TrdG   { Int_t nfit;  Float_t gamma, lkh; };
struct TrdV   {
                Int_t nvtx, ntrk, nhit;
                Float_t chi2, x, y, z;
};
struct TrdHit { Int_t  layer[25];
                Float_t plen[25], amp[25], plmc[25], trfr[25];
	      };
struct Rich   { Int_t status, nhits[3];
                Float_t beta, q, dist, prob, coo[3], npe[3];
                Int_t tile, pmts;
              };
struct Ecal   { Float_t cog[3], dir[3], s13r, bdt, energy[3]; // 0:D 1:C 2:E
                Float_t dtrk[3];       // 0:dx 1:dy 2:angle(deg) with TrTrack
                Float_t ql9m, dl9m[2]; // Charge, dx, dy of the closest L9 hit
                Int_t catl;
                Float_t enew[2];       // GetCorrectedEnergy
                Float_t dz[2];
              };

struct EcalH  { Int_t apex;
                Float_t csq, rrec, smax, rnt, emip, lmip;
                Float_t edep[10];
};
struct EcalX {
  Int_t ns; // number of shower, <0 hadronic
  Float_t etot[2], ll[2], bdt;
  Int_t stat[3];
  Float_t eng[3];
  Float_t start[3];
  Float_t locx[3], locy[3], locz[3];
  Float_t dirx[3], diry[3], dirz[3];
  Int_t nsat;
  Float_t satloss;
  Float_t cell[18][72];
};



struct Estim  {   
  Float_t ll1[3][2], ll2[3][2];//[iso][fit algo:A/gbl]
  Float_t chisq1[3][2], chisq2[3][2]; // [fit algo]
  Float_t pa1[3][2][4], pa2[3][2][4]; // [iso][fit algo][up/low]
  Float_t fDa1[3][2][2], fDa2[3][2][2]; // [iso][fit algo][up/low]
  Float_t test_pa1[3][2][2], test_pa2[3][2][2]; // [iso][up/low]
  Float_t res1_ll[3][2], res2_ll[3][2]; // [iso][fit algo]
  Float_t llr1[2],llr2[2];//[fit algo]
  Float_t LLR1[2], LLR2[2];//[fit algo]
  Float_t llr1_comb[2], llr2_comb[2]; // angle & residual combined likelihood values
  Float_t LLR1_comb[2], LLR2_comb[2]; // angle & residual combined likelihood values
  Float_t llr1_withres, llr2_withres; // LLR for 1,2
  Float_t LLR1_withres,LLR2_withres;
  Float_t llr1_nores, llr2_nores; // LLR for 1,2
  Float_t LLR1_nores, LLR2_nores;

  // Float_t x_A1[3][7], y_A1[3][7], z_A1[3][7]; // A track
  // Float_t x_G1[3][7], y_G1[3][7], z_G1[3][7]; //GBL track
  // Float_t theta_yz_A1[3][7], theta_yz_G1[3][7]; // A track, GBL track
  // Float_t theta_xz_A1[3][7], theta_xz_G1[3][7]; // A track, GBL track

  // Float_t theta_yz_A2[3][7], theta_yz_G2[3][7]; // A track, GBL track
  // Float_t theta_xz_A2[3][7], theta_xz_G2[3][7]; // A track, GBL track
  // Float_t x_A2[3][7], y_A2[3][7], z_A2[3][7]; // A track
  // Float_t x_G2[3][7], y_G2[3][7], z_G2[3][7]; //GBL track

};

struct Mass   { Float_t b1, b2, b3; // TrMass::GetBeta with flag= 1, 10, 100
                Float_t npk, m;     // TrMass::GetNpick, mass(Chikanian,BetaH)
                Float_t v[16], mql; // TrMass::GetMQL
                Int_t nh;           // Number of inner tracker hits
                Float_t ll[2];      // TrMass::GetLL with (1,1),(1,2)
                Float_t bl2, ru;    // TrMass::GetBL2, TrMass::Rndm
                Int_t zl[4];
                Float_t pz[4], betas;
              };
struct Refit  { Float_t tx[9], ty[9], tz[9]; 
                Float_t fx[9], fy[9], fz[9];
                Int_t nsx[9], nsy[9];};


struct Extra{
  Int_t ncls,Used,UsedM,NHits;
  Float_t Beta,EBeta,BetaW;
  Float_t NpCol, NpLkh, NpExp, NDist, Prob, Theta, ETheta;
  Float_t x[100],y[100],z[100],hpe[100],hb[100];
  Int_t NBeta[5];
  Float_t MBeta[5],SBeta[5];
  Float_t rll[5];
  Float_t rll_gbl[5];
  Float_t rgt0, rgtt;
  Float_t csx,csy,csq,ndx,ndy,ndf;
  Int_t scl;
  Float_t qyjin,qyjl1,qyjl9;
  //Float_t xtm,pl,rc[4];
  //Int_t nc[2];
};

struct Vertex { Float_t theta, phi, zenith;
                Float_t rgt[2], csqx[2], csqy[2], coo[7], ang[4]; };

struct MCinfo { Float_t coo[3], dir[3], rgt, mcin[2];
                Int_t   part[5]; // 0:L1 1:L2 2:L56 3:L78 4:L9
                Float_t trgt[5]; // 0:L1 1:L2 2:L56 3:L78 4:L9
                Int_t   patmc, npmc;
                Float_t trcx[9], trcy[9], trcz[9];
                Float_t dir1[3], dir9[3];
                Float_t dir3, dir5;
                Float_t rcoo[3], rdir[3], rmom, rchg, rmas;
                Int_t rpid; 
                Float_t mci2[7], cfw[2];
                Int_t rseed[2];
              };

struct EK     {
              Float_t Ek[2], rev_mass[2]; 
};

struct pred_Be_prob{
              Float_t pred_Be_prob[3]; // 0: Be7, 1: Be9, 2: Be10
};

struct Richcut_flag{
  Int_t richcut_flag[2];
};

struct Transformer_output{
  Float_t trans_pdf[3][3]; //1d for detector, 2d for Be7,9,10
};

class DST {
public:
  Status fStatus; Header fHeader; RTI    fRTI;
  Part   fPart;   Track  fTrack;  TrHit  fTrHit;  TrCls  fTrCls;
  Beta   fBeta;   BetaH  fBetaH;  BetaHs fBetaHs; TofHit fTofHit;
  Trd    fTrd;    TrdK   fTrdK;   TrdG   fTrdG;   TrdHit fTrdHit; TrdV fTrdV;
  Rich   fRich;   Ecal   fEcal;   EcalH  fEcalH;  Vertex fVertex; 
  MCinfo fMCinfo; Estim  fEstim;  Refit  fRefit; Mass   fMass;
  Extra fExtra; EcalX fEcalX; EK fEk; Richcut_flag fRichcut_flag; Transformer_output fTransformer_output;

  void Clear() { Int_t *ptr = (Int_t *)this;
                 Int_t size = sizeof(DST)/sizeof(Int_t);
		 for (Int_t i = 0; i < size; i++) ptr[i] = 0;
               }

  void Branch(TTree *tree, Int_t mode = 1) {
    tree->Branch("status", &fStatus, "status/i:ustat/i");
    tree->Branch("header", &fHeader, "run/i:event/i:ient/i:"
		                     "utime/i:phpat/i:jmpat/i:"
		                     "ntrdh/I:ntrds/I:ntrh/I:ntofc/I:"
		                     "nantic/I:nrichh/I:necalh/I:error/I");
    tree->Branch("rti",    &fRTI,    "theta/F:phi/F:zenith/F:"
		                     "cfi/F:lf/F:dl1/F:dl9/F:cf/F");
    if(0)
    tree->Branch("part",   &fPart,   "mom/F:emom/F:mass/F:emass/F:"
		                     "beta/F:ebeta/F");
    tree->Branch("track",  &fTrack,  "bith/I:bitx/I:bitr/I:"
		                     "rgt[4]/F:csqx[4]/F:csqy[4]/F:"
		                     "qin/F:ql1/F:ql9/F:"
		                     "p0x/F:p0y/F:theta/F:phi/F:"
                          "coox[9]/F:cooy[9]/F:cooz[9]/F:"
                          "cootheta[9]/F:coophi[9]/F:"
                          "xglob_JN[9]/F:yglob_JN[9]/F:zglob_JN[9]/F:"
                          "resx[2]/F:resy[2]/F:"
                          "hity[4]/F:qsta[2]/F:"
                          "qrms/F:rgt2[3]/F:csx2/F:csy2/F:"
                          "itp/I:addl/I:rec/I");


  if(1)
    tree->Branch("trhit",  &fTrHit,  "nhit[6]/I:qmax[6]/F:qli[7]/F:"
		                     "hdx[2]/F:hdy[4]/F:qsi[7]/I:"
		                     "nh10[7]/I:rata[7]/F:tkfd[7]/F:qs10[7][20]/F:"
                         "xglob[7]/F:yglob[7]/F:zglob[7]/F");
  if(1)
    tree->Branch("trcls",  &fTrCls,  "tkml[9]/I:"
		                     "xcog[9]/F:ycog[9]/F:mdist[9]/F:"
		                     "flag[9]/I:eta[9]/F:sig[9]/F");

  if(0)
tree->Branch("refit",  &fRefit,  "tx[9]/F:ty[9]/F:tz[9]/F:"
                        "fx[9]/F:fy[9]/F:fz[9]/F:"
                        "nsx[9]/I:nsy[9]/I");

    tree->Branch("beta",   &fBeta,   "pattern/I:beta/F");
    tree->Branch("betah",  &fBetaH,  "pattern/I:pbit/I:beta/F:"
		                     "chi2t/F:chi2c/F:q/F:ql[4]/F:"
		                     "clsn[4]/I:z[3]/I:p[3]/F:"
                         "type/I:qup/F:qlow/F");
    if(1)
    tree->Branch("betas",  &fBetaHs, "pattern/I:pbit/I:beta/F:"
		                     "chi2t/F:chi2c/F:q/F:ql[4]/F");
    if(0)
    tree->Branch("tofhit", &fTofHit, "nhit[4]/I:qsum[4]/F:qmax[4]/F:"
		                     "tmin[4]/F:tmax[4]/F");

    tree->Branch("trd",    &fTrd,    "coo[3]/F:theta/F:phi/F:q/F:"
		                     "dtrk[3]/F:ql1m/F:dl1m[2]/F:ntseg/I");
    tree->Branch("trdk",   &fTrdK,   "nhits/I:llr[3]/F:q/F:llre[3]/F:"
		                     "nhitt/I:llrt[3]/F:"
		                     "nhitd/I:llrd[3]/F");
    if(0)
    tree->Branch("trdg",   &fTrdG,   "nfit/I:gamma/F:lkh/F");

    if (0)
    tree->Branch("trdhit", &fTrdHit, "layer[25]/I:plen[25]/F:amp[25]/F:"
		                                 "plmc[25]/F:trfr[25]/F");
    if(0)
    tree->Branch("trdv",   &fTrdV, "nvtx/I:nvtrk/I:nvhit/I:chi2/F:"
                                    "trdvx/F:trdvy/F:trdvz/F");

    tree->Branch("rich",   &fRich,   "rstat/I:rnhit[3]/I:"
		                     "beta/F:q/F:dist/F:prob/F:rcoo[3]/F:"
		                     "npe[3]/F:tile/I:pmts/I");

    if(1)
    tree->Branch("ecal",   &fEcal,   "cog[3]/F:dir[3]/F:s13r/F:bdt/F:"
		                     "energy[3]/F:dtrk[3]/F:ql1m/F:dl1m[2]/F:"
		                     "catl/I:enew[2]/F:dz[2]/F");
    if(0)
        tree->Branch("ecalx", &fEcalX, "ns/I:etot[2]/F:ll[2]/F:bdt/F:"
    "stat[3]/I:eng[3]/F:start[3]/F:"
    "locx[3]/F:locy[3]/F:locz[3]/F:"
    "dirx[3]/F:diry[3]/F:dirz[3]/F:"
    "nsat/I:satloss/F:cell[18][72]/F"
  );
    if(0)
    tree->Branch("ecalh",  &fEcalH,  "apex/I:csq/F:rrec/F:smax/F:"
		                     "rnt/F:emip/F:lmip/F:"
		                     "edep[10]/F");

    if(1)
    tree->Branch("vertex", &fVertex, "vth/F:vph/F:vzen/F:vrgt[2]/F:"
		                     "vcsqx[2]/F:vcsqy[2]/F:"
		                     "vcoo[7]/F:vang[4]/F");

    tree->Branch("mcinfo", &fMCinfo, "mcoo[3]/F:mdir[3]/F:mrgt/F:mcin[2]/F:"
		                     "part[5]/I:trgt[5]/F:patmc/I:npmc/I:"
		                     "trcx[9]/F:trcy[9]/F:trcz[9]/F:"
		                     "dir1[3]/F:dir9[3]/F:"
                         "dir3/F:dir5/F:"
                         "rcoo[3]/F:rdir[3]/F:rmom/F:rchg/F:rmas/F:rpid/I:"
                         "mci2[7]/F:cfw[2]/F:rseed[2]/I");

    tree->Branch("mass",   &fMass,   "b1/F:b2/F:b3/F:npk/F:m/F:"
		                     "v[16]/F:mql/F:"
		                     "nh/I:ll[2]/F:bl2/F:ru/F:zl[4]/I:pz[4]/F:betas/F");
    tree->Branch("estim",  &fEstim,  "ll1[3][2]/F:ll2[3][2]/F:"
                         "chisq1[3][2]/F:chisq2[3][2]/F:"
                         "pa1[3][2][4]/F:pa2[3][2][4]/F:"
                         "fDa1[3][2][2]/F:fDa2[3][2][2]/F:"
                         "test_pa1[3][2][2]/F:test_pa2[3][2][2]/F:"
                         "res1_ll[3][2]/F:res2_ll[3][2]/F:"
                         "llr1[2]/F:llr2[2]/F:LLR1[2]/F:LLR2[2]/F:"
                         "llr1_comb[2]/F:llr2_comb[2]/F:"
                         "LLR1_comb[2]/F:LLR2_comb[2]/F:"
                         "llr1_withres/F:llr2_withres/F:"
                         "LLR1_withres/F:LLR2_withres/F:"
                         "llr1_nores/F:llr2_nores/F:"
                         "LLR1_nores/F:LLR2_nores/F");
                          
                        
    tree->Branch("extra", &fExtra,"ncls/I:Used/I:UsedM/I:NHits/I:"
        "Beta/F:EBeta/F:BetaW/F:NpCol/F:NpLkh/F:NpExp/F:NDist/F:"
        "Prob/F:Theta/F:ETheta/F:x[100]/F:y[100]/F:z[100]/F:hpe[100]/F:"
        "hb[100]/F:NBeta[5]/I:MBeta[5]/F:SBeta[5]/F:"
        "rll[5]/F:rll_gbl[5]/F:"
        "rgt0/F:rgtt/F:"
        "csx/F:csy/F:csq/F:ndx/F:ndy/F:ndf/F:scl/I:qyjin/F:qyjl1/F:qyjl9/F"
        );

    tree->Branch("ek", &fEk, "Ek[2]/F:rev_mass[2]/F");

  }

  void SetAddress(TTree *tree) {
    tree->SetBranchAddress("status", &fStatus);
    tree->SetBranchAddress("header", &fHeader);
    tree->SetBranchAddress("rti",    &fRTI);
if (tree->GetBranch       ("part")){
    tree->SetBranchAddress("part",   &fPart);}
    
    tree->SetBranchAddress("track",  &fTrack);
    tree->SetBranchAddress("trhit",  &fTrHit);


if (tree->GetBranch       ("trcls"))
    tree->SetBranchAddress("trcls",  &fTrCls);
    tree->SetBranchAddress("beta",   &fBeta);
    tree->SetBranchAddress("betah",  &fBetaH);
if (tree->GetBranch       ("betas"))
    tree->SetBranchAddress("betas",  &fBetaHs);
if (tree->GetBranch       ("tofhit"))
    tree->SetBranchAddress("tofhit", &fTofHit);
    tree->SetBranchAddress("trd",    &fTrd);
    tree->SetBranchAddress("trdk",   &fTrdK);
if (tree->GetBranch       ("trdg"))
    tree->SetBranchAddress("trdg",   &fTrdG);
if (tree->GetBranch       ("trdv"))
    tree->SetBranchAddress("trdv",   &fTrdV);
if (tree->GetBranch       ("trdhit"))
    tree->SetBranchAddress("trdhit", &fTrdHit);
    tree->SetBranchAddress("rich",   &fRich);
if (tree->GetBranch       ("tofhit"))
    tree->SetBranchAddress("ecal",   &fEcal);
  //tree->SetBranchAddress("ecalh",  &fEcalH);
if (tree->GetBranch       ("vertex"))
    tree->SetBranchAddress("vertex", &fVertex);
    tree->SetBranchAddress("mcinfo", &fMCinfo);
    tree->SetBranchAddress("mass",   &fMass);
if(tree->GetBranch       ("estim")){
    tree->SetBranchAddress("estim",  &fEstim);
}
else{cout<< "Warning: no estim branch in tree" << endl; }
if (tree->GetBranch       ("extra"))
    tree->SetBranchAddress("extra",  &fExtra);
if (tree->GetBranch       ("refit"))
    tree->SetBranchAddress("refit",  &fRefit);
if (tree->GetBranch       ("ecalx"))
    tree->SetBranchAddress("ecalx",  &fEcalX);
if(tree->GetBranch       ("ek")){
    tree->SetBranchAddress("ek",     &fEk);
  }



if(tree->GetBranch("richcut_flag")){
    tree->SetBranchAddress("richcut_flag", &fRichcut_flag);
  }
if(tree->GetBranch("transformer_output")){
    tree->SetBranchAddress("transformer_output", &fTransformer_output);
}
  }
};
