//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun  9 22:36:19 2011 by ROOT version 5.28/00d
// from TTree tree_MultiDie_CENT1/single
// found on file: Resultstakuv2c123456AnalysisResults_t.root
//////////////////////////////////////////////////////////

#ifndef ana_sgl_h
#define ana_sgl_h

#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TDirectory.h>
#include <vector>
#include <deque>
#include <cstdlib>
#include <TMath.h>

using namespace std;


class etrk;

struct etrk{ 
  etrk(double Cent=0, double PXv=0, double PYv=0, double PZv=0, 
       double Xv=0, double Yv=0, double Zv=0,
       double Px=0, double Py=0, double Pz=0, double Pt=0,
       double Eta=0, double Phi=0, double Theta=0, double Tpc=0, double Beta=0, 
       double E=0, double dPhi=0, double dEta=0,
       double Ntpc_ele1=0,double Ntpc_pio1=0,double Ntpc_kao1=0, double Ntpc_pro1=0,
       double Ntof_ele1=0,double Ntof_pio1=0,double Ntof_kao1=0, double Ntof_pro1=0,
       double Its=0, double Nits=0, double Ntpc=0, double Dcaxy=0, double Dcaz=0)
    :
    cent(Cent), pxv(PXv), pyv(PYv), pzv(PZv),
    xv(Xv), yv(Yv), zv(Zv),
    px(Px), py(Py), pz(Pz), pt(Pt),
    eta(Eta), phi(Phi), theta(Theta), tpc(Tpc), beta(Beta),
    e(E), dphi(dPhi), deta(dEta),
    ntpc_ele(Ntpc_ele1), ntpc_pio(Ntpc_pio1), ntpc_kao(Ntpc_kao1), ntpc_pro(Ntpc_pro1),
    ntof_ele(Ntof_ele1), ntof_pio(Ntof_pio1), ntof_kao(Ntof_kao1), ntof_pro(Ntof_pro1),
    its(Its), nits(Nits), ntpc(Ntpc), dcaxy(Dcaxy), dcaz(Dcaz)
  {
    conv_flag = 1;
    ghost_flag = 1;
  }

  double E(void){TMath::Sqrt(me*me+px*px+py*py+pz*pz);}
  double p(void){TMath::Sqrt(px*px+py*py+pz*pz);}

private:
  static const float me=0.000511; 
public:
  double cent;
  double pxv;
  double pyv;
  double pzv;
  double xv;
  double yv;
  double zv;
  double px;
  double py;
  double pz;
  double pt;
  double eta;
  double phi;
  double theta;
  double tpc;
  double beta;
  double e;
  double dphi;
  double deta;
  double ntpc_ele;
  double ntpc_pio;
  double ntpc_kao;
  double ntpc_pro;
  double ntof_ele;
  double ntof_pio;
  double ntof_kao;
  double ntof_pro;
  double its;
  double nits;
  double ntpc;
  double dcaxy;
  double dcaz;
  int conv_flag;
  int ghost_flag;
};


class ana_sgl {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           kNEvent;
   Double_t        kMag;
   TObjArray       *fkTriggerInfo;
   Double_t        kTriggerMask;
   Int_t           kTriggerCent;
   Double_t        fkNCut;
   Double_t        fkRunNumber;
   Double_t        fkCentrality;
   Double_t        fkXvPrim;
   Double_t        fkYvPrim;
   Double_t        fkZvPrim;
   Double_t        fkXRes;
   Double_t        fkYRes;
   Double_t        fkZRes;
   Double_t        fkNTrk;
   Double_t        fkTracks;
   Double_t        fkNacc;
   Double_t        fkNaccTrcklts;
   Double_t        fkNch;
   Double_t        fkZDCN1E;
   Double_t        fkZDCP1E;
   Double_t        fkZDCN2E;
   Double_t        fkZDCP2E;
   Double_t        fkV0A;
   Double_t        fkV0C;
   Int_t           fkNPar;
   Double_t        kPx[200];   //[fkNPar]
   Double_t        kPy[200];   //[fkNPar]
   Double_t        kPz[200];   //[fkNPar]
   Double_t        kPt[200];   //[fkNPar]
   Double_t        kP[200];   //[fkNPar]
   Double_t        kXv[200];   //[fkNPar]
   Double_t        kYv[200];   //[fkNPar]
   Double_t        kZv[200];   //[fkNPar]
   Double_t        kOneOverPt[200];   //[fkNPar]
   Double_t        kPhi[200];   //[fkNPar]
   Double_t        kTheta[200];   //[fkNPar]
   Double_t        kEta[200];   //[fkNPar]
   Double_t        kY[200];   //[fkNPar]
   Double_t        kE[200];   //[fkNPar]
   Double_t        kM[200];   //[fkNPar]
   Double_t        kCharge[200];   //[fkNPar]
   Double_t        kNclsITS[200];   //[fkNPar]
   Double_t        kNclsTPC[200];   //[fkNPar]
   Double_t        kNclsTPCiter1[200];   //[fkNPar]
   Double_t        kNFclsTPC[200];   //[fkNPar]
   Double_t        kNFclsTPCr[200];   //[fkNPar]
   Double_t        kNFclsTPCrFrac[200];   //[fkNPar]
   Double_t        kTPCsignalN[200];   //[fkNPar]
   Double_t        kTPCsignalNfrac[200];   //[fkNPar]
   Double_t        kTPCchi2Cl[200];   //[fkNPar]
   Double_t        kTrackStatus[200];   //[fkNPar]
   Double_t        kNclsTRD[200];   //[fkNPar]
   Double_t        kTRDntracklets[200];   //[fkNPar]
   Double_t        kTRDpidQuality[200];   //[fkNPar]
   Double_t        kTRDprobEle[200];   //[fkNPar]
   Double_t        kTRDprobPio[200];   //[fkNPar]
   Double_t        kImpactParXY[200];   //[fkNPar]
   Double_t        kImpactParZ[200];   //[fkNPar]
   Double_t        kTrackLength[200];   //[fkNPar]
   Double_t        kPdgCode[200];   //[fkNPar]
   Double_t        kPdgCodeMother[200];   //[fkNPar]
   Double_t        kPdgCodeGrandMother[200];   //[fkNPar]
   Double_t        kNumberOfDaughters[200];   //[fkNPar]
   Double_t        kHaveSameMother[200];   //[fkNPar]
   Double_t        kIsJpsiPrimary[200];   //[fkNPar]
   Double_t        kITSsignal[200];   //[fkNPar]
   Double_t        kITSsignalSSD1[200];   //[fkNPar]
   Double_t        kITSsignalSSD2[200];   //[fkNPar]
   Double_t        kITSsignalSDD1[200];   //[fkNPar]
   Double_t        kITSsignalSDD2[200];   //[fkNPar]
   Double_t        kITSclusterMap[200];   //[fkNPar]
   Double_t        kITSnSigmaEle[200];   //[fkNPar]
   Double_t        kITSnSigmaPio[200];   //[fkNPar]
   Double_t        kITSnSigmaMuo[200];   //[fkNPar]
   Double_t        kITSnSigmaKao[200];   //[fkNPar]
   Double_t        kITSnSigmaPro[200];   //[fkNPar]
   Double_t        kPIn[200];   //[fkNPar]
   Double_t        kTPCsignal[200];   //[fkNPar]
   Double_t        kTOFsignal[200];   //[fkNPar]
   Double_t        kTOFbeta[200];   //[fkNPar]
   Double_t        kTPCnSigmaEle[200];   //[fkNPar]
   Double_t        kTPCnSigmaPio[200];   //[fkNPar]
   Double_t        kTPCnSigmaMuo[200];   //[fkNPar]
   Double_t        kTPCnSigmaKao[200];   //[fkNPar]
   Double_t        kTPCnSigmaPro[200];   //[fkNPar]
   Double_t        kTOFnSigmaEle[200];   //[fkNPar]
   Double_t        kTOFnSigmaPio[200];   //[fkNPar]
   Double_t        kTOFnSigmaMuo[200];   //[fkNPar]
   Double_t        kTOFnSigmaKao[200];   //[fkNPar]
   Double_t        kTOFnSigmaPro[200];   //[fkNPar]
   Double_t        kKinkIndex0[200];   //[fkNPar]
   Double_t        kChi2NDF[200];   //[fkNPar]
   Double_t        kDecayLength[200];   //[fkNPar]
   Double_t        kR[200];   //[fkNPar]
   Double_t        kOpeningAngle[200];   //[fkNPar]
   Double_t        kThetaHE[200];   //[fkNPar]
   Double_t        kPhiHE[200];   //[fkNPar]
   Double_t        kThetaCS[200];   //[fkNPar]
   Double_t        kPhiCS[200];   //[fkNPar]
   Double_t        kLegDist[200];   //[fkNPar]
   Double_t        kLegDistXY[200];   //[fkNPar]
   Double_t        kDeltaEta[200];   //[fkNPar]
   Double_t        kDeltaPhi[200];   //[fkNPar]
   Double_t        kMerr[200];   //[fkNPar]
   Double_t        kDCA[200];   //[fkNPar]
   Double_t        kPairType[200];   //[fkNPar]
   Double_t        kPseudoProperTime[200];   //[fkNPar]
   Double_t        kXvPrim[200];   //[fkNPar]
   Double_t        kYvPrim[200];   //[fkNPar]
   Double_t        kZvPrim[200];   //[fkNPar]
   Double_t        kXRes[200];   //[fkNPar]
   Double_t        kYRes[200];   //[fkNPar]
   Double_t        kZRes[200];   //[fkNPar]
   Double_t        kNTrk[200];   //[fkNPar]
   Double_t        kTracks[200];   //[fkNPar]
   Double_t        kNacc[200];   //[fkNPar]
   Double_t        kNaccTrcklts[200];   //[fkNPar]
   Double_t        kNch[200];   //[fkNPar]
   Double_t        kCentrality[200];   //[fkNPar]
   Double_t        kNevents[200];   //[fkNPar]

   // List of branches
   TBranch        *b_kNEvent;   //!
   TBranch        *b_kMag;   //!
   TBranch        *b_fkTriggerInfo;   //!
   TBranch        *b_kTriggerMask;   //!
   TBranch        *b_kTriggerCent;   //!
   TBranch        *b_fkNCut;   //!
   TBranch        *b_fkRunNumber;   //!
   TBranch        *b_fkCentrality;   //!
   TBranch        *b_fkXvPrim;   //!
   TBranch        *b_fkYvPrim;   //!
   TBranch        *b_fkZvPrim;   //!
   TBranch        *b_fkXRes;   //!
   TBranch        *b_fkYRes;   //!
   TBranch        *b_fkZRes;   //!
   TBranch        *b_fkNTrk;   //!
   TBranch        *b_fkTracks;   //!
   TBranch        *b_fkNacc;   //!
   TBranch        *b_fkNaccTrcklts;   //!
   TBranch        *b_fkNch;   //!
   TBranch        *b_fkZDCN1E;   //!
   TBranch        *b_fkZDCP1E;   //!
   TBranch        *b_fkZDCN2E;   //!
   TBranch        *b_fkZDCP2E;   //!
   TBranch        *b_fkV0A;   //!
   TBranch        *b_fkV0C;   //!
   TBranch        *b_fkNPar;   //!
   TBranch        *b_kPx;   //!
   TBranch        *b_kPy;   //!
   TBranch        *b_kPz;   //!
   TBranch        *b_kPt;   //!
   TBranch        *b_kP;   //!
   TBranch        *b_kXv;   //!
   TBranch        *b_kYv;   //!
   TBranch        *b_kZv;   //!
   TBranch        *b_kOneOverPt;   //!
   TBranch        *b_kPhi;   //!
   TBranch        *b_kTheta;   //!
   TBranch        *b_kEta;   //!
   TBranch        *b_kY;   //!
   TBranch        *b_kE;   //!
   TBranch        *b_kM;   //!
   TBranch        *b_kCharge;   //!
   TBranch        *b_kNclsITS;   //!
   TBranch        *b_kNclsTPC;   //!
   TBranch        *b_kNclsTPCiter1;   //!
   TBranch        *b_kNFclsTPC;   //!
   TBranch        *b_kNFclsTPCr;   //!
   TBranch        *b_kNFclsTPCrFrac;   //!
   TBranch        *b_kTPCsignalN;   //!
   TBranch        *b_kTPCsignalNfrac;   //!
   TBranch        *b_kTPCchi2Cl;   //!
   TBranch        *b_kTrackStatus;   //!
   TBranch        *b_kNclsTRD;   //!
   TBranch        *b_kTRDntracklets;   //!
   TBranch        *b_kTRDpidQuality;   //!
   TBranch        *b_kTRDprobEle;   //!
   TBranch        *b_kTRDprobPio;   //!
   TBranch        *b_kImpactParXY;   //!
   TBranch        *b_kImpactParZ;   //!
   TBranch        *b_kTrackLength;   //!
   TBranch        *b_kPdgCode;   //!
   TBranch        *b_kPdgCodeMother;   //!
   TBranch        *b_kPdgCodeGrandMother;   //!
   TBranch        *b_kNumberOfDaughters;   //!
   TBranch        *b_kHaveSameMother;   //!
   TBranch        *b_kIsJpsiPrimary;   //!
   TBranch        *b_kITSsignal;   //!
   TBranch        *b_kITSsignalSSD1;   //!
   TBranch        *b_kITSsignalSSD2;   //!
   TBranch        *b_kITSsignalSDD1;   //!
   TBranch        *b_kITSsignalSDD2;   //!
   TBranch        *b_kITSclusterMap;   //!
   TBranch        *b_kITSnSigmaEle;   //!
   TBranch        *b_kITSnSigmaPio;   //!
   TBranch        *b_kITSnSigmaMuo;   //!
   TBranch        *b_kITSnSigmaKao;   //!
   TBranch        *b_kITSnSigmaPro;   //!
   TBranch        *b_kPIn;   //!
   TBranch        *b_kTPCsignal;   //!
   TBranch        *b_kTOFsignal;   //!
   TBranch        *b_kTOFbeta;   //!
   TBranch        *b_kTPCnSigmaEle;   //!
   TBranch        *b_kTPCnSigmaPio;   //!
   TBranch        *b_kTPCnSigmaMuo;   //!
   TBranch        *b_kTPCnSigmaKao;   //!
   TBranch        *b_kTPCnSigmaPro;   //!
   TBranch        *b_kTOFnSigmaEle;   //!
   TBranch        *b_kTOFnSigmaPio;   //!
   TBranch        *b_kTOFnSigmaMuo;   //!
   TBranch        *b_kTOFnSigmaKao;   //!
   TBranch        *b_kTOFnSigmaPro;   //!
   TBranch        *b_kKinkIndex0;   //!
   TBranch        *b_kChi2NDF;   //!
   TBranch        *b_kDecayLength;   //!
   TBranch        *b_kR;   //!
   TBranch        *b_kOpeningAngle;   //!
   TBranch        *b_kThetaHE;   //!
   TBranch        *b_kPhiHE;   //!
   TBranch        *b_kThetaCS;   //!
   TBranch        *b_kPhiCS;   //!
   TBranch        *b_kLegDist;   //!
   TBranch        *b_kLegDistXY;   //!
   TBranch        *b_kDeltaEta;   //!
   TBranch        *b_kDeltaPhi;   //!
   TBranch        *b_kMerr;   //!
   TBranch        *b_kDCA;   //!
   TBranch        *b_kPairType;   //!
   TBranch        *b_kPseudoProperTime;   //!
   TBranch        *b_kXvPrim;   //!
   TBranch        *b_kYvPrim;   //!
   TBranch        *b_kZvPrim;   //!
   TBranch        *b_kXRes;   //!
   TBranch        *b_kYRes;   //!
   TBranch        *b_kZRes;   //!
   TBranch        *b_kNTrk;   //!
   TBranch        *b_kTracks;   //!
   TBranch        *b_kNacc;   //!
   TBranch        *b_kNaccTrcklts;   //!
   TBranch        *b_kNch;   //!
   TBranch        *b_kCentrality;   //!
   TBranch        *b_kNevents;   //!

   ana_sgl(TTree *tree=0);
   virtual ~ana_sgl();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);


   //////////my own function/////////////////

   TH2D *hdedx_pt;
   TH2D *hdedx_tof_elec_pt;
   TH2D *hdedx_tof_all_pt;
   TH2D *hdedx_tof_elec_emc_pt;
   TH2D *hdedx_tof_all_emc_pt;
   TH2D *hdedx_emc_pt;

   TH2D *hbetatof_pt;
   TH2D *hbetatof_tof_elec_pt;
   TH2D *hbetatof_tof_all_pt;
   TH2D *hbetatof_tof_elec_emc_pt;
   TH2D *hbetatof_tof_all_emc_pt;
   TH2D *hbetatof_emc_pt;

   TH1D *hdedx[1000]; 
   TH1D *hdedx_tof_elec[1000]; 
   TH1D *hdedx_tof_all[1000]; 
   TH1D *hdedx_tof_elec_emc[1000]; 
   TH1D *hdedx_tof_all_emc[1000]; 
   TH1D *hdedx_emc[1000]; 


   TH1D *fEventStat;                  //! Histogram with event statistics
   TH1D *fEvent;
   TH2D *fdEdXvsPt;
   TH2D *fdEdXnSigmaElecvsPt;
   TH2D *fTOFbetavsPt;
   TH2D *fTOFnSigmaElecvsPt;
   TH1F *hCentrality;
   TH2F *hV0AC;
   TH2F *hV0AC_Ntrk;
   TH2F *hV0AC_NaccTrcklts;



   TH2D *fdEdXvsPtWithCut;
   TH1D *fPtElec[2];
   
   TH2D *fNelc_pos_cent[10]; 
   TH1D *fNelc_all_cent[10];
   
   TH2D *fNelc_pos; 
   TH1D *fNelc_all;
   TH2D *fNelc_all_pT;


   //// pair tree
   TH2D *hmasspt[7][11];
   TH2D *hmasspt_weight[7][11];

   TFile *fout; 
   Float_t fBinWidth ;
   Int_t nHistos ;

   void ana_init(char *filename);
   void ana_event(int ientry, int jentry);
   void ana_end(void);
   void loop_a_file(char *filename);
   void ana_set_simflag(bool a) { simflag = a;}

   bool kTOFcut(int itrk);
   bool GlobalTrackcut(int itrk);
   //bool PairTrackcut(int itrk);
   void fill_histograms(int itrk);
   void fill_to_tree_variables(void);
   void fill_to_tree_track_variables(int itrk);
   void add_histograms(TFile *fin);
   void select_trigger(int trig){ sel_trigger = trig ; }


   //// cut function
   void set_tpc_dedx_cuts(double low, double high){
     d_tpc_dedx_low = low;
     d_tpc_dedx_high = high;
   };
   void set_tof_cuts(double low, double high){
     d_flag_tof_cut = true;
     d_tof_low = low;
     d_tof_high = high;
   };
   
   void set_veto_for_kaon(double low, double high){
     d_flag_kaon_veto = true;
     d_dedx_kaon_veto_low = low;
     d_dedx_kaon_veto_high = high;
   }

   void set_veto_for_proton(double low, double high){
     d_flag_proton_veto = true;
     d_dedx_proton_veto_low = low;
     d_dedx_proton_veto_high = high;
   }

   void enable_pair_emc_cut(double low, double high){
     d_flag_emc_cut = true;
     d_emc_low = low;
     d_emc_high = high;
   };

   void enable_pair_phiv_cut(double low){
     d_flag_phiv = true;
     d_phiv_cut = low;
   };

   void enable_pait_pt_cut(double low, double high){
     d_flag_pt_cut = true;
     d_pt_cut_low = low;
     d_pt_cut_high = high;
   };

   void print_cuts(void);

   //// private functions 

   void calc_pair(vector<etrk> e1, vector<etrk> e2);
   void randomize_pool(vector<etrk> e1, vector<etrk> e2);
   //void fill_pair(etrk *e1, etrk *e2, int type);
   void fill_pair(vector<etrk>::iterator e1, vector<etrk>::iterator e2, int type);
   bool PairTrackcut(vector<etrk>::iterator e1, vector<etrk>::iterator e2);
   bool reject_conversion(bool val){d_conv_flag = val;}
   void check_conversion_pairs(vector<etrk> &e1, vector<etrk> &e2);
   void check_ghost_pairs(vector<etrk> &e1);
   void calc_vars(vector<etrk>::iterator e1, vector<etrk>::iterator e2, double &mass, double &phiv, double &px, double &py, double&pz,
		  double &pt, double &e, double &phi, double &eta, double &cos, double &psi);
   bool pair_cut(void);

   TTree *d_tree;
   TTree *d_ntpair;
   bool  simflag;
   int   sel_trigger;

 private:
   bool d_flag_tof_cut;
   bool d_flag_emc_cut;
   bool d_flag_phiv;
   bool d_flag_pt_cut;
   bool d_flag_kaon_veto;
   bool d_flag_proton_veto;
   double d_tpc_dedx_low;
   double d_tpc_dedx_high;
   double d_tof_low;
   double d_tof_high;
   double d_emc_low;
   double d_emc_high;
   double d_phiv_cut;
   double d_pt_cut_low ;
   double d_pt_cut_high ;
   double d_dedx_kaon_veto_low;   
   double d_dedx_kaon_veto_high;
   double d_dedx_proton_veto_low;   
   double d_dedx_proton_veto_high;

   bool magnetic_field_mm ;
   bool d_conv_flag;
   Int_t d_evt;
   Double_t d_cent;
   Double_t d_ntrk;
   Double_t d_xvprim;
   Double_t d_yvprim;
   Double_t d_zvprim;
   Double_t d_nacctrklets;
   Double_t d_xres;
   Double_t d_yres;
   Double_t d_zres;
   Int_t d_nelc;
   Double_t d_px[100];
   Double_t d_py[100];
   Double_t d_pz[100];
   Double_t d_p[100];
   Double_t d_pt[100];
   Double_t d_xv[100];
   Double_t d_yv[100];
   Double_t d_zv[100];
   Double_t d_phi[100];
   Double_t d_theta[100];
   Double_t d_eta[100];
   Double_t d_c[100];
   Double_t d_nclusITS[100];
   Double_t d_nclusTPC[100];
   Double_t d_nclusTPCiter[100];
   Double_t d_nfclusTPC[100];
   Double_t d_nfclusTPCr[100];
   Double_t d_nfclusTPCrFrac[100];
   Double_t d_TPCsignalN[100];
   Double_t d_TPCsignalNfrac[100];
   Double_t d_TPCchi2cl[100];
   Double_t d_trkstat[100];
   Double_t d_nclsTRD[100];
   Double_t d_TRDntracklets[100];
   Double_t d_TRDpidquality[100];
   Double_t d_TRDprobEle[100];
   Double_t d_TRDprobPio[100];
   Double_t d_impactXY[100];
   Double_t d_impactZ[100];
   Double_t d_tracklength[100];
   Double_t d_ITSsignal[100];
   Double_t d_ITSnsigmaEle[100];
   Double_t d_ITSnsigmaPio[100];
   Double_t d_ITSnsigmaMuo[100];
   Double_t d_ITSnsigmaKao[100];
   Double_t d_ITSnsigmaPro[100];
   Double_t d_PIn[100];
   Double_t d_TPCsignal[100];
   Double_t d_TOFsignal[100];
   Double_t d_TOFbeta[100];
   Double_t d_TPCnSigmaEle[100];
   Double_t d_TPCnSigmaPio[100];
   Double_t d_TPCnSigmaMuo[100];
   Double_t d_TPCnSigmaKao[100];
   Double_t d_TPCnSigmaPro[100];
   Double_t d_TOFnSigmaEle[100];
   Double_t d_TOFnSigmaPio[100];
   Double_t d_TOFnSigmaMuo[100];
   Double_t d_TOFnSigmaKao[100];
   Double_t d_TOFnSigmaPro[100];

   Double_t d_chi2ndf[100];
   Double_t d_E[100];
   Double_t d_dphi[100];
   Double_t d_deta[100];

   //////// pair variable 
   Double_t d_run;
   Int_t d_event;
   Double_t d_centrality;
   Double_t d_prim_xv;
   Double_t d_prim_yv;
   Double_t d_prim_zv;
   Double_t d_mass;
   Double_t d_pxpair;
   Double_t d_pypair;
   Double_t d_pzpair;
   Double_t d_ptpair;
   Double_t d_epair;
   Double_t d_etapair;
   Double_t d_phipair;
   Double_t d_cos;
   Double_t d_phiv;
   Double_t d_psi;
   Int_t d_pairtype;

   /////////////////////////
   Double_t d_cent1;
   Double_t d_xv1;
   Double_t d_yv1;
   Double_t d_zv1;
   Double_t d_px1;
   Double_t d_py1;
   Double_t d_pz1;
   Double_t d_pt1;
   Double_t d_eta1;
   Double_t d_phi1;
   Double_t d_theta1;
   Double_t d_tpc1;
   Double_t d_ntpc_ele1;
   Double_t d_ntpc_pio1;
   Double_t d_ntpc_kao1;
   Double_t d_ntpc_pro1;
   Double_t d_beta1;
   Double_t d_ntof_ele1;
   Double_t d_ntof_pio1;
   Double_t d_ntof_kao1;
   Double_t d_ntof_pro1;
   Double_t d_its1;
   Double_t d_nits1;
   Double_t d_ntpc1;
   Double_t d_e1;
   Double_t d_dphi1;
   Double_t d_deta1;
   Double_t d_dcaxy1;
   Double_t d_dcaz1;
   Int_t d_conv1;

   Double_t d_cent2;
   Double_t d_xv2;
   Double_t d_yv2;
   Double_t d_zv2;
   Double_t d_px2;
   Double_t d_py2;
   Double_t d_pz2;
   Double_t d_pt2;
   Double_t d_eta2;
   Double_t d_phi2;
   Double_t d_theta2;
   Double_t d_tpc2;
   Double_t d_ntpc_ele2;
   Double_t d_ntpc_pio2;
   Double_t d_ntpc_kao2;
   Double_t d_ntpc_pro2;
   Double_t d_beta2;
   Double_t d_ntof_ele2;
   Double_t d_ntof_pio2;
   Double_t d_ntof_kao2;
   Double_t d_ntof_pro2;
   Double_t d_its2;
   Double_t d_nits2;
   Double_t d_ntpc2;
   Double_t d_e2;
   Double_t d_dphi2;
   Double_t d_deta2;
   Double_t d_dcaxy2;
   Double_t d_dcaz2;
    Int_t d_conv2;

   Int_t nelec_pos[10]; //pT bin;

   vector<etrk> vem;
   vector<etrk> vep;
   vector<etrk> vem_tmp;
   vector<etrk> vep_tmp;


};

#endif

#ifdef ana_sgl_cxx
ana_sgl::ana_sgl(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Resultstakuv2c123456AnalysisResults_t.root");
      if (!f) {
         f = new TFile("Resultstakuv2c123456AnalysisResults_t.root");
         f->cd("Resultstakuv2c123456AnalysisResults_t.root:/PWG3_dielectron");
      }
      tree = (TTree*)gDirectory->Get("tree_MultiDie_CENT1");

   }
   Init(tree);
   d_conv_flag = false;
}

ana_sgl::~ana_sgl()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ana_sgl::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ana_sgl::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ana_sgl::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("kNEvent", &kNEvent, &b_kNEvent);
   fChain->SetBranchAddress("kMag", &kMag, &b_kMag);
   fChain->SetBranchAddress("fkTriggerInfo", &fkTriggerInfo, &b_fkTriggerInfo);
   fChain->SetBranchAddress("kTriggerMask", &kTriggerMask, &b_kTriggerMask);
   fChain->SetBranchAddress("kTriggerCent", &kTriggerCent, &b_kTriggerCent);
   fChain->SetBranchAddress("fkNCut", &fkNCut, &b_fkNCut);
   fChain->SetBranchAddress("fkRunNumber", &fkRunNumber, &b_fkRunNumber);
   fChain->SetBranchAddress("fkCentrality", &fkCentrality, &b_fkCentrality);
   fChain->SetBranchAddress("fkXvPrim", &fkXvPrim, &b_fkXvPrim);
   fChain->SetBranchAddress("fkYvPrim", &fkYvPrim, &b_fkYvPrim);
   fChain->SetBranchAddress("fkZvPrim", &fkZvPrim, &b_fkZvPrim);
   fChain->SetBranchAddress("fkXRes", &fkXRes, &b_fkXRes);
   fChain->SetBranchAddress("fkYRes", &fkYRes, &b_fkYRes);
   fChain->SetBranchAddress("fkZRes", &fkZRes, &b_fkZRes);
   fChain->SetBranchAddress("fkNTrk", &fkNTrk, &b_fkNTrk);
   fChain->SetBranchAddress("fkTracks", &fkTracks, &b_fkTracks);
   fChain->SetBranchAddress("fkNacc", &fkNacc, &b_fkNacc);
   fChain->SetBranchAddress("fkNaccTrcklts", &fkNaccTrcklts, &b_fkNaccTrcklts);
   fChain->SetBranchAddress("fkNch", &fkNch, &b_fkNch);
   fChain->SetBranchAddress("fkZDCN1E", &fkZDCN1E, &b_fkZDCN1E);
   fChain->SetBranchAddress("fkZDCP1E", &fkZDCP1E, &b_fkZDCP1E);
   fChain->SetBranchAddress("fkZDCN2E", &fkZDCN2E, &b_fkZDCN2E);
   fChain->SetBranchAddress("fkZDCP2E", &fkZDCP2E, &b_fkZDCP2E);
   fChain->SetBranchAddress("fkV0A", &fkV0A, &b_fkV0A);
   fChain->SetBranchAddress("fkV0C", &fkV0C, &b_fkV0C);
   fChain->SetBranchAddress("fkNPar", &fkNPar, &b_fkNPar);
   fChain->SetBranchAddress("kPx", kPx, &b_kPx);
   fChain->SetBranchAddress("kPy", kPy, &b_kPy);
   fChain->SetBranchAddress("kPz", kPz, &b_kPz);
   fChain->SetBranchAddress("kPt", kPt, &b_kPt);
   fChain->SetBranchAddress("kP", kP, &b_kP);
   fChain->SetBranchAddress("kXv", kXv, &b_kXv);
   fChain->SetBranchAddress("kYv", kYv, &b_kYv);
   fChain->SetBranchAddress("kZv", kZv, &b_kZv);
   fChain->SetBranchAddress("kOneOverPt", kOneOverPt, &b_kOneOverPt);
   fChain->SetBranchAddress("kPhi", kPhi, &b_kPhi);
   fChain->SetBranchAddress("kTheta", kTheta, &b_kTheta);
   fChain->SetBranchAddress("kEta", kEta, &b_kEta);
   fChain->SetBranchAddress("kY", kY, &b_kY);
   fChain->SetBranchAddress("kE", kE, &b_kE);
   fChain->SetBranchAddress("kM", kM, &b_kM);
   fChain->SetBranchAddress("kCharge", kCharge, &b_kCharge);
   fChain->SetBranchAddress("kNclsITS", kNclsITS, &b_kNclsITS);
   fChain->SetBranchAddress("kNclsTPC", kNclsTPC, &b_kNclsTPC);
   fChain->SetBranchAddress("kNclsTPCiter1", kNclsTPCiter1, &b_kNclsTPCiter1);
   fChain->SetBranchAddress("kNFclsTPC", kNFclsTPC, &b_kNFclsTPC);
   fChain->SetBranchAddress("kNFclsTPCr", kNFclsTPCr, &b_kNFclsTPCr);
   fChain->SetBranchAddress("kNFclsTPCrFrac", kNFclsTPCrFrac, &b_kNFclsTPCrFrac);
   fChain->SetBranchAddress("kTPCsignalN", kTPCsignalN, &b_kTPCsignalN);
   fChain->SetBranchAddress("kTPCsignalNfrac", kTPCsignalNfrac, &b_kTPCsignalNfrac);
   fChain->SetBranchAddress("kTPCchi2Cl", kTPCchi2Cl, &b_kTPCchi2Cl);
   fChain->SetBranchAddress("kTrackStatus", kTrackStatus, &b_kTrackStatus);
   fChain->SetBranchAddress("kNclsTRD", kNclsTRD, &b_kNclsTRD);
   fChain->SetBranchAddress("kTRDntracklets", kTRDntracklets, &b_kTRDntracklets);
   fChain->SetBranchAddress("kTRDpidQuality", kTRDpidQuality, &b_kTRDpidQuality);
   fChain->SetBranchAddress("kTRDprobEle", kTRDprobEle, &b_kTRDprobEle);
   fChain->SetBranchAddress("kTRDprobPio", kTRDprobPio, &b_kTRDprobPio);
   fChain->SetBranchAddress("kImpactParXY", kImpactParXY, &b_kImpactParXY);
   fChain->SetBranchAddress("kImpactParZ", kImpactParZ, &b_kImpactParZ);
   fChain->SetBranchAddress("kTrackLength", kTrackLength, &b_kTrackLength);
   fChain->SetBranchAddress("kPdgCode", kPdgCode, &b_kPdgCode);
   fChain->SetBranchAddress("kPdgCodeMother", kPdgCodeMother, &b_kPdgCodeMother);
   fChain->SetBranchAddress("kPdgCodeGrandMother", kPdgCodeGrandMother, &b_kPdgCodeGrandMother);
   fChain->SetBranchAddress("kNumberOfDaughters", kNumberOfDaughters, &b_kNumberOfDaughters);
   fChain->SetBranchAddress("kHaveSameMother", kHaveSameMother, &b_kHaveSameMother);
   fChain->SetBranchAddress("kIsJpsiPrimary", kIsJpsiPrimary, &b_kIsJpsiPrimary);
   fChain->SetBranchAddress("kITSsignal", kITSsignal, &b_kITSsignal);
   fChain->SetBranchAddress("kITSsignalSSD1", kITSsignalSSD1, &b_kITSsignalSSD1);
   fChain->SetBranchAddress("kITSsignalSSD2", kITSsignalSSD2, &b_kITSsignalSSD2);
   fChain->SetBranchAddress("kITSsignalSDD1", kITSsignalSDD1, &b_kITSsignalSDD1);
   fChain->SetBranchAddress("kITSsignalSDD2", kITSsignalSDD2, &b_kITSsignalSDD2);
   fChain->SetBranchAddress("kITSclusterMap", kITSclusterMap, &b_kITSclusterMap);
   fChain->SetBranchAddress("kITSnSigmaEle", kITSnSigmaEle, &b_kITSnSigmaEle);
   fChain->SetBranchAddress("kITSnSigmaPio", kITSnSigmaPio, &b_kITSnSigmaPio);
   fChain->SetBranchAddress("kITSnSigmaMuo", kITSnSigmaMuo, &b_kITSnSigmaMuo);
   fChain->SetBranchAddress("kITSnSigmaKao", kITSnSigmaKao, &b_kITSnSigmaKao);
   fChain->SetBranchAddress("kITSnSigmaPro", kITSnSigmaPro, &b_kITSnSigmaPro);
   fChain->SetBranchAddress("kPIn", kPIn, &b_kPIn);
   fChain->SetBranchAddress("kTPCsignal", kTPCsignal, &b_kTPCsignal);
   fChain->SetBranchAddress("kTOFsignal", kTOFsignal, &b_kTOFsignal);
   fChain->SetBranchAddress("kTOFbeta", kTOFbeta, &b_kTOFbeta);
   fChain->SetBranchAddress("kTPCnSigmaEle", kTPCnSigmaEle, &b_kTPCnSigmaEle);
   fChain->SetBranchAddress("kTPCnSigmaPio", kTPCnSigmaPio, &b_kTPCnSigmaPio);
   fChain->SetBranchAddress("kTPCnSigmaMuo", kTPCnSigmaMuo, &b_kTPCnSigmaMuo);
   fChain->SetBranchAddress("kTPCnSigmaKao", kTPCnSigmaKao, &b_kTPCnSigmaKao);
   fChain->SetBranchAddress("kTPCnSigmaPro", kTPCnSigmaPro, &b_kTPCnSigmaPro);
   fChain->SetBranchAddress("kTOFnSigmaEle", kTOFnSigmaEle, &b_kTOFnSigmaEle);
   fChain->SetBranchAddress("kTOFnSigmaPio", kTOFnSigmaPio, &b_kTOFnSigmaPio);
   fChain->SetBranchAddress("kTOFnSigmaMuo", kTOFnSigmaMuo, &b_kTOFnSigmaMuo);
   fChain->SetBranchAddress("kTOFnSigmaKao", kTOFnSigmaKao, &b_kTOFnSigmaKao);
   fChain->SetBranchAddress("kTOFnSigmaPro", kTOFnSigmaPro, &b_kTOFnSigmaPro);
   fChain->SetBranchAddress("kKinkIndex0", kKinkIndex0, &b_kKinkIndex0);
   fChain->SetBranchAddress("kChi2NDF", kChi2NDF, &b_kChi2NDF);
   fChain->SetBranchAddress("kDecayLength", kDecayLength, &b_kDecayLength);
   fChain->SetBranchAddress("kR", kR, &b_kR);
   fChain->SetBranchAddress("kOpeningAngle", kOpeningAngle, &b_kOpeningAngle);
   fChain->SetBranchAddress("kThetaHE", kThetaHE, &b_kThetaHE);
   fChain->SetBranchAddress("kPhiHE", kPhiHE, &b_kPhiHE);
   fChain->SetBranchAddress("kThetaCS", kThetaCS, &b_kThetaCS);
   fChain->SetBranchAddress("kPhiCS", kPhiCS, &b_kPhiCS);
   fChain->SetBranchAddress("kLegDist", kLegDist, &b_kLegDist);
   fChain->SetBranchAddress("kLegDistXY", kLegDistXY, &b_kLegDistXY);
   fChain->SetBranchAddress("kDeltaEta", kDeltaEta, &b_kDeltaEta);
   fChain->SetBranchAddress("kDeltaPhi", kDeltaPhi, &b_kDeltaPhi);
   fChain->SetBranchAddress("kMerr", kMerr, &b_kMerr);
   fChain->SetBranchAddress("kDCA", kDCA, &b_kDCA);
   fChain->SetBranchAddress("kPairType", kPairType, &b_kPairType);
   fChain->SetBranchAddress("kPseudoProperTime", kPseudoProperTime, &b_kPseudoProperTime);
   fChain->SetBranchAddress("kXvPrim", kXvPrim, &b_kXvPrim);
   fChain->SetBranchAddress("kYvPrim", kYvPrim, &b_kYvPrim);
   fChain->SetBranchAddress("kZvPrim", kZvPrim, &b_kZvPrim);
   fChain->SetBranchAddress("kXRes", kXRes, &b_kXRes);
   fChain->SetBranchAddress("kYRes", kYRes, &b_kYRes);
   fChain->SetBranchAddress("kZRes", kZRes, &b_kZRes);
   fChain->SetBranchAddress("kNTrk", kNTrk, &b_kNTrk);
   fChain->SetBranchAddress("kTracks", kTracks, &b_kTracks);
   fChain->SetBranchAddress("kNacc", kNacc, &b_kNacc);
   fChain->SetBranchAddress("kNaccTrcklts", kNaccTrcklts, &b_kNaccTrcklts);
   fChain->SetBranchAddress("kNch", kNch, &b_kNch);
   fChain->SetBranchAddress("kCentrality", kCentrality, &b_kCentrality);
   fChain->SetBranchAddress("kNevents", kNevents, &b_kNevents);
   Notify();
}

Bool_t ana_sgl::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ana_sgl::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ana_sgl::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ana_sgl_cxx
