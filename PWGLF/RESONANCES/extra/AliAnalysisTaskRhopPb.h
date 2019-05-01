#ifndef AliAnalysisTaskRhopPb_cxx
#define AliAnalysisTaskRhopPb_cxx

#include <AliPID.h>
#include <AliPIDResponse.h>
#include <AliPIDCombined.h>
#include "AliEventCuts.h"

// example of an analysis task creating trees for rho in pPb analysis

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH3F;
class AliESDtrackCuts;
class AliESDEvent ;
class fPIDResponse;
class fPIDCombined;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskRhopPb : public AliAnalysisTaskSE {

public:
  AliAnalysisTaskRhopPb(const char *name = "AliAnalysisTaskRhopPb");
  virtual ~AliAnalysisTaskRhopPb();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Bool_t UserNotify();

  AliEventCuts fEventCuts; /// Event cuts

protected:  
  TH1I* hStats; //!
  TH2F* PIDQA; //!
  TH2F* phitr_all; //!
  TH2F* phitr_tpc; //!
  TH2F* phitr_tof; //!
  TH2F* phitr_pid; //!
  TH2F* cosdip; //!

private:
  AliAnalysisTaskRhopPb(const AliAnalysisTaskRhopPb&); // not implemented
  AliAnalysisTaskRhopPb& operator=(const AliAnalysisTaskRhopPb&); // not implemented
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram with name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram with name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram with name key
  Bool_t TestLambda(Double_t pt,Double_t l1,Double_t l2) ;
  Double_t InvMass( Double_t m1, Double_t* p1, Double_t m2, Double_t* p2 );
  Double_t InvMass3( Double_t m1, Double_t* p1, Double_t m2, Double_t* p2, Double_t m3, Double_t* p3 );
  Float_t  TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge);

private:
  AliESDtrackCuts *fESDtrackCuts; // Track cut
  TList * fOutputContainer;       //final histogram container
  TTree*  fMyTree; //!
 
  Int_t  fEventCounter = 0, fEventCounter2 = 0, fEventCounter2_1 = 0, fEventCounter2_2 = 0;  // number of analyzed events
  Int_t  fnCINT1B = 0, fnCINT1A = 0, fnCINT1C = 0,fnCINT1E = 0;
  Int_t  fEventCounter1 = 0;         // number of analyzed events
  Int_t  fnCINT1B1 = 0, fnCINT1A1 = 0, fnCINT1C1 = 0,fnCINT1E1 = 0;
  
  Float_t  centr = 0., centr2 = 0.;
  Int_t    rnumb = 0, trig_fl = 0, trig_mc = 0, trig_sc = 0, trig_mb = 0, cbin = 0;
  Int_t    trMult = 0;
  Int_t    ntr = 0, ntr_all = 0, ntr_esdall = 0, ng = 0, nMult = 0;
  float    vtx51[3];

  Int_t VtxNCTrack = 0, VtxNCSPD = 0, VtxFlag = 0;
  Int_t VtxNCTrack_def = 0, VtxNCSPD_def = 0;

  Float_t chi2PerClusterITS = 0.;
  Float_t chi2PerClusterTPC = 0.;
  Int_t   nClustersITS = 0;
  Int_t   nClustersTPC = 0;
  Int_t   nClustersTPCShared = 0;
  Float_t fracClustersTPCShared = 0.;
  Float_t nCrossedRowsTPC = 0.;
  Float_t ratioCrossedRowsOverFindableClustersTPC = 0.;
 
  Float_t chi2TPC[5000], chi2ITS[5000], fracTPCS[5000], nTPCC[5000], ratTPCCF[5000];
  Int_t   nclITS[5000], nclTPC[5000], nclTPCS[5000];

  Int_t    nCellG[5000], nmodG[5000], nxG[5000], nzG[5000];
  Float_t  dxG[5000], dzG[5000], enG[5000], momG[5000][3], M20G[5000], M02G[5000]; 
  //Int_t    decay;
  Float_t  xyImp_pip = 0., zImp_pip = 0., xyImp_pim = 0., zImp_pim = 0.;
  Int_t    nClTr[5000], chTr[5000], nITSTr[5000], nSPDTr[5000], dchG[5000], TPCref[5000], ITSref[5000], TOF_match = 0;
  Float_t  nChiTr[5000], xyzImpTr[5000][2], momTr[5000][3], cpvG[5000], dmomG[5000];
  //Int_t    disp1; 
  Float_t  dispG[5000]; //, pt_pim, mPM;
  Float_t  pp1 = 0., pp2 = 0., pt1 = 0., pt2 = 0., pz1 = 0., pz2 = 0., theta = 0., ppp = 0., phi_pip = 0., phi = 0.;
  
  // PID
  Int_t    skip_PID = 0;
  Float_t  pid_sigma[5000][6], pid_sigma1[5000], pid_sigma2[5000], pid_sigma3[5000];  // 0->TPC-pi, 1->TPC-k, 2->TPC-p, 3->TOF-pi, 4->TOF-k, 4->TOF-p
  Float_t  pid_comb[5000][5];
  // *************************

  const AliPIDResponse *fPIDResponse;     //! PID response object
  AliPIDCombined       *fPIDCombined;     //! combined PID object
  

  ClassDef(AliAnalysisTaskRhopPb, 1);
};

#endif


