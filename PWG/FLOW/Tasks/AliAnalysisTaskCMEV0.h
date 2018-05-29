/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskCMEV0_H
#define AliAnalysisTaskCMEV0_H

/////////////////////////////////////////////////
// AliAnalysisTaskCMEV0:
// analysis task for ZDC gain Equalization
// and CME analysis using VZERO Detector
// Author: Rihan Haque (mhaque@cern.ch)
/////////////////////////////////////////////////

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TString.h"
//#include "TRandom3.h"
#include "AliAnalysisTaskSE.h"


class AliAODEvent;
class AliAODVZERO;
class AliVVertex;
class AliFlowEventSimple;
class AliMultSelection;
class AliAnalysisUtils;
class AliPIDResponse;    
class TList;

class AliAnalysisTaskCMEV0 : public AliAnalysisTaskSE {

public:

  AliAnalysisTaskCMEV0();
//AliAnalysisTaskCMEV0(const char *name);
  AliAnalysisTaskCMEV0(const TString name);
  virtual ~AliAnalysisTaskCMEV0();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);


  void    SetInputListNUA(TList *finputNUA)           {this->fListNUACorr       =    finputNUA;}
  void    SetGainCorrZDNP(TList *finputZDN)           {this->fListZDNCorr       =    finputZDN;}
  void    SetFBEfficiencyList(TList *fFBlist)         {this->fListFBHijing      =      fFBlist;}
  void    SetInputListforV0M(TList *finputV0M)        {this->fListV0MCorr       =    finputV0M;}
  void    SetRejectPileUp(Bool_t  pileup)             {this->fRejectPileUp      =       pileup;}
  void    SetRejectPileUpTight(Bool_t  pileupt8)      {this->fRejectPileUpTight =     pileupt8;}
  void    SetStoreTPCQnAvg(Bool_t bstoreTPCQn)        {this->bFillAvgTPCQn      =  bstoreTPCQn;} 
  void    SetFillNUAHist(Bool_t bfillNUAhist)         {this->bFillEtaPhiNUA     =  bfillNUAhist;} 
  void    SetApplyNUACorr(Bool_t fUseNUACorr)         {this->bApplyNUACorr      =  fUseNUACorr;}
  void    SetApplyZDCCorr(Bool_t fUseZDCCorr)         {this->bApplyZDCCorr      =  fUseZDCCorr;}
  void    SetDataSet(TString fdataset)                {this->sDataSet           =     fdataset;}
  void    SetAnalysisSet(TString fanalysisSet)        {this->sAnalysisSet       = fanalysisSet;}
  void    SetCentEstimator(TString centEstim)         {this->sCentEstimator     =    centEstim;}
  void    SetHarmonicN(Int_t harmonic1)               {this->fHarmonicN         =    harmonic1;}
  void    SetHarmonicM(Int_t harmonic2)               {this->fHarmonicM         =    harmonic2;}
  void    SetPsiHarmonic(Int_t nforpsi)               {this->fHarmonicPsi       =      nforpsi;}
  void    SetApplyNUAinEP(Bool_t bApplyNUAEP)         {this->bApplyNUAforEP     =  bApplyNUAEP;}
  void    SetSourceFileNUA(TString sfilenua)          {this->sFileNUA           =     sfilenua;}
  void    SetFillZDNCalHist(Bool_t bfillZDC)          {this->bFillZDCinfo       =     bfillZDC;}
  void    SetSkipNestedLoop(Bool_t bskipNest)         {this->bSkipNestedTrk     =    bskipNest;}
  void    SetMCEffiDimension(TString mcDimen)         {this->sMCdimension       =      mcDimen;}
  void    SetRemoveNegTrkRndm(Bool_t remRndm)         {this->bRemNegTrkRndm     =      remRndm;}
  void    SetApplyV0MCorr(Bool_t  bV0Mcorr)           {this->bApplyV0MCorr      =     bV0Mcorr;}
  void    SetPileUpCutParam(Float_t m,Float_t c)      {this->fPileUpSlopeParm   = m;  this->fPileUpConstParm = c;}
  void    SetTrackFilterBit(Int_t gf)                 {this->gFilterBit         =           gf;}
  void    SetHBTcutParameter(Float_t hb)              {this->fHBTCutValue       =           hb;}
  

private:

  AliAnalysisTaskCMEV0(const AliAnalysisTaskCMEV0& aAnalysisTask);
  AliAnalysisTaskCMEV0& operator=(const AliAnalysisTaskCMEV0& aAnalysisTask);

  Bool_t CheckEventIsPileUp(AliAODEvent* faod);
  Bool_t PileUpMultiVertex(const AliAODEvent* faod);
  double GetWDist(const AliVVertex* v0, const AliVVertex* v1);

  void OpenInfoCalbration(Int_t  run, Float_t fHarmonic);
  void GetV0QvectAndMult(const AliAODVZERO *aodV0,Float_t fHarmonic,Double_t& Qxan,Double_t& Qyan,Double_t& sumMa,Double_t& Qxcn,Double_t& Qycn,Double_t& sumMc);
  void DefineHistograms();
  void GetNUACorrectionHist(Int_t run, TString sfileNUA);
  void GetZDCCorrectionHist(Int_t run);
  void GetV0MCorrectionHist(Int_t run);

  void InitializeRunArray(TString sPeriod);
  Int_t GetCurrentRunIndex(Int_t  run);
  Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign);



  AliFlowEventSimple*         fEvent;         //! input event
  AliMultSelection*   fMultSelection;         //! MultSelection (RUN2 centrality estimator)
  AliAnalysisUtils*    fAnalysisUtil;         //! Event selection

  TList*                 fListHistos;         //! collection of output
  TList*                 fListCalibs;         //! collection of Calib Histos
  TList*                fListNUAHist;         //! collection of NUA Histograms

  TList*               fListFBHijing;         // Hijing Efficiency list
  TList*                fListNUACorr;         // NUA Correction List
  TList*                fListZDNCorr;         // ZDC gain Correction Wgts
  TList*                fListV0MCorr;

  Bool_t               fRejectPileUp;         //
  Bool_t          fRejectPileUpTight;         //
  Bool_t               bFillAvgTPCQn;         //
  Bool_t              bFillEtaPhiNUA;         //
  Bool_t               bApplyNUACorr;         //
  Bool_t               bApplyZDCCorr;         //
  Bool_t              bApplyNUAforEP;         //
  Bool_t                bFillZDCinfo;         //
  Bool_t              bSkipNestedTrk;         //
  Bool_t              bRemNegTrkRndm;         //
  Bool_t               bApplyV0MCorr;         //

  TString                   sDataSet;         // Dataset: 2010, 2011, or 2015.
  TString               sAnalysisSet;         // Values: recenter1,recenter2,analysis1
  TString             sCentEstimator;         // Centrality Estimator
  TString                   sFileNUA;         // NUA source file.
  TString               sMCdimension;         // Dimention for MC effi Correction



  Int_t                  runNums[90];         //!  array of runnumbers
  Int_t                     fRunFlag;         //!  number of total run
  Int_t                   fOldRunNum;         //!
  Int_t                      fievent;         //!  counter of event for cout
  Float_t                    EvtCent;         //!  Event centrality 
  Int_t                   fHarmonicN;         //   Harmonic n
  Int_t                   fHarmonicM;         //   Harmonic m
  Int_t                 fHarmonicPsi;         //   Harmonic psi
  Float_t           fPileUpSlopeParm;         //
  Float_t           fPileUpConstParm;         //
  Int_t                   gFilterBit;         //
  Float_t               fHBTCutValue;         //
  Float_t               fRefMultCorr;         //!
  Float_t                fRefMultRaw;         //!


  TH1F            *fHist_Event_count;         //!  event count with different cuts
  TH1F          *fPileUpMultSelCount;         //!
  TH1F                 *fPileUpCount;         //!
  TH1F              *fTaskConfigParm;         //!  Cut parameters which were used

  //runtime v0 calibration info histograms:
  TH1D*                      fMultV0;         //! profile from V0 multiplicity
  TH1D*                     fQxnmV0A;         //! <Qx2> V0A
  TH1D*                     fQynmV0A;         //! <Qy2> V0A
  TH1D*                     fQxnsV0A;         //! sigma Qx2 V0A
  TH1D*                     fQynsV0A;         //! sigma Qy2 V0A
  TH1D*                     fQxnmV0C;         //! <Qx2> V0C
  TH1D*                     fQynmV0C;         //! <Qy2> V0C
  TH1D*                     fQxnsV0C;         //! sigma Qx2 V0C
  TH1D*                     fQynsV0C;         //! sigma Qy2 V0C
  //---------------------------------

  //--- profiles for TPC Q-vector recenter. (They are filled, not read)
  TProfile2D     *fHCos1nPosChEtaPosVzPos;   //!
  TProfile2D     *fHCos1nPosChEtaNegVzPos;   //!
  TProfile2D     *fHCos1nPosChEtaPosVzNeg;   //!
  TProfile2D     *fHCos1nPosChEtaNegVzNeg;   //!
  TProfile2D     *fHSin1nPosChEtaPosVzPos;   //!
  TProfile2D     *fHSin1nPosChEtaNegVzPos;   //!
  TProfile2D     *fHSin1nPosChEtaPosVzNeg;   //!
  TProfile2D     *fHSin1nPosChEtaNegVzNeg;   //!

  TProfile2D     *fHCos1nNegChEtaPosVzPos;   //!
  TProfile2D     *fHCos1nNegChEtaNegVzPos;   //!
  TProfile2D     *fHCos1nNegChEtaPosVzNeg;   //!
  TProfile2D     *fHCos1nNegChEtaNegVzNeg;   //!
  TProfile2D     *fHSin1nNegChEtaPosVzPos;   //!
  TProfile2D     *fHSin1nNegChEtaNegVzPos;   //!
  TProfile2D     *fHSin1nNegChEtaPosVzNeg;   //!
  TProfile2D     *fHSin1nNegChEtaNegVzNeg;   //!

  TH2F        *fHV0AEventPlaneVsCent;   //!
  TH2F        *fHV0CEventPlaneVsCent;   //!
  TH2F        *fHTPCEventPlaneVsCent;   //!


  TH2F   *fHEnergyZNCvsCent;    //!
  TH2F   *fHEnergyZNAvsCent;    //!
  TH2F   *fHEnergyZPCvsCent;    //!
  TH2F   *fHEnergyZPAvsCent;    //!

  TProfile2D  *fHEnergyZNCvsCentRun;   //!
  TProfile2D  *fHEnergyZNAvsCentRun;   //!
  TProfile2D  *fHEnergyZPCvsCentRun;   //!
  TProfile2D  *fHEnergyZPAvsCentRun;   //!

  TH2F *fHEnergyZPCvsZPA;  //!
  TH2F *fHEnergyZNCvsZNA;  //!

  TH1F *hUnderOverBinNUApos; //!   //temporary Debug, remove for stable code
  TH1F *hUnderOverBinNUAneg; //!   //temporary Debug, remove for stable code

  TH1F   *fHCentBinTrkRecenter; //!

  TH2D          *fHCorrectZDNP; //!

  TProfile2D  *fV0AQ2xVsCentRun; //!
  TProfile2D  *fV0AQ2yVsCentRun; //!
  TProfile2D  *fV0CQ2xVsCentRun; //!
  TProfile2D  *fV0CQ2yVsCentRun; //!

  TProfile2D  *fV0AQ3xVsCentRun; //!
  TProfile2D  *fV0AQ3yVsCentRun; //!
  TProfile2D  *fV0CQ3xVsCentRun; //!
  TProfile2D  *fV0CQ3yVsCentRun; //!

  TProfile2D  *fTPCQ2xVsCentRun; //!
  TProfile2D  *fTPCQ2yVsCentRun; //!
  TProfile2D  *fTPCQ3xVsCentRun; //!
  TProfile2D  *fTPCQ3yVsCentRun; //!

  TList           *mListNUAPos; //!
  TList           *mListNUANeg; //!
  TFile            *fileNUApos; //!
  TFile            *fileNUAneg; //!

  TProfile2D     *fAvgMultCentRun;   //!
  TProfile2D  *fAvgWgtMultCentRun;   //!
  TProfile2D   *fAvgPOIposCentRun;   //!
  TProfile2D   *fAvgPOInegCentRun;   //!
  TProfile2D    *fAvgPOIPPCentRun;   //! same sign Pos-Pos pairs
  TProfile2D    *fAvgPOINNCentRun;   //! same sign Neg-Neg pairs
  TProfile2D    *fAvgPOIOSCentRun;   //! opposite sign pairs

  TProfile2D      *fV0MultChVsRun;   //!

  TH1F           *fEventStatvsRun;   //!
  TH1F        *fEtaBinFinderForQA;   //!
  TH1F        *fVzBinFinderForNUA;   //!

  TH2F            *fHistVtxZvsRun; //!
  TH2F            *fHistVtxXvsRun; //!
  TH2F            *fHistVtxYvsRun; //!

  TProfile2D    *fRejectRatioVsCR; //!

  TH2F            *fCentDistvsRun; //! 
  TH1D              *fHCorrectV0M; //!   for V0-Mult Gain Correction per channel.
  TH2D           *fHAvgerageQnV0A; //!   V0A Average <Qn>, n=2,3
  TH2D           *fHAvgerageQnV0C; //!   V0C Average <Qn>, n=2,3



  TProfile2D     *fCentV0MvsVzRun; //!
  TProfile2D      *fCent3pvsVzRun; //!

  TH2F         *fRefMultCorrvsRaw; //!
  TH2F           *fTPCvsGlobalTrk; //!
  TH2F         *fTPCuncutvsGlobal; //!
  TH1F             *fGlobalTracks; //!
  TH2F             *fTPCvsITSfb96; //!
  TH2F             *fTPCvsITSfb32; //!
  TH2F           *fTPCFEvsITSfb96; //!
  TH3F           *fCentCL1vsVzRun; //!
  TH1F             *fdPhiFemtoCut; //!
  TH1F            *fVzDistribuion; //!
  //TRandom3                fRand; //!













  //  [ Arrays of Histrograms here: ]


  TH3D        *fHCorrectNUApos[5]; //! 5 centrality bin
  TH3D        *fHCorrectNUAneg[5]; //! 5 centrality bin


  TH3F         *fHist3DEtaPhiVz_Pos_Run[5][90];  //! 5 centrality bin 90 Bins for Run. NUA
  TH3F         *fHist3DEtaPhiVz_Neg_Run[5][90];  //! 5 centrality bin 90 Bins for Run. NUA

  TH1D         *fFB_Efficiency_Cent[10];  //!

  TH3F          *fFB_Efficiency_Pos[10];  //!  3d correction Map
  TH3F          *fFB_Efficiency_Neg[10];  //!


  //CME using scalar product method:
  TProfile     *fHist_Corr3p_SP_Norm_PN[2][3];  //! Norm = 10 centrality bins along X
  TProfile     *fHist_Corr3p_SP_Norm_PP[2][3];  //!
  TProfile     *fHist_Corr3p_SP_Norm_NN[2][3];  //!
  TProfile     *fHist_Reso2n_SP_Norm_Det[2][3]; //! 

  //CME Using Event plane method:
  TProfile     *fHist_Corr3p_EP_Norm_PN[2][3];  //! 
  TProfile     *fHist_Corr3p_EP_Norm_PP[2][3];  //!
  TProfile     *fHist_Corr3p_EP_Norm_NN[2][3];  //!
  TProfile     *fHist_Reso2n_EP_Norm_Det[2][3]; //! 

  //CME(EP) vs Refmult:
  TProfile     *fHist_Corr3p_EP_Refm_PN[2][3];  //! 
  TProfile     *fHist_Corr3p_EP_Refm_PP[2][3];  //!
  TProfile     *fHist_Corr3p_EP_Refm_NN[2][3];  //!
  TProfile     *fHist_Reso2n_EP_Refm_Det[2][3]; //! 



  TProfile2D      *fHist_Corr3p_vsRun_EP_PN[2];  // 0=V0A,1=V0C
  TProfile2D      *fHist_Corr3p_vsRun_EP_PP[2];  //
  TProfile2D      *fHist_Corr3p_vsRun_EP_NN[2];  //

  //CME ZDN correlator: Spectator neutron
  TProfile2D     *fHist_Corr3p_ZDN_SP_PN[3];  //! Norm = 10 centrality bins along X
  TProfile2D     *fHist_Corr3p_ZDN_SP_PP[3];  //!
  TProfile2D     *fHist_Corr3p_ZDN_SP_NN[3];  //!
  TProfile2D     *fHist_Reso2n_ZDN_SP_Det[3]; //! 

  //CME pT differential Histograms:

  //(pT_A + pT_B)/2.0
  TProfile     *fHist_Corr3p_pTSum_EP_V0A_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_pTSum_EP_V0A_PP[2][6]; //!
  TProfile     *fHist_Corr3p_pTSum_EP_V0A_NN[2][6]; //!
  TProfile     *fHist_Corr3p_pTSum_EP_V0C_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_pTSum_EP_V0C_PP[2][6]; //!
  TProfile     *fHist_Corr3p_pTSum_EP_V0C_NN[2][6]; //!

  // |(pT_A - pT_B)|
  TProfile     *fHist_Corr3p_pTDiff_EP_V0A_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_pTDiff_EP_V0A_PP[2][6]; //!
  TProfile     *fHist_Corr3p_pTDiff_EP_V0A_NN[2][6]; //!
  TProfile     *fHist_Corr3p_pTDiff_EP_V0C_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_pTDiff_EP_V0C_PP[2][6]; //!
  TProfile     *fHist_Corr3p_pTDiff_EP_V0C_NN[2][6]; //!
 
  // |(Eta_A - Eta_B)|
  TProfile     *fHist_Corr3p_EtaDiff_EP_V0A_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_EtaDiff_EP_V0A_PP[2][6]; //!
  TProfile     *fHist_Corr3p_EtaDiff_EP_V0A_NN[2][6]; //!
  TProfile     *fHist_Corr3p_EtaDiff_EP_V0C_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_EtaDiff_EP_V0C_PP[2][6]; //!
  TProfile     *fHist_Corr3p_EtaDiff_EP_V0C_NN[2][6]; //!
 

  //QA: pos/neg ratio:
  TH3F         *fHistChPosvsEtaPtRun[10];  //! 10 Centrality Bin
  TH3F         *fHistChNegvsEtaPtRun[10];  //!

  //QA eta dependence:
  TProfile  *fHist_Corr3p_QAEta_SP_V0A_PN[2]; //!  Pos and Neg Mag field.
  TProfile  *fHist_Corr3p_QAEta_SP_V0A_PP[2]; //!
  TProfile  *fHist_Corr3p_QAEta_SP_V0A_NN[2]; //!

  TProfile  *fHist_Corr3p_QAEta_SP_V0C_PN[2]; //!  Pos and Neg Mag field.
  TProfile  *fHist_Corr3p_QAEta_SP_V0C_PP[2]; //!
  TProfile  *fHist_Corr3p_QAEta_SP_V0C_NN[2]; //!

  //--- profiles for TPC <Q> (They are filled)
  TProfile2D     *fHCos1nPosChEtaVz[4];   //!
  TProfile2D     *fHCos2nPosChEtaVz[4];   //!
  TProfile2D     *fHCos3nPosChEtaVz[4];   //!
  TProfile2D     *fHCos4nPosChEtaVz[4];   //!
  TProfile2D     *fHSin1nPosChEtaVz[4];   //!
  TProfile2D     *fHSin2nPosChEtaVz[4];   //!
  TProfile2D     *fHSin3nPosChEtaVz[4];   //!
  TProfile2D     *fHSin4nPosChEtaVz[4];   //!

  TProfile2D     *fHCos1nNegChEtaVz[4];   //!
  TProfile2D     *fHCos2nNegChEtaVz[4];   //!
  TProfile2D     *fHCos3nNegChEtaVz[4];   //!
  TProfile2D     *fHCos4nNegChEtaVz[4];   //!
  TProfile2D     *fHSin1nNegChEtaVz[4];   //!
  TProfile2D     *fHSin2nNegChEtaVz[4];   //!
  TProfile2D     *fHSin3nNegChEtaVz[4];   //!
  TProfile2D     *fHSin4nNegChEtaVz[4];   //!

  TProfile2D     *fHCos2nDWPosChEtaVz[4];   //!
  TProfile2D     *fHSin2nDWPosChEtaVz[4];   //!
  TProfile2D     *fHCos2nDWNegChEtaVz[4];   //!
  TProfile2D     *fHSin2nDWNegChEtaVz[4];   //!

  //Store Non-isotropic terms:
  TProfile2D   *fHist_NonIso_SP_PP_Mag0[2];  //! Mag0 = B < 0
  TProfile2D   *fHist_NonIso_SP_NN_Mag0[2];  //!
  TProfile2D   *fHist_NonIso_SP_PP_Mag1[2];  //! Mag1 = B > 0
  TProfile2D   *fHist_NonIso_SP_NN_Mag1[2];  //!


  //2particle correlation vs Cent:
  TProfile2D  *fHist_Corr2p_EP_Norm_PN[2];  //!  Two magnetic fields
  TProfile2D  *fHist_Corr2p_EP_Norm_PP[2];  //! 
  TProfile2D  *fHist_Corr2p_EP_Norm_NN[2];  //! 

  //2particle correlation vs Refm:
  TProfile2D  *fHist_Corr2p_EP_Refm_PN[2];  //!  Two magnetic fields
  TProfile2D  *fHist_Corr2p_EP_Refm_PP[2];  //! 
  TProfile2D  *fHist_Corr2p_EP_Refm_NN[2];  //! 




  ClassDef(AliAnalysisTaskCMEV0, 1); // 
};

#endif

