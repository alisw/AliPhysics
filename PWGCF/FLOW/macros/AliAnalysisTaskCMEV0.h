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
#include "TProfile2D.h"
#include "TProfile.h"
#include "TString.h"
#include "AliAnalysisTaskSE.h"


class AliAODEvent;
class AliAODVZERO;
class AliVVertex;
class AliFlowEventSimple;
class AliMultSelection;
class AliAnalysisUtils;
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
  void    SetRejectPileUp(Bool_t  pileup)             {this->fRejectPileUp      =       pileup;}
  void    SetRejectPileUpTight(Bool_t  pileupt8)      {this->fRejectPileUpTight =     pileupt8;}
  void    SetStoreTPCQnAvg(Bool_t bstoreTPCQn)        {this->bFillAvgTPCQn      =  bstoreTPCQn;} 
  void    SetFillNUAHist(Bool_t bfillNUAhist)         {this->bFillEtaPhiNUA     =  bfillNUAhist;} 
  void    SetApplyNUACorr(Bool_t fUseNUACorr)         {this->bApplyNUACorr      =  fUseNUACorr;}
  void    SetApplyZDCCorr(Bool_t fUseZDCCorr)         {this->bApplyZDCCorr      =  fUseZDCCorr;}
  void    SetDataSet(TString fdataset)                {this->sDataSet           =     fdataset;}
  void    SetAnalysisSet(TString fanalysisSet)        {this->sAnalysisSet       = fanalysisSet;}
  void    SetCentEstimator(TString centEstim)         {this->sCentEstimator     =    centEstim;}
  void    SetHarmonic(Float_t harmonic)               {this->fHarmonic          =     harmonic;}
  void    SetApplyNUAinEP(Bool_t bApplyNUAEP)         {this->bApplyNUAforEP     =  bApplyNUAEP;}






private:

  AliAnalysisTaskCMEV0(const AliAnalysisTaskCMEV0& aAnalysisTask);
  AliAnalysisTaskCMEV0& operator=(const AliAnalysisTaskCMEV0& aAnalysisTask);

  Bool_t CheckEventIsPileUp(AliAODEvent* faod);
  Bool_t PileUpMultiVertex(const AliAODEvent* faod);
  double GetWDist(const AliVVertex* v0, const AliVVertex* v1);

  void OpenInfoCalbration(Int_t  run);
  void GetV0QvectAndMult(const AliAODVZERO *aodV0,Double_t& Qxan,Double_t& Qyan,Double_t& sumMa,Double_t& Qxcn,Double_t& Qycn,Double_t& sumMc);
  void DefineHistograms();
  void GetNUACorrectionHist(Int_t run, Float_t cent);
  void GetZDCCorrectionHist(Int_t run);

  void InitializeRunArray(TString sPeriod);
  Int_t GetCurrentRunIndex(Int_t  run);


  AliFlowEventSimple*         fEvent;         //! input event
  AliMultSelection*   fMultSelection;         //! MultSelection (RUN2 centrality estimator)
  AliAnalysisUtils*    fAnalysisUtil;         //! Event selection

  TList*                 fListHistos;         //! collection of output
  TList*                 fListCalibs;         //! collection of Calib Histos
  TList*               fListFBHijing;         // Hijing Efficiency list
  TList*                fListNUACorr;         // NUA Correction List
  TList*                fListZDNCorr;         // ZDC gain Correction Wgts

  Bool_t               fRejectPileUp;         //
  Bool_t          fRejectPileUpTight;         //
  Bool_t               bFillAvgTPCQn;         //
  Bool_t              bFillEtaPhiNUA;         //
  Bool_t               bApplyNUACorr;         //
  Bool_t               bApplyZDCCorr;         //
  Bool_t              bApplyNUAforEP;         //


  TString                   sDataSet;         // Dataset: 2010, 2011, or 2015.
  TString               sAnalysisSet;         // Values: recenter1,recenter2,analysis1
  TString             sCentEstimator;         // Centrality Estimator


  Int_t                  runNums[90];         //!  array of runnumbers
  Int_t                     fRunFlag;         //!  number of total run
  Int_t                   fOldRunNum;         //!
  Int_t                      fievent;         //!  counter of event for cout
  Float_t                    EvtCent;         //!  Event centrality 
  Float_t                  fHarmonic;         //!  Harmonic

  TH1F            *fHist_Event_count;         //!  event count with different cuts
  TH1F          *fPileUpMultSelCount;         //!
  TH1F                 *fPileUpCount;         //!
  //runtime v0 calibration info histograms:
  TH1D*                      fMultV0;            //! profile from V0 multiplicity
  TH1D*                     fQxnmV0A;            //! <Qx2> V0A
  TH1D*                     fQynmV0A;            //! <Qy2> V0A
  TH1D*                     fQxnsV0A;            //! sigma Qx2 V0A
  TH1D*                     fQynsV0A;            //! sigma Qy2 V0A
  TH1D*                     fQxnmV0C;            //! <Qx2> V0C
  TH1D*                     fQynmV0C;            //! <Qy2> V0C
  TH1D*                     fQxnsV0C;            //! sigma Qx2 V0C
  TH1D*                     fQynsV0C;            //! sigma Qy2 V0C
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

  TH3D        *fHCorrectNUApos; //!
  TH3D        *fHCorrectNUAneg; //!

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

  TProfile2D  *fV0AQnxVsCentRun; //!
  TProfile2D  *fV0AQnyVsCentRun; //!
  TProfile2D  *fV0CQnxVsCentRun; //!
  TProfile2D  *fV0CQnyVsCentRun; //!
  TProfile2D  *fTPCQnxVsCentRun; //!
  TProfile2D  *fTPCQnyVsCentRun; //!

  TList           *mListNUAPos; //!
  TList           *mListNUANeg; //!
  TFile            *fileNUApos; //!
  TFile            *fileNUAneg; //!







  //  [ Arrays of Histrograms here: ]
  TH3F         *fHist3DEtaPhiVz_Pos_Run[4][90];  //! 4 centrality bin 90 Bins for Run. NUA
  TH3F         *fHist3DEtaPhiVz_Neg_Run[4][90];  //! 4 centrality bin 90 Bins for Run. NUA

  TH1D         *fFB_Efficiency_Cent[10];  //!

  //CME using scalar product method:
  TProfile     *fHist_Corr3p_SP_Norm_PN[3];  //! Norm = 10 centrality bins along X
  TProfile     *fHist_Corr3p_SP_Norm_PP[3];  //!
  TProfile     *fHist_Corr3p_SP_Norm_NN[3];  //!
  TProfile     *fHist_Reso2n_SP_Norm_Det[3]; //! 

  //CME Using Event plane method:
  TProfile     *fHist_Corr3p_EP_Norm_PN[3];  //! 
  TProfile     *fHist_Corr3p_EP_Norm_PP[3];  //!
  TProfile     *fHist_Corr3p_EP_Norm_NN[3];  //!
  TProfile     *fHist_Reso2n_EP_Norm_Det[3]; //! 

  //CME ZDN correlator: Spectator neutron
  TProfile2D     *fHist_Corr3p_ZDN_SP_PN[3];  //! Norm = 10 centrality bins along X
  TProfile2D     *fHist_Corr3p_ZDN_SP_PP[3];  //!
  TProfile2D     *fHist_Corr3p_ZDN_SP_NN[3];  //!
  TProfile2D     *fHist_Reso2n_ZDN_SP_Det[3]; //! 


  //CME pT differential Histograms:

  //(pT_A + pT_B)/2.0
  TProfile     *fHist_Corr3p_pTSum_EP_V0A_PN[6]; //! 
  TProfile     *fHist_Corr3p_pTSum_EP_V0A_PP[6]; //!
  TProfile     *fHist_Corr3p_pTSum_EP_V0A_NN[6]; //!
  TProfile     *fHist_Corr3p_pTSum_EP_V0C_PN[6]; //! 
  TProfile     *fHist_Corr3p_pTSum_EP_V0C_PP[6]; //!
  TProfile     *fHist_Corr3p_pTSum_EP_V0C_NN[6]; //!

  // |(pT_A - pT_B)|
  TProfile     *fHist_Corr3p_pTDiff_EP_V0A_PN[6]; //! 
  TProfile     *fHist_Corr3p_pTDiff_EP_V0A_PP[6]; //!
  TProfile     *fHist_Corr3p_pTDiff_EP_V0A_NN[6]; //!
  TProfile     *fHist_Corr3p_pTDiff_EP_V0C_PN[6]; //! 
  TProfile     *fHist_Corr3p_pTDiff_EP_V0C_PP[6]; //!
  TProfile     *fHist_Corr3p_pTDiff_EP_V0C_NN[6]; //!
 
  // |(Eta_A - Eta_B)|
  TProfile     *fHist_Corr3p_EtaDiff_EP_V0A_PN[6]; //! 
  TProfile     *fHist_Corr3p_EtaDiff_EP_V0A_PP[6]; //!
  TProfile     *fHist_Corr3p_EtaDiff_EP_V0A_NN[6]; //!
  TProfile     *fHist_Corr3p_EtaDiff_EP_V0C_PN[6]; //! 
  TProfile     *fHist_Corr3p_EtaDiff_EP_V0C_PP[6]; //!
  TProfile     *fHist_Corr3p_EtaDiff_EP_V0C_NN[6]; //!
 

  //QA: pos/neg ratio:
  TH3F         *fHistChPosvsEtaPtRun[10];  //! 10 Centrality Bin
  TH3F         *fHistChNegvsEtaPtRun[10];  //!




  ClassDef(AliAnalysisTaskCMEV0, 1); // 
};

#endif

