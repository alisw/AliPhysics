/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id: $ */

/////////////////////////////////////////////////
// AliAnalysisTaskCMEMC:
// analysis task for CME analysis on MC
// Author: Rihan Haque (mhaque@cern.ch)
// and R Sultana
/////////////////////////////////////////////////

#ifndef ALIANALYSISTASKCMEMC_H
#define ALIANALYSISTASKCMEMC_H

#include "AliAnalysisTaskSE.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TProfile.h"
#include "TProfile2D.h"
//#include "TStopwatch.h"

class    AliVEvent;      
class    AliVVertex;    
class    AliESDEvent;       
class    AliAODEvent;      
class    AliPIDResponse;    
class    AliMultSelection;    
class    AliAnalysisUtils;




class AliAnalysisTaskCMEMC : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskCMEMC();
  AliAnalysisTaskCMEMC(const char *name);
  virtual ~AliAnalysisTaskCMEMC();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * /*option*/);
    
  //User Defined Functions:
  void SetFilterBit(Int_t b)                  {this->fFilterBit = b;}
  void SetNSigmaCutTPC(Double_t b)            {this->fNSigmaCut = b;}
  void SetPtRangeMin(Double_t b)              {this->fMinPtCut  = b;}
  void SetPtRangeMax(Double_t b)              {this->fMaxPtCut  = b;}
  void SetEtaRangeMin(Double_t b)             {this->fMinEtaCut = b;}
  void SetEtaRangeMax(Double_t b)             {this->fMaxEtaCut = b;}

  void SetVzRangeMin(Double_t b)              {this->fMinVzCut  = b;}
  void SetVzRangeMax(Double_t b)              {this->fMaxVzCut  = b;}

  void SetTrackCutdEdxMin(Float_t   dEd)      {this->fdEdxMin = dEd; }
  void SetTrackCutDCAxyMax(Double_t  dc)      {this->fDCAxyMax = dc; }
  void SetTrackCutDCAzMax(Double_t   dc)      {this->fDCAzMax  = dc; }
  void SetTrackCutChi2Min(Double_t  chi)      {this->fTrkChi2Min = chi; }
  void SetTrackCutNclusterMin(Int_t ncl)      {this->fTPCclustMin = ncl; }
  void SetFlagUseKinkTracks(Bool_t kink)      {this->bUseKinkTracks = kink; }


  void SetEventPlaneHarmonic(Int_t pn)            {this->gPsiN     = pn;}
  void SetHarmonicsFor3Particle(Int_t n, Int_t m) {this->gN = n; this->gM = m;}




  void SetCentralityEstimator(TString sCent)  {this->sCentrEstimator = sCent;} 
  void SetCollisionSystem(TString s)          {this->sNucleiTP  = s;}
  void SetCentralityPercentileMin(Double_t b) {this->fCentralityPercentMin = b;}
  void SetCentralityPercentileMax(Double_t b) {this->fCentralityPercentMax = b;}
  void SetPileUpCutParam(Float_t m,Float_t c) {this->fPileUpSlopeParm = m;  this->fPileUpConstParm = c;}

  void SetFlagForMCcorrection(Bool_t b)       {this->bApplyMCcorr    = b;}
  void SetFlagV0MGainCorr(Bool_t b)           {this->bV0MGainCorr    = b;}
  void SetFlagSkipPileUpCuts(Bool_t b)        {this->bSkipPileUpCut  = b;}
  void SetFlagFillNUAforPID(Bool_t b)         {this->bFillNUAHistPID = b;}




  void SetFBEfficiencyFilePath(TString path)  {this->sPathOfMCFile   = path;}
  //void SetFBEfficiencyList(TList *flist)    {this->fListFBHijing   = flist;}
  //void SetListForNUACorr(TList *fln)        {this->fListNUACorr    = fln;}
  //void SetListForV0MCorr(TList *flv)        {this->fListV0MCorr    = flv;}
  void SetFBEfficiencyList(TList *flist)      {this->fListFBHijing   = (TList *) flist->Clone(); }
  void SetListForNUACorr(TList *flist)        {this->fListNUACorr    = (TList *) flist->Clone(); }
  void SetListForV0MCorr(TList *flist)        {this->fListV0MCorr    = (TList *) flist->Clone(); }


 protected:


 private:

  AliVEvent             *fVevent;             //! event
  AliESDEvent           *fESD;                //! esd
  AliAODEvent           *fAOD;                //! aod
  AliPIDResponse        *fPIDResponse;        //! PID response Handler
  AliMultSelection      *fMultSelection;      //!
  AliAnalysisUtils      *fAnalysisUtil;       //! Event selection
  TList                 *fListHist;           //!
  TFile                 *mfileFBHijing;       //!
  TList                 *fListFBHijing;       //
  TList                 *fListNUACorr;        //
  TList                 *fListV0MCorr;        //

  //histograms:
  TH1F         *fHistTaskConfigParameters;   //! Task input parameters FB / cut values etc.
  TH1F                  *fHistPileUpCount;   //!
  TH1F               *fHistMultSelPUCount;   //!
  TH2F                  *fHistEtaPtBefore;   //! Eta-Pt acceptance
  TH2F                   *fHistEtaPtAfter;   //! Eta-Pt acceptance
  TH2F        *fHistTPCvsGlobalMultBefore;   //!
  TH2F         *fHistTPCvsGlobalMultAfter;   //!

  TH2F             *fHistTPCdEdxvsPBefore;   //!
  TH2F              *fHistTPCdEdxvsPAfter;   //!
  TH2F             *fHistTOFBetavsPBefore;   //!
  TH2F              *fHistTOFBetavsPAfter;   //!
  TH2F            *fHistTOFMassvsPtBefore;   //!

  TH2F                *fHistTOFMatchCount;   //!    

  TH2F            *fHistTPCVsESDTrkBefore;   //!  
  TH2F             *fHistTPCVsESDTrkAfter;   //!    

  TH2F           *fHistTPConlyVsCL1Before;   //!  
  TH2F           *fHistTPConlyVsV0MBefore;   //!    
  TH2F            *fHistTPConlyVsCL1After;   //!    
  TH2F            *fHistTPConlyVsV0MAfter;   //!   
  TH2F            *fHistGlobalVsV0MBefore;   //!   
  TH2F             *fHistGlobalVsV0MAfter;   //!   

  TH2F              *fHistRawVsCorrMultFB;   //!  
  TH2F                *hCentvsTPCmultCuts;   //! 

  TH2F             *fHV0AEventPlaneVsCent;   //!
  TH2F             *fHV0CEventPlaneVsCent;   //!
  TH2F            *fHTPCAEventPlaneVsCent;   //!
  TH2F            *fHTPCCEventPlaneVsCent;   //!
  TH2F             *fHTPCEventPlaneVsCent;   //!   Full Event plane
  TH2F                    *fV0MultChVsRun;   //!   To fill VOM multiplicity 
  
  TH1F                   *fCentDistBefore;   //!   without PileUp cut
  TH1F                    *fCentDistAfter;   //!   with PileUp cut
  TH1D                      *fHCorrectV0M;   //!   To read Gain Correction file.
  TH2D                   *fHAvgerageQnV0A;   //!   V0A Average <Qn>, n=2,3
  TH2D                   *fHAvgerageQnV0C;   //!   V0C Average <Qn>, n=2,3

  TH1D                *fHCentWeightForRun;   //!   Event weights for non-flat centrality 

  TH2F                 *fQAEtaPhiAfterNUA;   //!
  TH2F             *fQAEtaPhiAfterNUAPion;   //!
  TH2F             *fQAEtaPhiAfterNUAKaon;   //!
  TH2F           *fQAEtaPhiAfterNUAProton;   //!






  TProfile              *fV0AQ2xVsCentRun; //!
  TProfile              *fV0AQ2yVsCentRun; //!
  TProfile              *fV0CQ2xVsCentRun; //!
  TProfile              *fV0CQ2yVsCentRun; //!
  TProfile              *fV0AQ3xVsCentRun; //!
  TProfile              *fV0AQ3yVsCentRun; //!
  TProfile              *fV0CQ3xVsCentRun; //!
  TProfile              *fV0CQ3yVsCentRun; //!


  TProfile              *fTPCAQ2xVsCentRun; //!
  TProfile              *fTPCAQ2yVsCentRun; //!
  TProfile              *fTPCCQ2xVsCentRun; //!
  TProfile              *fTPCCQ2yVsCentRun; //!
  TProfile              *fTPCAQ3xVsCentRun; //!
  TProfile              *fTPCAQ3yVsCentRun; //!
  TProfile              *fTPCCQ3xVsCentRun; //!
  TProfile              *fTPCCQ3yVsCentRun; //!
  TProfile              *fTPCAQ4xVsCentRun; //!
  TProfile              *fTPCAQ4yVsCentRun; //!
  TProfile              *fTPCCQ4xVsCentRun; //!
  TProfile              *fTPCCQ4yVsCentRun; //!

  TProfile              *fTPCFQ2xVsCentRun; //!
  TProfile              *fTPCFQ2yVsCentRun; //!
//TProfile              *fTPCFQ3xVsCentRun; //!
//TProfile              *fTPCFQ3yVsCentRun; //!


  //TStopwatch                 watch;  //!

  Int_t                 fFilterBit;  //
  Int_t                         gN;  //
  Int_t                         gM;  //
  Int_t                      gPsiN;  //
  Int_t                 fOldRunNum;  //!
  Int_t                fEventCount;  //!
  Int_t               fTPCclustMin;  // 


  Double_t               fNSigmaCut;  //
  Double_t                fMinPtCut;  //
  Double_t                fMaxPtCut;  //
  Double_t               fMinEtaCut;  //
  Double_t               fMaxEtaCut;  //
  Double_t              fTrkChi2Min;  //  
  Double_t                fDCAxyMax;  //  
  Double_t                 fDCAzMax;  //  
  Float_t                  fdEdxMin;  //  

  Double_t    fCentralityPercentMin;  //
  Double_t    fCentralityPercentMax;  //
  Double_t         fPileUpSlopeParm;  //
  Double_t         fPileUpConstParm;  //

  Float_t                 fMinVzCut;  //
  Float_t                 fMaxVzCut;  //


  Bool_t               bApplyMCcorr;  //
  Bool_t               bV0MGainCorr;  //
  Bool_t             bSkipPileUpCut;  //
  Bool_t            bFillNUAHistPID;  //
  Bool_t             bUseKinkTracks;  // 

  TString             sPathOfMCFile;  //
  TString                 sNucleiTP;  //
  TString           sCentrEstimator;  //
 


  TH1F            *fHistEventCount;   //!    last in the list












  //-------- Arrays ----------
  TH1F           *fHistPtwithoutPIDcut[4][10];   //!
  TH1F           *fHistPtwithTPCNsigma[4][10];   //!
  TH1F        *fHistPtwithTPCTOFNsigma[4][10];   //!
  TH1F           *fHistPtwithNoPIDCuts[4][10];   //!


  TH1F           *fHistPtwithTOFSignal[3];   //!
  TH2F        *fHistTOFnSigmavsPtAfter[3];   //!
  TH2F        *fHistTPCnSigmavsPtAfter[3];   //!
  TH3F     *fHistTPCTOFnSigmavsPtAfter[3];   //!
  TH2F       *fHistTPCdEdxvsPtPIDAfter[3];   //!

  TProfile                 *V2IntProQC[3]; //!
  TProfile          *fResolution_TPCEP[3];  //!

  TProfile                    *fDnnPro[8]; //!
  TProfile                    *fCMEPro[8]; //!


  
  
  TH3D                  *fHCorrectNUApos[5];   //! 5 centrality bin, read NUA from file
  TH3D                  *fHCorrectNUAneg[5];   //! 5 centrality bin, read NUA from file

  TH3D              *fHCorrectNUAposPion[5];   //! 5 centrality bin, read NUA from file
  TH3D              *fHCorrectNUAnegPion[5];   //! 5 centrality bin, read NUA from file

  TH3D              *fHCorrectNUAposKaon[5];   //! 5 centrality bin, read NUA from file
  TH3D              *fHCorrectNUAnegKaon[5];   //! 5 centrality bin, read NUA from file

  TH3D            *fHCorrectNUAposProton[5];   //! 5 centrality bin, read NUA from file
  TH3D            *fHCorrectNUAnegProton[5];   //! 5 centrality bin, read NUA from file





  ///////////// CORRELATORS ////////////////
  




  // 3p correlator vs Centrality, EP method:
  //Charge:
  TProfile     *fHist_Corr3p_EP_Norm_PN[2][4];  //! 
  TProfile     *fHist_Corr3p_EP_Norm_PP[2][4];  //!
  TProfile     *fHist_Corr3p_EP_Norm_NN[2][4];  //!
  TProfile     *fHist_Reso2n_EP_Norm_Det[2][4]; //! 
  //Pion:
  TProfile     *fHist_Corr3p_Pion_EP_Norm_PN[2][4];  //! 
  TProfile     *fHist_Corr3p_Pion_EP_Norm_PP[2][4];  //!
  TProfile     *fHist_Corr3p_Pion_EP_Norm_NN[2][4];  //!
  //Kaon:
  TProfile     *fHist_Corr3p_Kaon_EP_Norm_PN[2][4];  //! 
  TProfile     *fHist_Corr3p_Kaon_EP_Norm_PP[2][4];  //!
  TProfile     *fHist_Corr3p_Kaon_EP_Norm_NN[2][4];  //!
  //Proton:
  TProfile     *fHist_Corr3p_Proton_EP_Norm_PN[2][4];  //! 
  TProfile     *fHist_Corr3p_Proton_EP_Norm_PP[2][4];  //!
  TProfile     *fHist_Corr3p_Proton_EP_Norm_NN[2][4];  //!



  // 3p correlator vs RefMult, EP method:
  //Charge:
  TProfile     *fHist_Corr3p_EP_Refm_PN[2][4];  //! 
  TProfile     *fHist_Corr3p_EP_Refm_PP[2][4];  //!
  TProfile     *fHist_Corr3p_EP_Refm_NN[2][4];  //!
  TProfile     *fHist_Reso2n_EP_Refm_Det[2][4]; //! 
  //Pion:
  TProfile     *fHist_Corr3p_Pion_EP_Refm_PN[2][4];  //! 
  TProfile     *fHist_Corr3p_Pion_EP_Refm_PP[2][4];  //!
  TProfile     *fHist_Corr3p_Pion_EP_Refm_NN[2][4];  //!
  //Kaon:
  TProfile     *fHist_Corr3p_Kaon_EP_Refm_PN[2][4];  //! 
  TProfile     *fHist_Corr3p_Kaon_EP_Refm_PP[2][4];  //!
  TProfile     *fHist_Corr3p_Kaon_EP_Refm_NN[2][4];  //!
  //Proton:
  TProfile     *fHist_Corr3p_Proton_EP_Refm_PN[2][4];  //! 
  TProfile     *fHist_Corr3p_Proton_EP_Refm_PP[2][4];  //!
  TProfile     *fHist_Corr3p_Proton_EP_Refm_NN[2][4];  //!










  //   2p correlator vs Centrality  EP method:
  TProfile     *fHist_Corr2p_EP_Norm_PN[2][4];  //! 
  TProfile     *fHist_Corr2p_EP_Norm_PP[2][4];  //!
  TProfile     *fHist_Corr2p_EP_Norm_NN[2][4];  //!
  //Pion
  TProfile     *fHist_Corr2p_Pion_EP_Norm_PN[2][4];  //! 
  TProfile     *fHist_Corr2p_Pion_EP_Norm_PP[2][4];  //!
  TProfile     *fHist_Corr2p_Pion_EP_Norm_NN[2][4];  //!
  //Kaon
  TProfile     *fHist_Corr2p_Kaon_EP_Norm_PN[2][4];  //! 
  TProfile     *fHist_Corr2p_Kaon_EP_Norm_PP[2][4];  //!
  TProfile     *fHist_Corr2p_Kaon_EP_Norm_NN[2][4];  //!
  //Proton
  TProfile     *fHist_Corr2p_Proton_EP_Norm_PN[2][4];  //! 
  TProfile     *fHist_Corr2p_Proton_EP_Norm_PP[2][4];  //!
  TProfile     *fHist_Corr2p_Proton_EP_Norm_NN[2][4];  //!


  // 2p correlator vs RefMult, EP method:
  //Charge:
  TProfile     *fHist_Corr2p_EP_Refm_PN[2][4];  //! 
  TProfile     *fHist_Corr2p_EP_Refm_PP[2][4];  //!
  TProfile     *fHist_Corr2p_EP_Refm_NN[2][4];  //!
  //Pion:
  TProfile     *fHist_Corr2p_Pion_EP_Refm_PN[2][4];  //! 
  TProfile     *fHist_Corr2p_Pion_EP_Refm_PP[2][4];  //!
  TProfile     *fHist_Corr2p_Pion_EP_Refm_NN[2][4];  //!
  //Kaon:
  TProfile     *fHist_Corr2p_Kaon_EP_Refm_PN[2][4];  //! 
  TProfile     *fHist_Corr2p_Kaon_EP_Refm_PP[2][4];  //!
  TProfile     *fHist_Corr2p_Kaon_EP_Refm_NN[2][4];  //!
  //Proton:
  TProfile     *fHist_Corr2p_Proton_EP_Refm_PN[2][4];  //! 
  TProfile     *fHist_Corr2p_Proton_EP_Refm_PP[2][4];  //!
  TProfile     *fHist_Corr2p_Proton_EP_Refm_NN[2][4];  //!








  //3particle Differential Histograms:
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
 

  //CME PID  differential Histograms:
  //(pT_A + pT_B)/2.0
  TProfile     *fHist_Corr3p_Pion_pTSum_EP_V0A_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Pion_pTSum_EP_V0A_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Pion_pTSum_EP_V0A_NN[2][6]; //!
  TProfile     *fHist_Corr3p_Pion_pTSum_EP_V0C_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Pion_pTSum_EP_V0C_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Pion_pTSum_EP_V0C_NN[2][6]; //!
  // |(pT_A - pT_B)|
  TProfile     *fHist_Corr3p_Pion_pTDiff_EP_V0A_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Pion_pTDiff_EP_V0A_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Pion_pTDiff_EP_V0A_NN[2][6]; //!
  TProfile     *fHist_Corr3p_Pion_pTDiff_EP_V0C_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Pion_pTDiff_EP_V0C_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Pion_pTDiff_EP_V0C_NN[2][6]; //!
  // |(Eta_A - Eta_B)|
  TProfile     *fHist_Corr3p_Pion_EtaDiff_EP_V0A_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Pion_EtaDiff_EP_V0A_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Pion_EtaDiff_EP_V0A_NN[2][6]; //!
  TProfile     *fHist_Corr3p_Pion_EtaDiff_EP_V0C_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Pion_EtaDiff_EP_V0C_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Pion_EtaDiff_EP_V0C_NN[2][6]; //!
 

  //(pT_A + pT_B)/2.0
  TProfile     *fHist_Corr3p_Kaon_pTSum_EP_V0A_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Kaon_pTSum_EP_V0A_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Kaon_pTSum_EP_V0A_NN[2][6]; //!
  TProfile     *fHist_Corr3p_Kaon_pTSum_EP_V0C_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Kaon_pTSum_EP_V0C_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Kaon_pTSum_EP_V0C_NN[2][6]; //!
  // |(pT_A - pT_B)|
  TProfile     *fHist_Corr3p_Kaon_pTDiff_EP_V0A_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Kaon_pTDiff_EP_V0A_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Kaon_pTDiff_EP_V0A_NN[2][6]; //!
  TProfile     *fHist_Corr3p_Kaon_pTDiff_EP_V0C_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Kaon_pTDiff_EP_V0C_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Kaon_pTDiff_EP_V0C_NN[2][6]; //!
  // |(Eta_A - Eta_B)|
  TProfile     *fHist_Corr3p_Kaon_EtaDiff_EP_V0A_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Kaon_EtaDiff_EP_V0A_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Kaon_EtaDiff_EP_V0A_NN[2][6]; //!
  TProfile     *fHist_Corr3p_Kaon_EtaDiff_EP_V0C_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Kaon_EtaDiff_EP_V0C_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Kaon_EtaDiff_EP_V0C_NN[2][6]; //!
 

  //(pT_A + pT_B)/2.0
  TProfile     *fHist_Corr3p_Proton_pTSum_EP_V0A_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Proton_pTSum_EP_V0A_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Proton_pTSum_EP_V0A_NN[2][6]; //!
  TProfile     *fHist_Corr3p_Proton_pTSum_EP_V0C_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Proton_pTSum_EP_V0C_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Proton_pTSum_EP_V0C_NN[2][6]; //!
  // |(pT_A - pT_B)|
  TProfile     *fHist_Corr3p_Proton_pTDiff_EP_V0A_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Proton_pTDiff_EP_V0A_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Proton_pTDiff_EP_V0A_NN[2][6]; //!
  TProfile     *fHist_Corr3p_Proton_pTDiff_EP_V0C_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Proton_pTDiff_EP_V0C_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Proton_pTDiff_EP_V0C_NN[2][6]; //!
  // |(Eta_A - Eta_B)|
  TProfile     *fHist_Corr3p_Proton_EtaDiff_EP_V0A_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Proton_EtaDiff_EP_V0A_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Proton_EtaDiff_EP_V0A_NN[2][6]; //!
  TProfile     *fHist_Corr3p_Proton_EtaDiff_EP_V0C_PN[2][6]; //! 
  TProfile     *fHist_Corr3p_Proton_EtaDiff_EP_V0C_PP[2][6]; //!
  TProfile     *fHist_Corr3p_Proton_EtaDiff_EP_V0C_NN[2][6]; //!
 



  //2 particle Differential Histograms:


  //(pT_A + pT_B)/2.0
  TProfile     *fHist_Corr2p_pTSum_EP_Harm1_PN[2][6]; //! 
  TProfile     *fHist_Corr2p_pTSum_EP_Harm1_PP[2][6]; //!
  TProfile     *fHist_Corr2p_pTSum_EP_Harm1_NN[2][6]; //!
  TProfile     *fHist_Corr2p_pTSum_EP_Harm2_PN[2][6]; //! 
  TProfile     *fHist_Corr2p_pTSum_EP_Harm2_PP[2][6]; //!
  TProfile     *fHist_Corr2p_pTSum_EP_Harm2_NN[2][6]; //!
  // |(pT_A - pT_B)|
  TProfile     *fHist_Corr2p_pTDiff_EP_Harm1_PN[2][6]; //! 
  TProfile     *fHist_Corr2p_pTDiff_EP_Harm1_PP[2][6]; //!
  TProfile     *fHist_Corr2p_pTDiff_EP_Harm1_NN[2][6]; //!
  TProfile     *fHist_Corr2p_pTDiff_EP_Harm2_PN[2][6]; //! 
  TProfile     *fHist_Corr2p_pTDiff_EP_Harm2_PP[2][6]; //!
  TProfile     *fHist_Corr2p_pTDiff_EP_Harm2_NN[2][6]; //!
  // |(Eta_A - Eta_B)|
  TProfile     *fHist_Corr2p_EtaDiff_EP_Harm1_PN[2][6]; //! 
  TProfile     *fHist_Corr2p_EtaDiff_EP_Harm1_PP[2][6]; //!
  TProfile     *fHist_Corr2p_EtaDiff_EP_Harm1_NN[2][6]; //!
  TProfile     *fHist_Corr2p_EtaDiff_EP_Harm2_PN[2][6]; //! 
  TProfile     *fHist_Corr2p_EtaDiff_EP_Harm2_PP[2][6]; //!
  TProfile     *fHist_Corr2p_EtaDiff_EP_Harm2_NN[2][6]; //!
 


   




  // NUA histograms:
  TH3F        *fHist3DEtaPhiVz_Pos_Run[4][5];  //! 4 particle 5 centrality bin 
  TH3F        *fHist3DEtaPhiVz_Neg_Run[4][5];  //! 4 particle 5 centrality bin 



  TH1D           *fFB_Efficiency_Cent[10];   //!   for charge
  TH1D      *fFB_Efficiency_Pion_Cent[10];   //!   
  TH1D      *fFB_Efficiency_Kaon_Cent[10];   //!   
  TH1D    *fFB_Efficiency_Proton_Cent[10];   //!   


  //--------- PileUp Functions -----------
  Bool_t CheckEventIsPileUp(AliAODEvent* faod);
  Bool_t PileUpMultiVertex(const AliAODEvent* faod);
  Bool_t IsLabelUsed(TArrayI labelArray, Int_t label);
  double GetWDist(const AliVVertex* v0, const AliVVertex* v1);

  //for NUA and gain corrections:
  void  GetNUACorrectionHist(Int_t run);
  void  GetV0MCorrectionHist(Int_t run);


  
  //----------- other functions ----------
  void  SetUpCentralityOutlierCut();
  void  SetupEventAndTaskConfigInfo();
  void  SetupMCcorrectionMap();
  Int_t GetCentralityScaled0to10(Double_t fCent);




  
  AliAnalysisTaskCMEMC(const AliAnalysisTaskCMEMC &other);
  AliAnalysisTaskCMEMC &operator=(const AliAnalysisTaskCMEMC &other);    
  ClassDef(AliAnalysisTaskCMEMC, 2) 

};

#endif
