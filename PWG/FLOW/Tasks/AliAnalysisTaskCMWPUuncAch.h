/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id: $ */

////////////////////////////////////////////////
// AliAnalysisTaskCVE:
// Simple CVE AnalysisTask
// PA: Rihan Haque (mhaque@cern.ch, rihanphys@gmail.com)
// Date: Jan 13, 2020.
// Last Mod: Jan 13, 2020.
/////////////////////////////////////////////////



#ifndef ALIANALYSISTASKCMWPUuncAch_H
#define ALIANALYSISTASKCMWPUuncAch_H

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

class    AliVEvent;      
class    AliVVertex;    
class    AliESDEvent;       
class    AliAODEvent;      
class    AliPIDResponse;    
class    AliMultSelection;    
class    AliAnalysisUtils;





class AliAnalysisTaskCMWPUuncAch : public AliAnalysisTaskSE {

 public:

  //-----> Mandatory Functions:
  AliAnalysisTaskCMWPUuncAch();
  AliAnalysisTaskCMWPUuncAch(const char *name);
  virtual ~AliAnalysisTaskCMWPUuncAch();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * /*option*/);
  
  //-----> User Defined Functions:
  Int_t GetCentralityScaled0to10(Double_t fCent);
  void  SetupEventAndTaskConfigInfo();
  void  SetListForTrkCorr(TList *flist)      {this->fListTRKCorr = (TList *) flist->Clone(); }
  void  SetListForNUACorr(TList *flist)      {this->fListNUACorr = (TList *) flist->Clone(); }
  void  SetListForV0MCorr(TList *flist)      {this->fListV0MCorr = (TList *) flist->Clone(); }
  void  SetParticle(Int_t part)              {this->fParticle    =  part;}
  void  SetCumulantHarmonic(Int_t harm)      {this->gHarmonic    =  harm;}

  
  /******* Event Cut Ranges ******/
  void SetCentralityPercentileMin(Double_t centMin) {this->fCentralityMin  = centMin;}
  void SetCentralityPercentileMax(Double_t centMax) {this->fCentralityMax  = centMax;}
  void SetCentralityEstimator(TString sEstim)       {this->sCentrEstimator = sEstim;}
  void SetVzRangeMin(Double_t vzMin)                {this->fMinVzCut       =  vzMin;}
  void SetVzRangeMax(Double_t vzMax)                {this->fMaxVzCut       =  vzMax;}

  void SetPileUpCutParam(Float_t m,Float_t c) {this->fPileUpSlopeParm = m;  this->fPileUpConstParm = c;}
  void SetFlagSkipPileUpCuts(Bool_t b, Bool_t a)        {this->bSkipPileUpCut  = b; this->usecorrectedach = a;}

  /******* Track Cut Ranges ******/
  void SetNSigmaCutTPC(Double_t     nSigTPC)     {this->fNSigmaTPCCut  =  nSigTPC;}
  void SetNSigmaCutTOF(Double_t     nSigTOF)     {this->fNSigmaTOFCut  =  nSigTOF;}
  void SetTrackCutdEdxMin(Float_t  ndEdxMin)     {this->fdEdxMin    = ndEdxMin;}
  void SetTrackCutChi2Min(Double_t  chi2Min)     {this->fTrkChi2Min =  chi2Min;}
  void SetFlagUseKinkTracks(Bool_t    bKink)     {this->bUseKinkTracks = bKink;}
  
  void SetTrackCutNclusterMin(Int_t ncl)         {this->fTPCclustMin = ncl;}
  void SetFilterBit(Int_t fb)                    {this->fFilterBit   =  fb;}
  void SetEtaRangeMin(Double_t emn)              {this->fMinEtaCut   = emn;}
  void SetEtaRangeMax(Double_t emx)              {this->fMaxEtaCut   = emx;}
  void SetPtRangeMin(Double_t ptL)               {this->fMinPtCut    = ptL;}
  void SetPtRangeMax(Double_t ptH)               {this->fMaxPtCut    = ptH;}
  void SetDCAXYRangeMax(Double_t dcaxy)          {this->fDCAxyMax    = dcaxy;}
  void SetDCAZRangeMax(Double_t dcaz)            {this->fDCAzMax    =  dcaz;}
  void SetChi2Range(Double_t chi2)               {this->fChi2    =  chi2;}
  void SetEtaNeg(Double_t etaL)                  {this->fEtaGapNeg   = etaL;}
  void SetEtaPos(Double_t etaH)                  {this->fEtaGapPos   = etaH;}

  //------ End of user defined function -------







  
 protected:

 private:


  AliVEvent             *fVevent;             //! event
  AliESDEvent           *fESD;                //! esd Event
  AliAODEvent           *fAOD;                //! aod Event
  AliPIDResponse        *fPIDResponse;        //! PID Handler
  AliMultSelection      *fMultSelection;      //! For Centrality 
  AliAnalysisUtils      *fAnalysisUtil;       //! Event Selection Options
  TList                 *fListHist;           //! OutputList

  TList                 *fListTRKCorr;        //  Supplied from Task
  TList                 *fListNUACorr;        //  Supplied from Task
  TList                 *fListV0MCorr;        //  Supplied from Task

  //Set from the AddTask:
  Float_t         fCentralityMin;  //
  Float_t         fCentralityMax;  //
  
  //Track Variables to be used:
  Int_t                 gHarmonic;  //
  Int_t                 fParticle;  //
  Int_t                fFilterBit;  //
  Int_t              fTPCclustMin;  //
  Bool_t           bUseKinkTracks;  //
  
  Float_t           fNSigmaTPCCut;  //
  Float_t           fNSigmaTOFCut;  //
  Float_t               fMinPtCut;  //
  Float_t               fMaxPtCut;  //
  Float_t               fDCAxyMax;  //                                                                                                        
  Float_t               fDCAzMax;  // 
  Float_t                 fChi2; 
  Double_t         fPileUpSlopeParm;  //
  Double_t         fPileUpConstParm;  //
  Bool_t             bSkipPileUpCut;  //
  Bool_t             usecorrectedach;  //
  Double_t              fEtaGapNeg;  //
  Double_t              fEtaGapPos;  //
  Float_t              fMinEtaCut;  //
  Float_t              fMaxEtaCut;  //
  Float_t             fTrkChi2Min;  //
  Float_t                fdEdxMin;  // 
  //Event Variables to be used:
  Float_t               fMinVzCut;  //
  Float_t               fMaxVzCut;  //
  
  TString         sCentrEstimator;  //


  

  /////////// Default HISTOGRAMS /////////
  
  //QA histograms:
  TH1F                  *fCentDistBeforCut;    //! Cent before Any Cut
  TH1F                  *fCentDistAfterCut;    //! Cent After All Cuts
  TH2F                  *fHistEtaPtBeforCut;   //! Eta-Pt before Any Cut
  TH2F                  *fHistEtaPhiBeforCut;  //! Eta-Phi before Any Cut
  TH2F                  *fHistEtaPhiAfterCut;  //! Eta-Phi after trk Cut    
 


  
  //User Defined Histograms: 


  ///Used For Tracking Efficiency:
  TH1D          *fHCorrectMCposChrg;
  TH1D          *fHCorrectMCposPion;
  TH1D          *fHCorrectMCposKaon;
  TH1D          *fHCorrectMCposProt;
  TH1D          *fHCorrectMCnegChrg;
  TH1D          *fHCorrectMCnegPion;
  TH1D          *fHCorrectMCnegKaon;
  TH1D          *fHCorrectMCnegProt;





  
  //QA and Stepcount 
  TH2F            *fHistAChrgVsCent;  //!
  TH2F            *fHistAChrgVsCentunc;  //!

  TH2F            *fHistTPConlyVsCL1Before;   //!  
  TH2F            *fHistTPConlyVsV0MBefore;   //!    
  TH2F            *fHistTPConlyVsCL1After;   //!    
  TH2F            *fHistTPConlyVsV0MAfter;   //!   
  TH2F            *fHistGlobalVsV0MBefore;   //!   
  TH2F            *fHistGlobalVsV0MAfter;   //
  TH2F            *fTPCvsGlobalTrkBefore; //!  Global vs TPC tracks for QA
  TH2F            *fTPCvsGlobalTrkAfter; //!  Global vs TPC tracks for QA
  TH2F            *fHistTPCVsESDTrkBefore;   //!  
  TH2F            *fHistTPCVsESDTrkAfter;   //
  TH1F            *fHistPileUpCount;   //!

  TH1F            *fHistEventCount;   //!
  TH1F  	  *fHCorrectEVNTWGTChrg;   //!   //eventwgt for charge



  ///v2 vs Ach (Results)
  TProfile     *fHistv2AchChrgPos[2][10];
  TProfile     *fHistv2AchPionPos[2][10]; //! [1st] = method, [2nd] = centrality.
  TProfile     *fHistv2AchKaonPos[2][10]; //!
  TProfile     *fHistv2AchProtPos[2][10]; //!
  TProfile     *fHistv2AchChrgNeg[2][10];
  TProfile     *fHistv2AchPionNeg[2][10]; //! [1st] = method, [2nd] = centrality.
  TProfile     *fHistv2AchKaonNeg[2][10]; //!
  TProfile     *fHistv2AchProtNeg[2][10]; //!

  TProfile     *fHistv2AchChrgPosChrgNeg[2][10];
  TProfile     *fHistv2AchPionPosPionNeg[2][10]; //! [1st] = method, [2nd] = centrality.
  TProfile     *fHistv2AchKaonPosKaonNeg[2][10]; //! [1st] = method, [2nd] = centrality.
  TProfile     *fHistv2AchProtPosProtNeg[2][10]; //! [1st] = method, [2nd] = centrality.
    
  TProfile     *fHistv2AchChrgNegChrgPos[2][10];
  TProfile     *fHistv2AchPionNegPionPos[2][10]; //! [1st] = method, [2nd] = centrality.
  TProfile     *fHistv2AchKaonNegKaonPos[2][10];
  TProfile     *fHistv2AchProtNegProtPos[2][10]; //! [1st] = method, [2nd] = centrality.
    


  ///Used For NUA Corrections:
  TH3F          *fHCorrectNUAposChrg[5];   //! [5] = centrality bins
  TH3F          *fHCorrectNUAnegChrg[5];   //! 
  TH3F          *fHCorrectNUAposPion[5];   //! 
  TH3F          *fHCorrectNUAnegPion[5];   //! 
  TH3F          *fHCorrectNUAposKaon[5];   //! 
  TH3F          *fHCorrectNUAnegKaon[5];   //! 
  TH3F          *fHCorrectNUAposProt[5];   //! 
  TH3F          *fHCorrectNUAnegProt[5];   //! 


  /// TO fill NUA for new Cut:
  TH3F          *fHFillNUAPosPID[5];    //! 
  TH3F          *fHFillNUANegPID[5];    //! 

 
  
  TProfile      *fHistEPResolution[2];       //! EP resolution vs Cent
  TProfile      *fHistEPResolutionAch[10];   //! EP resolution vs Ach 
  TProfile      *fHistv2cumAchChrgAll[10];  //! Charge inclusive
  TProfile      *fHistAch[10];

  
  

  ///Custom Functions:
  void  GetNUACorrectionHist(Int_t run=0,Int_t kParticleID=0);
  void  GetEVNTWGTCorrectionHist(Int_t run=0,Int_t kParticleID=0);
  //void  GetV0MCorrectionHist(Int_t run=0);
  void  GetV0MCorrectionHist(Int_t run=0,Int_t kParticleID=0);
  //void  GetMCCorrectionHist(Int_t run=0);
  void  GetMCCorrectionHist(Int_t run=0,Float_t centr=0);
  
  Bool_t CheckEventIsPileUp(AliAODEvent* faod);
  Bool_t PileUpMultiVertex(const AliAODEvent* faod);
  double GetWDist(const AliVVertex* v0, const AliVVertex* v1);


  AliAnalysisTaskCMWPUuncAch(const AliAnalysisTaskCMWPUuncAch &other);
  AliAnalysisTaskCMWPUuncAch &operator=(const AliAnalysisTaskCMWPUuncAch &other);    
  ClassDef(AliAnalysisTaskCMWPUuncAch, 1) 

};

#endif
