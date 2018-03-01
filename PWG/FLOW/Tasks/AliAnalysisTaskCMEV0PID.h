/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id: $ */

/////////////////////////////////////////////////
// AliAnalysisTaskCMEV0PID:
// analysis task for ZDC gain Equalization
// and CME analysis for identified particles 
// using VZERO (and ZDC) Detector
// Author: Rihan Haque (mhaque@cern.ch)
// and R Sultana
/////////////////////////////////////////////////

#ifndef AliAnalysisTaskCMEV0PID_H
#define AliAnalysisTaskCMEV0PID_H

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




class AliAnalysisTaskCMEV0PID : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskCMEV0PID();
  AliAnalysisTaskCMEV0PID(const char *name);
  virtual ~AliAnalysisTaskCMEV0PID();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * /*option*/);
    
  //User Defined Functions:
  void SetFilterBit(Int_t b)                 {this->fFilterBit = b;}
  void SetNSigmaCutTPC(Float_t b)            {this->fNSigmaCut = b;}
  void SetPtRangeMin(Float_t b)              {this->fMinPtCut  = b;}
  void SetPtRangeMax(Float_t b)              {this->fMaxPtCut  = b;}
  void SetEtaRangeMin(Float_t b)             {this->fMinEtaCut = b;}
  void SetEtaRangeMax(Float_t b)             {this->fMaxEtaCut = b;}
  void SetCollisionSystem(TString s)         {this->sNucleiTP  = s;}

  void SetCentralityPercentileMin(Float_t b) {this->fCentralityPercentMin = b;}
  void SetCentralityPercentileMax(Float_t b) {this->fCentralityPercentMax = b;}
  void SetFBEfficiencyList(TList *flist)     {this->fListFBHijing  =  flist;}
  void SetFlagForMCcorrection(Bool_t b)      {this->bApplyMCcorr   = b;}
  void SetFBEfficiencyFilePath(TString path) {this->sPathOfMCFile  =   path;}
  void SetPileUpCutParam(Float_t m,Float_t c){this->fPileUpSlopeParm = m;  this->fPileUpConstParm = c;}


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
  TList                 *fListFBHijing;       //!

  //histograms:
  TH1F         *fHistTaskConfigParameters;   //! Task input parameters FB / cut values etc.
  TH1F                  *fHistPileUpCount;   //!
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

  Int_t            fSkipOutlierCut;  //!
  Int_t                 fFilterBit;  //
  Float_t               fNSigmaCut;  //
  Float_t                fMinPtCut;  //
  Float_t                fMaxPtCut;  //
  Float_t               fMinEtaCut;  //
  Float_t               fMaxEtaCut;  //
  Float_t    fCentralityPercentMin;  //
  Float_t    fCentralityPercentMax;  //
  Float_t         fPileUpSlopeParm;  //
  Float_t         fPileUpConstParm;  //
  Bool_t              bApplyMCcorr;  //
  TString            sPathOfMCFile;  //
  TString                sNucleiTP;  //
 



  TH1F    *fHistEventCount;   //!    last in the list

  //-------- Arrays ----------
  TH1F           *fHistPtwithTPCNsigma[3];   //!
  TH1F          *fHistPtwithTOFmasscut[3];   //!
  TH1F           *fHistPtwithTOFSignal[3];   //!
  TH2F        *fHistTOFnSigmavsPtAfter[3];   //!
  TH2F        *fHistTPCnSigmavsPtAfter[3];   //!
  TH3F     *fHistTPCTOFnSigmavsPtAfter[3];   //!
  TH2F       *fHistTPCdEdxvsPtPIDAfter[3];   //!

  TH3F        *fHist3DEtaPhiVz_Pos_Run[3][5];  //! 3 particle 5 centrality bin 
  TH3F        *fHist3DEtaPhiVz_Neg_Run[3][5];  //! 3 particle 5 centrality bin 



  TH1D           *fFB_Efficiency_Cent[10];   //!

  //--------- PileUp Functions -----------
  Bool_t CheckEventIsPileUp(AliAODEvent* faod);
  Bool_t PileUpMultiVertex(const AliAODEvent* faod);
  double GetWDist(const AliVVertex* v0, const AliVVertex* v1);
  //----------- other functions ----------
  void  SetUpCentralityOutlierCut();
  void  SetupEventAndTaskConfigInfo();
  void  SetupMCcorrectionMap(TString sMCfilePath);
  Int_t GetCentralityScaled0to10(Float_t fCent);

  AliAnalysisTaskCMEV0PID(const AliAnalysisTaskCMEV0PID &other);
  AliAnalysisTaskCMEV0PID& operator=(const AliAnalysisTaskCMEV0PID &other);    
  ClassDef(AliAnalysisTaskCMEV0PID,1) 

};

#endif
