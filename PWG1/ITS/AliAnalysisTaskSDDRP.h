#ifndef ALIANALYSISTASKSDDRP
#define ALIANALYSISTASKSDDRP

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysiTaskSDDRP
// AliAnalysisTaskSE to extract from ESD + ESDfreinds + ITS rec points
// performance plots for SDD detector
//
// Author: F. Prino, prino@to.infn.it
//*************************************************************************

class TList;
class TH1F;
class TH2F;
class TTree;
class TString;
class AliESDEvent;
class AliESDfriend;
class AliITSresponseSDD;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskSDDRP : public AliAnalysisTaskSE {

 public:
  
  AliAnalysisTaskSDDRP();
  virtual ~AliAnalysisTaskSDDRP();
  virtual void   UserExec(Option_t *option);
  virtual void   UserCreateOutputObjects();
  virtual void   Terminate(Option_t *option);

  void SetUseITSstandaloneTracks(Bool_t use){
    fUseITSsaTracks=use;
  }
  void SetMinITSPoints(Int_t minp=3){
    fMinITSpts=minp;
  }
  void SetMinTPCPoints(Int_t minp=70){
    fMinTPCpts=minp;
  }
  void SetUseOnlyCINT1BTriggers(Bool_t use=kTRUE){
    fOnlyCINT1BTrig=use;
  }
  void SetMinPfordEdx(Float_t minp=0.5){
    fMinPfordEdx=minp;
  }
  void SetExcludeBadModules(Bool_t opt=kTRUE){
    fExcludeBadMod=opt;
  }
  Bool_t CheckModule(Int_t lay, Int_t lad, Int_t det) const;
  

 private:
  AliAnalysisTaskSDDRP(const AliAnalysisTaskSDDRP &source);
  AliAnalysisTaskSDDRP& operator=(const AliAnalysisTaskSDDRP &source);
  
  TList*  fOutput;          //! QA histos
  TH1F*   fHistNEvents;     //! histo with N of events  
  TH1F*   fHistAllPMod;     //! histo of tracks crossing SDD modules
  TH1F*   fHistGoodPMod;    //! histo of tracks with good point in SDD module
  TH1F*   fHistBadRegMod;   //! histo of tracks crossing bad region of SDD mod.
  TH1F*   fHistMissPMod;    //! histo of tracks with missing point in SDD mod.
  TH1F*   fHistSkippedMod;  //! histo of tracks skipping an SDD module
  TH1F*   fHistOutAccMod;   //! histo of tracks out of accept. in SDD module
  TH1F*   fHistNoRefitMod;  //! histo of points rejected in refit vs. SDD mod.
  TH2F*   fHistdEdxL3VsP;   //! 2D histo of dE/dx vs. momentum -- layer 3
  TH2F*   fHistdEdxL4VsP;   //! 2D histo of dE/dx vs. momentum -- layer 4
  TH2F*   fHistdEdxVsMod;   //! 2D histo of dE/dx vs. module number
  TH1F*   fRecPMod;         //! histo with module occupancy (RecP) 
  TH1F*   fTrackPMod;       //! histo with module occupancy (TrP)
  TH1F*   fGoodAnMod;       //! histo good anodes per module 
  TH1F*   fRecPLadLay3;     //! histo with ladder occupancy on layer3 (RecP) 
  TH1F*   fRecPLadLay4;     //! histo with ladder occupancy on layer4 (RecP)
  TH1F*   fTrackPLadLay3;   //! histo with ladder occupancy on layer3 (TrP)
  TH1F*   fTrackPLadLay4;   //! histo with ladder occupancy on layer4 (TrP)
  TH1F*   fGoodAnLadLay3;   //! histo good anodes per ladder on layer3 
  TH1F*   fGoodAnLadLay4;   //! histo good anodes per ladder on layer4 
  TH1F*   fDriftTimeRP;     //! histo with drift time distribution (RecP)
  TH1F*   fDriftTimeTPAll;  //! histo with drift time distribution (TrP)
  TH1F*   fDriftTimeTPNoExtra; //! histo with drift time distribution (TrP)
  TH1F*   fDriftTimeTPExtra;//! histo with drift time distribution (TrP)
  TH1F*   fSignalTime[8];   //! histos of dE/dx in time windows
  AliESDEvent  *fESD;       // ESD object
  AliESDfriend *fESDfriend; // ESD friend object
  AliITSresponseSDD* fResp; // ResponseSDD object
  Bool_t  fUseITSsaTracks;   // Flag for using standalone ITS tracs
  Int_t   fMinITSpts;       // Minimum number of ITS points per track
  Int_t   fMinTPCpts;       // Minimum number of TPC points per track
  Float_t fMinPfordEdx;     // Minimum momentum for dE/dx
  Bool_t  fOnlyCINT1BTrig;  // Flag for using all events or only intections
  Bool_t  fInitialised;     // True if initialised
  Bool_t  fExcludeBadMod;   // Flag to reject bad modules
 
  ClassDef(AliAnalysisTaskSDDRP,3);  
};


#endif
