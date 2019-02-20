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
class TProfile;
class TString;
class AliESDEvent;
class AliESDfriend;
class AliITSresponseSDD;
class AliTriggerConfiguration;

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
  void SetTriggerClass(TString trclass) {
    fTriggerClass=trclass;
  }
  void SetUseOnlyEventsWithSDD(Bool_t use=kTRUE){
    fOnlyEventsWithSDD=use;
  }
  void SetMinPfordEdx(Float_t minp=0.5){
    fMinPfordEdx=minp;
  }
  void SetExcludeBadModules(Bool_t opt=kTRUE){
    fExcludeBadMod=opt;
  }
  void SetReadCDB(Bool_t opt=kTRUE){
    fReadCDB=opt;
  }
  Bool_t CheckModule(Int_t lay, Int_t lad, Int_t det) const;
 

 private:
  AliAnalysisTaskSDDRP(const AliAnalysisTaskSDDRP &source);
  AliAnalysisTaskSDDRP& operator=(const AliAnalysisTaskSDDRP &source);
  
  TList*  fOutput;          //! QA histos
  TH1F*   fHistNEvents;     //! histo with N of events  
  TH1F*   fHistCluInLay;    //! histo with number of tracks per layer
  TH1F*   fHistAllPMod;     //! histo of tracks crossing SDD modules
  TH1F*   fHistGoodPMod;    //! histo of tracks with good point in SDD module
  TH1F*   fHistBadRegMod;   //! histo of tracks crossing bad region of SDD mod.
  TH1F*   fHistMissPMod;    //! histo of tracks with missing point in SDD mod.
  TH1F*   fHistSkippedMod;  //! histo of tracks skipping an SDD module
  TH1F*   fHistOutAccMod;   //! histo of tracks out of accept. in SDD module
  TH1F*   fHistNoRefitMod;  //! histo of points rejected in refit vs. SDD mod.

  TH1F*   fHistAllPXloc;    //! histo of xlocal for track hit points
  TH1F*   fHistGoodPXloc;   //! histo of xlocal for track hit points + good clu
  TH1F*   fHistBadRegXloc;  //! histo of xlocal for track hit points + bad reg.
  TH1F*   fHistMissPXloc;   //! histo of xlocal for track hit points + miss clu
  TH1F*   fHistAllPZloc;    //! histo of zlocal for track hit points
  TH1F*   fHistGoodPZloc;   //! histo of zlocal for track hit points + good clu
  TH1F*   fHistBadRegZloc;  //! histo of zlocal for track hit points + bad reg.
  TH1F*   fHistMissPZloc;   //! histo of zlocal for track hit points + miss clu

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
  TH2F*   fEtaPhiTracks;    //! eta,phi distribution of tracks
  TH2F*   fEtaPhiTracksLay3; //! eta,phi distrib. of tracks with point on Lay3
  TH2F*   fEtaPhiTracksLay4; //! eta,phi distrib. of tracks with point on Lay4
  TH1F*   fDriftTimeRP;     //! histo with drift time distribution (RecP)
  TH1F*   fDriftTimeTPAll;  //! histo with drift time distribution (TrP)
  TH2F*   fDriftTimeTPAllMod; //! histo with drift time distribution (TrP) per module
  TH1F*   fDriftTimeTPNoExtra; //! histo with drift time distribution (TrP)
  TH1F*   fDriftTimeTPExtra;//! histo with drift time distribution (TrP)
  TH1F*   fSignalTime[8];   //! histos of dE/dx in time windows
  TH2F*   fCluSizAnVsTime;  //! Histo with anode cluster size vs. time
  TH2F*   fCluSizTbVsTime;  //! Histo with time-bin cluster size vs. time
  TProfile* fProfRecPtsLay3VsTime; //! Profile of occupancy vs. time
  TProfile* fProfRecPtsLay4VsTime; //! Profile of occupancy vs. time
  TProfile* fProfTrPtsLay3VsTime;  //! Profile of occupancy vs. time
  TProfile* fProfTrPtsLay4VsTime;  //! Profile of occupancy vs. time
  TProfile* fProfFracTrRecLay3VsTime;  //! Profile of occupancy vs. time
  TProfile* fProfFracTrRecLay4VsTime;  //! Profile of occupancy vs. time
  TProfile* fProfFracTrkWithPntLay3VsTime;  //! Profile of occupancy vs. time
  TProfile* fProfFracTrkWithPntLay4VsTime;  //! Profile of occupancy vs. time
  AliITSresponseSDD* fResp; // ResponseSDD object
  AliTriggerConfiguration* fTrigConfig; // trigger configuration object
  Bool_t  fUseITSsaTracks;   // Flag for using standalone ITS tracs
  Int_t   fMinITSpts;       // Minimum number of ITS points per track
  Int_t   fMinTPCpts;       // Minimum number of TPC points per track
  Float_t fMinPfordEdx;     // Minimum momentum for dE/dx
  TString fTriggerClass;    // Name of selected trigger class
  Bool_t  fOnlyEventsWithSDD; // Flag to use only trigger cluster with SDD
  Bool_t  fExcludeBadMod;   // Flag to reject bad modules
  Bool_t  fReadCDB;         // Flag to Switch on/off OCDB access
  Bool_t  fInitCalib;       // Flag to check calib initiatization
 
  ClassDef(AliAnalysisTaskSDDRP,9);
};


#endif
