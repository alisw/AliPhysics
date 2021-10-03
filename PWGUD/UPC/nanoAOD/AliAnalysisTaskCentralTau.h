/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKCENTRALTAU_H
#define ALIANALYSISTASKCENTRALTAU_H

class TH1;
class TH2;
class TTree;
class TList;
class TFile;
class AliTOFTriggerMask;
class TBits;

#include "AliPIDResponse.h"
#include "AliESDtrackCuts.h"
#include "AliTimeRangeCut.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCentralTau : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskCentralTau();
  AliAnalysisTaskCentralTau(const char *name);
  virtual ~AliAnalysisTaskCentralTau();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
  void SetParameters(Float_t cutE){cutEta = cutE;}
  Int_t TestPIDhypothesis(AliESDtrack *trk);
  void TPCandTOFsignalInfo(AliESDtrack *trk, Int_t trkID);
  void SetCrossed(Int_t spd[4], TBits &crossed);
  Int_t GetChipId(Int_t index, Int_t &chipId2, Bool_t debug=0);
  Bool_t IsSTGFired(TBits bits, Int_t dphiMin=4, Int_t dphiMax=10, Bool_t tolerance = 1);
  void FillTree(TTree *t, TLorentzVector v);
 private:
 
  AliPIDResponse *fPIDResponse;
  AliESDtrackCuts *fTrackCutsBit0;
  AliESDtrackCuts *fTrackCutsBit1;
  AliESDtrackCuts *fTrackCutsBit4;
  Float_t cutEta;

  TList *fOutputList;		//<
  TList *fOutputPID;   //<
  TH2I *hTriggerCounter;	//!
  TH1I *hParticleTypeCounter; //
  TTree *tTwoTracks;		//!
  TTree *tPID;    //!
  
  Float_t fPtDaughter[2], fPt, fY, fM, fPhi, fZNAenergy, fZNCenergy, fZNAtime[4], fZNCtime[4];
  TLorentzVector fVectDaughter[2];
  Int_t fSignDaughter[2], fChannel, fSign, fRunNumber, fADAdecision, fADCdecision, fV0Adecision, fV0Cdecision;
  Bool_t fTriggers[10], fTriggerClass[3];
  
  // PID analysis
  Float_t fPIDpt[2],fTPCsignal[2],fTOFsignal[2];
  Int_t fTPCmostProbableTrackType[2], fTOFmostProbableTrackType[2];

  AliAnalysisTaskCentralTau(const AliAnalysisTaskCentralTau&); //not implemented
  AliAnalysisTaskCentralTau& operator =(const AliAnalysisTaskCentralTau&); //not implemented
  
  ClassDef(AliAnalysisTaskCentralTau, 25);
};

#endif
