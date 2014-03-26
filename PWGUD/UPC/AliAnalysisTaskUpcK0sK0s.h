/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKUPCK0SK0S_H
#define ALIANALYSISTASKUPCK0SK0S_H

class TClonesArray;
class TTree;
class TH1;
class TH2;
class TList;

#define ntrg 17
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskUpcK0sK0s : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskUpcK0sK0s();
  AliAnalysisTaskUpcK0sK0s(const char *name);
  virtual ~AliAnalysisTaskUpcK0sK0s();

  virtual void Init();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void RunAODhist();
  virtual void RunAODtree();
  virtual void RunESD();
  virtual void Terminate(Option_t *);
  void SetRunTree(Bool_t runTree){fRunTree = runTree;}
  void SetRunHist(Bool_t runHist){fRunHist = runHist;}
  void SortArray(Int_t *dArray);

 private:
  Int_t fType; // 0 - ESD, 1 - AOD
  Bool_t fRunTree; 
  Bool_t fRunHist;
  
  //event tree
  TTree *fK0sTree;
  //tree variables
  Int_t fRunNum;
  UInt_t fPerNum, fOrbNum;
  //trigger
  Bool_t fTrigger[ntrg];
  UInt_t fL0inputs, fL1inputs;
  Int_t fVtxContrib;
  UShort_t fBCrossNum, fNtracklets;
  //vzero, zdc
  Double_t fZDCAenergy, fZDCCenergy;
  Int_t fV0Adecision, fV0Cdecision;
  //input data
  TString fDataFilnam;
  Short_t fRecoPass;
  Long64_t fEvtNum;
  //vertices
  TClonesArray *fK0sAODv0s;
  TClonesArray *fK0sESDv0s;
  //tracks
  TClonesArray *fK0sAODTracks;
  TClonesArray *fK0sESDTracks; 
 
  
  TList *fListHist;
  TH1D *fHistTriggersPerRun;
  TH2D *fHistK0sMass;
  
  AliAnalysisTaskUpcK0sK0s(const AliAnalysisTaskUpcK0sK0s&); //not implemented
  AliAnalysisTaskUpcK0sK0s& operator =(const AliAnalysisTaskUpcK0sK0s&); //not implemented
  
  ClassDef(AliAnalysisTaskUpcK0sK0s, 1); 
};

#endif









