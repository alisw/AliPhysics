/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKUPCPSI2S_H
#define ALIANALYSISTASKUPCPSI2S_H

class TClonesArray;
class TTree;
class TH1;

#define ntrg 17
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskUpcPsi2s : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskUpcPsi2s();
  AliAnalysisTaskUpcPsi2s(const char *name);
  virtual ~AliAnalysisTaskUpcPsi2s();

  virtual void Init();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void RunAOD();
  virtual void RunESD();
  virtual void Terminate(Option_t *);

 private:
 
  Int_t fType; // 0 - ESD, 1 - AOD
  
  //counter
  TH1I *hCounter; //!

  //event tree
  TTree *fJPsiTree;
  TTree *fPsi2sTree;
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
  TObjString *fDataFilnam;
  Short_t fRecoPass;
  Long64_t fEvtNum;
  //tracks
  TClonesArray *fJPsiAODTracks;
  TClonesArray *fJPsiESDTracks; 
  TClonesArray *fPsi2sAODTracks;
  TClonesArray *fPsi2sESDTracks;
  
  AliAnalysisTaskUpcPsi2s(const AliAnalysisTaskUpcPsi2s&); //not implemented
  AliAnalysisTaskUpcPsi2s& operator =(const AliAnalysisTaskUpcPsi2s&); //not implemented
  
  ClassDef(AliAnalysisTaskUpcPsi2s, 1); 
};

#endif









