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
class AliPIDResponse;

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
  virtual void RunAOD();
  virtual void RunESDtrig();
  virtual void RunESDtree();
  virtual void Terminate(Option_t *);

 private:
  Int_t fType; // 0 - ESD, 1 - AOD
  
  AliPIDResponse *fPIDResponse;
  
  //event tree
  TTree *fK0sTree;
  //tree variables
  Int_t fRunNum;
  UInt_t fPerNum, fOrbNum;
  //trigger
  Bool_t fTrigger[ntrg];
  Bool_t fTriggerInputsMC[ntrg];
  UInt_t fL0inputs, fL1inputs;
  
  Double_t fPIDTPCMuon[4];
  Double_t fPIDTPCElectron[4];
  Double_t fPIDTPCPion[4];
  Double_t fPIDTPCKaon[4];
  Double_t fPIDTPCProton[4];
  
  Double_t fPIDTOFMuon[4];
  Double_t fPIDTOFElectron[4];
  Double_t fPIDTOFPion[4];
  Double_t fPIDTOFKaon[4];
  Double_t fPIDTOFProton[4];
  
  Int_t fVtxContrib;
  Double_t fVtxPos[3];
  Double_t fMCVtxPos[3];
  Double_t fVtxErr[3];
  Double_t fVtxChi2,fVtxNDF;
  Int_t fSpdVtxContrib;
  Double_t fSpdVtxPos[3];
  
  Bool_t fIsVtxContributor[4];
  
  UShort_t fBCrossNum, fNtracklets, fNLooseTracks;
  //vzero, zdc
  Double_t fZNAenergy, fZNCenergy;
  Double_t fZPAenergy, fZPCenergy;
  Double_t fZDCAtime, fZDCCtime;
  Int_t fV0Adecision, fV0Cdecision;
  //input data
  TObjString *fDataFilnam;
  Short_t fRecoPass;
  Long64_t fEvtNum;
  //vertices
  TClonesArray *fK0sAODv0s;
  //tracks
  TClonesArray *fK0sAODTracks;
  TClonesArray *fK0sESDTracks; 
  
  TList *fListTrig;
  TH1D *fHistCcup4TriggersPerRun;
  TH1D *fHistCcup7TriggersPerRun;
  TH1D *fHistCcup2TriggersPerRun;
  TH1D *fHistCint1TriggersPerRun;
  TH1D *fHistCint6TriggersPerRun;
  TH1D *fHistC0tvxAndCint1TriggersPerRun;
  
  AliAnalysisTaskUpcK0sK0s(const AliAnalysisTaskUpcK0sK0s&); //not implemented
  AliAnalysisTaskUpcK0sK0s& operator =(const AliAnalysisTaskUpcK0sK0s&); //not implemented
  
  ClassDef(AliAnalysisTaskUpcK0sK0s, 1); 
};

#endif









