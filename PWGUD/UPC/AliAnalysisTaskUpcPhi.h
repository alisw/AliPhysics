/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKUPCPHI_H
#define ALIANALYSISTASKUPCPHI_H

class TClonesArray;
class TTree;
class TH1;
class TH2;
class TList;
class AliPIDResponse;
class AliAODEvent;
class AliESDEvent;
class TBits;

#define ntrg 17
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskUpcPhi : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskUpcPhi();
  AliAnalysisTaskUpcPhi(const char *name);
  virtual ~AliAnalysisTaskUpcPhi();

  virtual void Init();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void RunAODtrig();
  virtual void RunAODhist();
  virtual void RunAODtree();
  virtual void RunAODMC(AliAODEvent *aod);
  virtual void RunESDtrig();
  virtual void RunESDhist();
  virtual void RunESDtree();
  virtual void RunESDMC(AliESDEvent *esd);
  virtual void Terminate(Option_t *);
  void SetRunTree(Bool_t runTree){fRunTree = runTree;}
  void SetRunHist(Bool_t runHist){fRunHist = runHist;}
  void SetRunSyst(Bool_t runSyst){fRunSystematics = runSyst;}
  void SetIsMC(Bool_t MC){isMC = MC;}
  void InitSystematics();

 private:
  Int_t fType; // 0 - ESD, 1 - AOD
  Bool_t isMC;
  Bool_t fRunTree; 
  Bool_t fRunHist;
  Bool_t fRunSystematics;
  
  AliPIDResponse *fPIDResponse;
  
  //event tree
  TTree *fITSTree;
  TTree *fTPCTree;
  
  //tree variables
  Int_t fRunNum;
  UInt_t fPerNum, fOrbNum;
  //trigger
  Bool_t fTrigger[ntrg];
  Bool_t fTriggerInputsMC[ntrg];
  UInt_t fL0inputs, fL1inputs;
  
  Double_t fPIDITSMuon[2];
  Double_t fPIDITSElectron[2];
  Double_t fPIDITSPion[2];
  Double_t fPIDITSKaon[2];
  Double_t fPIDITSProton[2];
  
  Double_t fPIDTPCMuon[2];
  Double_t fPIDTPCElectron[2];
  Double_t fPIDTPCPion[2];
  Double_t fPIDTPCKaon[2];
  Double_t fPIDTPCProton[2];
  
  Int_t fVtxContrib;
  Double_t fVtxPos[3];
  Double_t fVtxErr[3];
  Double_t fVtxChi2,fVtxNDF;
  Double_t fKfVtxPos[3];
  Double_t fSpdVtxPos[3];
  TBits fFastOrMap;
  UShort_t fBCrossNum, fNtracklets, fNLooseITSTracks, fNLooseTPCTracks;
  //vzero, zdc
  Double_t fZDCAenergy, fZDCCenergy;
  Int_t fV0Adecision, fV0Cdecision;
  Int_t fADAdecision, fADCdecision;
  //input data
  TObjString *fDataFilnam;
  Short_t fRecoPass;
  Long64_t fEvtNum;
  //tracks
  TClonesArray *fPhiAODTracks;
  TClonesArray *fPhiESDTracks;
    //mc
  TClonesArray *fGenPart;
  
  TList *fListTrig;
  TH1D *fHistCcup4TriggersPerRun;
  TH1D *fHistCcup7TriggersPerRun;
  TH1D *fHistCcup2TriggersPerRun;
  TH1D *fHistCint1TriggersPerRun;
  TH1D *fHistC0tvxAndCint1TriggersPerRun;
  TH1D *fHistZedTriggersPerRun;
  TH1D *fHistCvlnTriggersPerRun;
  TH1D *fHistMBTriggersPerRun;
  TH1D *fHistCentralTriggersPerRun;
  TH1D *fHistSemiCentralTriggersPerRun;
  
  TH1D *fHistCTest58TriggersPerRun;
  TH1D *fHistCTest59TriggersPerRun;
  TH1D *fHistCTest60TriggersPerRun;
  TH1D *fHistCTest61TriggersPerRun;
  
  TH1D *fHistCcup8TriggersPerRun;
  TH1D *fHistCcup9TriggersPerRun;
  
  TList *fListHist;
  
  AliAnalysisTaskUpcPhi(const AliAnalysisTaskUpcPhi&); //not implemented
  AliAnalysisTaskUpcPhi& operator =(const AliAnalysisTaskUpcPhi&); //not implemented
  
  ClassDef(AliAnalysisTaskUpcPhi, 3); 
};

#endif









