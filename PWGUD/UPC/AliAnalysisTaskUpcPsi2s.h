/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKUPCPSI2S_H
#define ALIANALYSISTASKUPCPSI2S_H

class TClonesArray;
class TTree;
class TH1;
class TH2;
class TList;
class AliPIDResponse;
class AliAODEvent;
class AliESDEvent;

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
  virtual void RunAODtrig();
  virtual void RunAODhist();
  virtual void RunAODtree();
  virtual void RunAODMC(AliAODEvent *aod);
  virtual void RunAODsystematics(AliAODEvent *aod);
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
  Double_t GetMedian(Double_t *daArray);

 private:
  Int_t fType; // 0 - ESD, 1 - AOD
  Bool_t isMC;
  Bool_t fRunTree; 
  Bool_t fRunHist;
  Bool_t fRunSystematics;
  
  AliPIDResponse *fPIDResponse;
  
  //event tree
  TTree *fJPsiTree;
  TTree *fPsi2sTree;
  //tree variables
  Int_t fRunNum;
  UInt_t fPerNum, fOrbNum;
  //trigger
  Bool_t fTrigger[ntrg];
  Bool_t fTriggerInputsMC[4];
  UInt_t fL0inputs, fL1inputs;
  Bool_t fTOFtrig1, fTOFtrig2;
  Double_t fTOFphi[4];
  
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
  Double_t fVtxErr[3];
  Double_t fVtxChi2,fVtxNDF;
  Double_t fKfVtxPos[3];
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
  //tracks
  TClonesArray *fJPsiAODTracks;
  TClonesArray *fJPsiESDTracks; 
  TClonesArray *fPsi2sAODTracks;
  TClonesArray *fPsi2sESDTracks;
    //mc
  TClonesArray *fGenPart;
  
  TList *fListTrig;
  TH1D *fHistCcup4TriggersPerRun;
  TH1D *fHistCcup7TriggersPerRun;
  TH1D *fHistCcup2TriggersPerRun;
  TH1D *fHistCint1TriggersPerRun;
  TH1D *fHistCint6TriggersPerRun;
  TH1D *fHistC0tvxAndCint1TriggersPerRun;
  TH1D *fHistZedTriggersPerRun;
  TH1D *fHistCvlnTriggersPerRun;
  TH1D *fHistMBTriggersPerRun;
  TH1D *fHistCentralTriggersPerRun;
  TH1D *fHistSemiCentralTriggersPerRun;
  
  TList *fListHist;
  TH1D *fHistNeventsJPsi; 
  TH2D *fHistTPCsignalJPsi;
  TH2D *fHistDiLeptonPtJPsi;
  TH1D *fHistDiElectronMass;
  TH1D *fHistDiMuonMass;
  TH1D *fHistDiLeptonMass;
  
  TH1D *fHistNeventsPsi2s;
  TH2D *fHistPsi2sMassVsPt;
  TH1D *fHistPsi2sMassCoherent;
  
  TList *fListSystematics;
  TList *fListJPsiLoose;
  TList *fListJPsiTight;
  TList *fListPsi2sLoose;
  TList *fListPsi2sTight;
  
  AliAnalysisTaskUpcPsi2s(const AliAnalysisTaskUpcPsi2s&); //not implemented
  AliAnalysisTaskUpcPsi2s& operator =(const AliAnalysisTaskUpcPsi2s&); //not implemented
  
  ClassDef(AliAnalysisTaskUpcPsi2s, 1); 
};

#endif









