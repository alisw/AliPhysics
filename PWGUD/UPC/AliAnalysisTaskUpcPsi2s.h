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
class AliTOFTriggerMask;
class TBits;
class TFile;

#define ntrgMB 20
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
  void SetTracking(Int_t tracking){fTracking = tracking;}
  void SetIsMC(Bool_t MC){isMC = MC;}
  void InitSystematics();
  Double_t GetMedian(Double_t *daArray);

 private:
  Int_t fType; // 0 - ESD, 1 - AOD
  Int_t fTracking; //0 - Global, 1 - ITSsa
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
  Bool_t fTrigger[ntrgMB];
  Bool_t fTriggerInputsMC[ntrgMB];
  UInt_t fL0inputs, fL1inputs;
  AliTOFTriggerMask *fTOFmask;
  Bool_t fIsPhysicsSelected;
  
  Float_t fPIDTPCMuon[4];
  Float_t fPIDTPCElectron[4];
  Float_t fPIDTPCPion[4];
  Float_t fPIDTPCKaon[4];
  Float_t fPIDTPCProton[4];
  
  Float_t fPIDTOFMuon[4];
  Float_t fPIDTOFElectron[4];
  Float_t fPIDTOFPion[4];
  Float_t fPIDTOFKaon[4];
  Float_t fPIDTOFProton[4];
  
  Int_t fVtxContrib;
  Float_t fVtxPos[3];
  Float_t fMCVtxPos[3];
  Float_t fVtxErr[3];
  Float_t fVtxChi2,fVtxNDF;
  Float_t fKfVtxPos[3];
  Int_t fSpdVtxContrib;
  Float_t fSpdVtxPos[3];
  
  Bool_t fIsVtxContributor[4];
  
  UShort_t fBCrossNum, fNtracklets, fNLooseTracks;
  //vzero, zdc
  Float_t fZNAenergy, fZNCenergy;
  Float_t fZPAenergy, fZPCenergy;
  Float_t fZNATDCm[4];
  Float_t fZNCTDCm[4];
  Float_t fZPATDCm[4];
  Float_t fZPCTDCm[4];
  Int_t fV0Adecision, fV0Cdecision;
  Int_t fADAdecision, fADCdecision;
  //input data
  TObjString *fDataFilnam;
  Short_t fRecoPass;
  Long64_t fEvtNum;
  //spd
  TBits fFOFiredChips;
  //PF protection
  TBits fIR1Map;
  TBits fIR2Map;
  //tracks
  TClonesArray *fJPsiAODTracks;
  TClonesArray *fJPsiESDTracks; 
  TClonesArray *fPsi2sAODTracks;
  TClonesArray *fPsi2sESDTracks;
    //mc
  TClonesArray *fGenPart;
  
  //EVE tree
  TTree *fEveTree;
  Float_t fPt, fY, fM, fDiLeptonM, fDiLeptonPt, fPIDsigma;
  Int_t fChannel;
  
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

  TH1D *fHistCcup8TriggersPerRun;
  TH1D *fHistCcup9TriggersPerRun;
  TH1D *fHistCcup10TriggersPerRun;
  TH1D *fHistCcup11TriggersPerRun;
  TH1D *fHistCcup12TriggersPerRun;
  TH1D *fHistCcup25TriggersPerRun;
  TH1D *fHistCcup26TriggersPerRun;
  TH1D *fHistCcup27TriggersPerRun;
  TH1D *fHistCcup29TriggersPerRun;
  TH1D *fHistCcup30TriggersPerRun;
  TH1D *fHistCcup31TriggersPerRun;
  TH1D *fHistCtrueTriggersPerRun;
  
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
  
  TFile *fSPDfile;
  TH1D *hBCmod4;
  TH2D *hSPDeff;
  
  AliAnalysisTaskUpcPsi2s(const AliAnalysisTaskUpcPsi2s&); //not implemented
  AliAnalysisTaskUpcPsi2s& operator =(const AliAnalysisTaskUpcPsi2s&); //not implemented
  
  ClassDef(AliAnalysisTaskUpcPsi2s, 11); 
};

#endif









