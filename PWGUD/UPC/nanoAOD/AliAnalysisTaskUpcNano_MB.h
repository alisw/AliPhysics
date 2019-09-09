/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKUPCNano_MB_H
#define ALIANALYSISTASKUPCNano_MB_H

class TH1;
class TH2;
class TTree;
class TList;
class TFile;
class AliTOFTriggerMask;
class TBits;

#include "AliTimeRangeCut.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskUpcNano_MB : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskUpcNano_MB();
  AliAnalysisTaskUpcNano_MB(const char *name);
  virtual ~AliAnalysisTaskUpcNano_MB();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void RunMC(AliVEvent *aod);
  virtual void Terminate(Option_t *);
  
  void SetIsMC(Bool_t MC){isMC = MC;}
  void SetIsESD(Bool_t ESD){isESD = ESD;}
  void SetCutEta(Float_t cut){cutEta = cut;}
  Double_t GetMedian(Double_t *daArray);
  void SetCrossed(Int_t spd[4], TBits &crossed);
  Int_t GetChipId(Int_t index, Int_t &chipId2, Bool_t debug=0);
  Bool_t IsSTGFired(TBits bits, Int_t dphiMin=4, Int_t dphiMax=10, Bool_t tolerance = 1);
  void FillTree(TTree *t, TLorentzVector v);
 private:
 
  AliPIDResponse *fPIDResponse;
  AliTimeRangeCut fTimeRangeCut;
  AliESDtrackCuts *fTrackCutsBit0;
  AliESDtrackCuts *fTrackCutsBit1;
  AliESDtrackCuts *fTrackCutsBit5;
  Bool_t isMC; 
  Bool_t isESD;
  Float_t cutEta;

  TList *fOutputList;		//<
  TH1D *fHistEvents;		//!
  TH1D *fHistMCTriggers;	//!
   
  TTree *fTreePhi;		//!
  TTree *fTreeJPsi;		//!
  TTree *fTreePsi2s;		//!
  TTree *fTreeRho;		//!
  TTree *fTreeGen;		//!

  TH2D *hTPCPIDMuonCorr; 	//!
  TH1D *hTPCPIDMuon;		//!
  TH2D *hTPCPIDElectronCorr;	//!
  TH1D *hTPCPIDElectron;	//!
  TH1D *hTPCPIDPion;		//!
  TH2D *hTPCPIDPionCorr;	//!
  TH1D *hTOFPIDProton;		//!
  TH2D *hTOFPIDProtonCorr;	//!
  TH1D *hITSPIDKaon;		//!
  TH2D *hITSPIDKaonCorr;	//!
  TH2D *hTPCdEdxCorr;		//!
  TH2I *hTriggerCounter;	//!
  
  Float_t fPt, fY, fM, fDiLeptonM, fDiLeptonPt, fZNAenergy, fZNCenergy, fZNAtime[4], fZNCtime[4], fPIDsigma;
  Int_t fChannel, fSign, fRunNumber;
  Bool_t fTriggerInputsMC[11], fTriggers[10], fInEtaGen, fInEtaRec;
  
  TFile *fSPDfile;
  TFile *fTOFfile;
  Int_t fLoadedRun;
  TH2F *hTOFeff;
  TH1D *hSPDeff;
  AliTOFTriggerMask *fTOFmask;
  TBits fFOCrossFiredChips;
  
  AliAnalysisTaskUpcNano_MB(const AliAnalysisTaskUpcNano_MB&); //not implemented
  AliAnalysisTaskUpcNano_MB& operator =(const AliAnalysisTaskUpcNano_MB&); //not implemented
  
  ClassDef(AliAnalysisTaskUpcNano_MB, 26); 
};

#endif
