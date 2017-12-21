/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKUPCNano_MB_H
#define ALIANALYSISTASKUPCNano_MB_H

class TH1;
class TTree;
class TList;
class TBits;
class AliTOFTriggerMask;
class AliTriggerBCMask;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskUpcNano_MB : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskUpcNano_MB();
  AliAnalysisTaskUpcNano_MB(const char *name);
  virtual ~AliAnalysisTaskUpcNano_MB();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void RunMC(AliAODEvent *aod);
  virtual void Terminate(Option_t *);
  
  void SetIsMC(Bool_t MC){isMC = MC;}
  void SetCutEta(Bool_t toCut){cutEta = toCut;}
  Double_t GetMedian(Double_t *daArray);
  void FillTree(TTree *t, TLorentzVector v);
 private:
 
  AliPIDResponse *fPIDResponse;
  Bool_t isMC, cutEta;

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
  TH1D *hNLooseTracks;		//!
  TH1D *hNV0BB;			//!
  TH1D *hNV0BG;			//!
  TH1D *hNADBB;			//!
  TH1D *hNADBG;			//!
  
  Double_t fPt, fY, fM, fDiLeptonM, fDiLeptonPt, fZNAenergy, fZNCenergy, fZNAtime, fZNCtime, fPIDsigma;
  Int_t fChannel, fSign, fRunNumber, fOldRun, fClosestIR1, fClosestIR2,fNV0BB,fNV0BG,fNADBB,fNADBG,fGoodBC;
  Bool_t fTriggerInputsMC[10];
  TBits fFOFiredChips;
  AliTriggerBCMask *fBCmask;
  AliTOFTriggerMask *fTOFmask;
  
  AliAnalysisTaskUpcNano_MB(const AliAnalysisTaskUpcNano_MB&); //not implemented
  AliAnalysisTaskUpcNano_MB& operator =(const AliAnalysisTaskUpcNano_MB&); //not implemented
  
  ClassDef(AliAnalysisTaskUpcNano_MB, 9); 
};

#endif
