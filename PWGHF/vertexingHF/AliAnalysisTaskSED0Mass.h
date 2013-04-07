#ifndef ALIANALYSISTASKSED0MASS_H
#define ALIANALYSISTASKSED0MASS_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// Class AliAnalysisTaskSED0Mass
// AliAnalysisTaskSE for D0 candidates invariant mass histogram
// and comparison to MC truth (kinematics stored in the AOD) and cut variables
// distributions
// Authors: A.Dainese, andrea.dainese@ln.infn.it
// and C.Bianchin, chiara.bianchin@pd.infn.it
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1F.h>
#include <THnSparse.h>

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliNormalizationCounter.h"

class AliAODEvent;

class AliAnalysisTaskSED0Mass : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSED0Mass();
  AliAnalysisTaskSED0Mass(const char *name,AliRDHFCutsD0toKpi* cuts);
  virtual ~AliAnalysisTaskSED0Mass();


  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void SetArray(Int_t type=AliAnalysisTaskSED0Mass::kD0){fArray=type;}
  enum{kD0,kLS};

  void SetReadMC(Bool_t readMC=kFALSE){fReadMC=readMC;}
  void SetCutOnDistr(Bool_t cutondistr=kFALSE){fCutOnDistr=cutondistr;}
  void SetUsePid4Distr(Bool_t usepid=kTRUE){fUsePid4Distr=usepid;}
  void SetFillOnlyD0D0bar(Int_t flagfill){fFillOnlyD0D0bar=flagfill;}
  void SetFillVarHists(Bool_t flag) {fFillVarHists=flag;}
  void SetFillPtHistos(Bool_t flag) {fFillPtHist=flag;}
  void SetFillImpactParameterHistos(Bool_t flag) {fFillImpParHist=flag;}
  void SetSystem(Int_t sys){fSys=sys; if(fSys==1) SetFillVarHists(kFALSE);}
  void SetRejectSDDClusters(Bool_t flag) { fIsRejectSDDClusters=flag; }
  void SetUseSelectionBit(Bool_t flag) { fUseSelectionBit=flag; }
  void SetWriteVariableTree(Bool_t flag) { fWriteVariableTree=flag; }
  void SetDrawDetSignal(Bool_t flag) { fDrawDetSignal=flag; }
  void SetPIDCheck(Bool_t flag) { fPIDCheck=flag; }


  Bool_t GetCutOnDistr() const {return fCutOnDistr;}
  Bool_t GetUsePid4Distr() const {return fUsePid4Distr;}
  Int_t  GetFillOnlyD0D0bar() const {return fFillOnlyD0D0bar;}
  Bool_t GetFillVarHists() const {return fFillVarHists;}
  Bool_t GetFillPtHistos() const {return fFillPtHist;}
  Bool_t GetFillImpactParameterHistos() const {return fFillImpParHist;}
  Int_t  GetSystem() const {return fSys;}
  Bool_t GetRejectSDDClusters() const { return fIsRejectSDDClusters; }
  Bool_t GetUseSelectionBit() const { return fUseSelectionBit; }
  Bool_t GetWriteVariableTree() const {return fWriteVariableTree;}
  Bool_t GetDrawDetSignal() const {return fDrawDetSignal;}
  Bool_t GetPIDCheck() const {return fPIDCheck;}

 private:

  AliAnalysisTaskSED0Mass(const AliAnalysisTaskSED0Mass &source);
  AliAnalysisTaskSED0Mass& operator=(const AliAnalysisTaskSED0Mass& source); 
  void	   DrawDetSignal(AliAODRecoDecayHF2Prong *part, TList *ListDetSignal);

  void     FillMassHists(AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliAODMCHeader *mcHeader, AliRDHFCutsD0toKpi *cuts, TList *listout);
  void     FillVarHists(AliAODEvent *aodev,AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi *cuts, TList *listout);
  AliAODVertex* GetPrimaryVtxSkipped(AliAODEvent *aodev);
  void CreateImpactParameterHistos();
  Int_t CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const;
  Float_t GetTrueImpactParameter(AliAODMCHeader *mcHeader, TClonesArray* arrayMC, AliAODMCParticle *partD0) const ;

  TList    *fOutputMass;          //! list send on output slot 1
  TList    *fOutputMassPt;        //! list send on output slot 6
  TList    *fDistr;               //! list send on output slot 2
  TH1F     *fNentries;            //! histogram with number of events on output slot 3
  AliRDHFCutsD0toKpi *fCuts;      //  Cuts - sent to output slot 4
  THnSparseF *fHistMassPtImpParTC[5];   //! histograms for impact paramter studies
  Int_t     fArray;               //  can be D0 or Like Sign candidates
  Bool_t    fReadMC;              //  flag for MC array: kTRUE = read it, kFALSE = do not read it
  Bool_t    fCutOnDistr;          //  flag to decide if apply cut also on distributions: 0 no cuts, 1 looser cuts, 2 tighter cuts 
  Bool_t    fUsePid4Distr;        //  flag to use the particle identification to fill the signal histograms of distributions. It has effect only with fReadMC=kFALSE
  AliNormalizationCounter *fCounter;//!AliNormalizationCounter on output slot 5
  Int_t     fNPtBins;             //  number of pt bins
  Double_t  fLsNormalization;     //  normalization
  Int_t     fFillOnlyD0D0bar;     // flag to fill mass histogram with D0/D0bar only (0 = fill with both, 1 = fill with D0 only, 2 = fill with D0bar only)
  TObjArray fDaughterTracks;      // keeps the daughter tracks
  Int_t     fIsSelectedCandidate; // selection outcome
  Bool_t    fFillVarHists;        // flag to enable filling variable histos
  Int_t     fSys;                 // fSys=0 -> p-p; fSys=1 ->PbPb (in this case fFillVarHists=kFALSE by default: set it to kTRUE *after* if needed)
  Bool_t    fIsRejectSDDClusters; // flag to reject events with SDD clusters
  Bool_t    fFillPtHist;          // flag to fill Pt and Impact Parameter Histograms
  Bool_t    fFillImpParHist;      // flag to fill Pt and Impact Parameter Histograms
  Bool_t    fUseSelectionBit;     // flag to check or not the selection bit

  Bool_t    fWriteVariableTree;       // flag to decide whether to write the candidate variables on a tree variables
  TTree    *fVariablesTree;           //! tree of the candidate variables after track selection on output slot 7
  Double_t *fCandidateVariables;      //!  variables to be written to the tree
  Bool_t	fPIDCheck;			// flag to decide whether to fill "PID = x" bins in fNentrie
  Bool_t    fDrawDetSignal;		// flag to decide whether to draw the TPC dE/dx and TOF signal before/after PID
  TList	   *fDetSignal;		//!Detector signal histograms (on output slot 8)

  ClassDef(AliAnalysisTaskSED0Mass,18); // AliAnalysisTaskSE for D0->Kpi
};

#endif

