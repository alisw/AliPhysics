#ifndef AliAnalysisTaskSED0MassTMVA_H
#define AliAnalysisTaskSED0MassTMVA_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///*************************************************************************
/// \class Class AliAnalysisTaskSED0MassTMVA
/// \brief AliAnalysisTaskSE for D0 candidates invariant mass histogram
/// and comparison to MC truth (kinematics stored in the AOD) and cut variables
/// distributions
/// \author Authors: A.Dainese, andrea.dainese@ln.infn.it
/// \author and C.Bianchin, chiara.bianchin@pd.infn.it
///*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1F.h>
#include <THnSparse.h>

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliNormalizationCounter.h"
#include "AliEventCuts.h"

class AliAODEvent;

class AliAnalysisTaskSED0MassTMVA : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSED0MassTMVA();
  AliAnalysisTaskSED0MassTMVA(const char *name,AliRDHFCutsD0toKpi* cuts);
  virtual ~AliAnalysisTaskSED0MassTMVA();


  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void CreateMCAcceptanceHistos();
  Bool_t CheckAcc(TClonesArray* arrayMC,Int_t nProng, Int_t *labDau);
  void FillMCAcceptanceHistos(TClonesArray *arrayMC, AliAODMCHeader *mcHeader);

  void NormIPvar(AliAODEvent *aod, AliAODRecoDecayHF2Prong *part,TClonesArray *arrMC);
  void SetArray(Int_t type=AliAnalysisTaskSED0MassTMVA::kD0){fArray=type;}
  enum{kD0,kLS};

  void SetReadMC(Bool_t readMC=kFALSE){fReadMC=readMC;}
  void SetDoMCAcceptanceHistos(Bool_t doMCAcc=kTRUE){fStepMCAcc=doMCAcc;}
  void SetCutOnDistr(Bool_t cutondistr=kFALSE){fCutOnDistr=cutondistr;}
  void SetUsePid4Distr(Bool_t usepid=kTRUE){fUsePid4Distr=usepid;}
  void SetFillOnlyD0D0bar(Int_t flagfill){fFillOnlyD0D0bar=flagfill;}
  void SetFillVarHists(Bool_t flag) {fFillVarHists=flag;}
  void SetFillPtHistos(Bool_t flag) {fFillPtHist=flag;}
  void SetFillYHistos(Bool_t flag) {fFillYHist=flag;}
  void SetFillImpactParameterHistos(Bool_t flag) {fFillImpParHist=flag;}
  void SetFillSparses(Bool_t flag) {fFillSparses=flag;}
  void SetUseRejectionMethod(Bool_t flag=kFALSE, Float_t factor=0.01) {fUseRejectionMethod=flag; fRejectionFactor=factor;}
  void SetSystem(Int_t sys){fSys=sys; if(fSys==1) SetFillVarHists(kFALSE);}
  void SetRejectSDDClusters(Bool_t flag) { fIsRejectSDDClusters=flag; }
  void SetUseSelectionBit(Bool_t flag) { fUseSelectionBit=flag; }
  void SetWriteVariableTree(Bool_t flag) { fWriteVariableTree=flag; }
  void SetWriteProtosgnVar(Bool_t flag) { fWriteProtosgnVar=flag; }
  void SetSelectTrueD0(Bool_t flag) { fSelectTrueD0 = flag; }
  void SetUseMassWindow(Bool_t flag) { fUsedMassWindow=flag; }
  void SetDrawDetSignal(Bool_t flag) { fDrawDetSignal=flag; }
  void SetPIDCheck(Bool_t flag) { fPIDCheck=flag; }
  void SetUseQuarkLevelTag(Bool_t opt){fUseQuarkTagInKine=opt;}
  void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}
  void SetPileupRejectionVZEROTPCout(Bool_t flag) {fEnablePileupRejVZEROTPCout=flag;}
  void SetPileupRejectionVZEROTPCcls(Bool_t flag, Bool_t rejpileup) {fEnablePileupRejVZEROTPCcls=flag; fRejectOutOfBunchPileUp=rejpileup;}
  void SetFillSubSampleHist(Bool_t flag) {fFillSubSampleHist=flag;}


  void SetEnableCentralityCorrCutsPbPb(Bool_t flag=kFALSE, Int_t year=2018) {
    fEnableCentralityCorrCuts=flag;
    if(year==2018){
      fEventCuts.SetupPbPb2018();
      fEventCuts.SetManualMode();
    }else{
      fEventCuts.SetupRun2PbPb();
    }
  }


  Bool_t GetCutOnDistr() const {return fCutOnDistr;}
  Bool_t GetUsePid4Distr() const {return fUsePid4Distr;}
  Int_t  GetFillOnlyD0D0bar() const {return fFillOnlyD0D0bar;}
  Bool_t GetFillVarHists() const {return fFillVarHists;}
  Bool_t GetFillPtHistos() const {return fFillPtHist;}
  Bool_t GetFillYHistos() const {return fFillYHist;}
  Bool_t GetFillImpactParameterHistos() const {return fFillImpParHist;}
  Int_t  GetSystem() const {return fSys;}
  Bool_t GetRejectSDDClusters() const { return fIsRejectSDDClusters; }
  Bool_t GetUseSelectionBit() const { return fUseSelectionBit; }
  Bool_t GetWriteVariableTree() const {return fWriteVariableTree;}
  Bool_t GetDrawDetSignal() const {return fDrawDetSignal;}
  Bool_t GetPIDCheck() const {return fPIDCheck;}
  Bool_t GetFillSubSampleHist() const {return fFillSubSampleHist;}

 private:

  AliAnalysisTaskSED0MassTMVA(const AliAnalysisTaskSED0MassTMVA &source);
  AliAnalysisTaskSED0MassTMVA& operator=(const AliAnalysisTaskSED0MassTMVA& source);
  void	   DrawDetSignal(AliAODRecoDecayHF2Prong *part, TList *ListDetSignal);
  Double_t GetBeautyMotherY(TClonesArray* arrayMC, AliAODMCParticle *mcPart);
  void     FillMassHists(AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliAODMCHeader *mcHeader, AliRDHFCutsD0toKpi *cuts, TList *listout);
  void     FillVarHists(AliAODEvent *aodev,AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi *cuts, TList *listout);
  void     FillCandVariables(AliAODEvent *aodev, AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliAODMCHeader *mcHeader, AliRDHFCutsD0toKpi *cuts);
  AliAODVertex* GetPrimaryVtxSkipped(AliAODEvent *aodev);
  void CreateImpactParameterHistos();
  Int_t CheckOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const;
  Float_t GetTrueImpactParameter(AliAODMCHeader *mcHeader, TClonesArray* arrayMC, AliAODMCParticle *partD0) const ;
  Float_t ComputeTopomatic(AliAODEvent *aodev,AliAODRecoDecayHF2Prong *part);

  TList    *fOutputMass;          //!<! list send on output slot 1
  TList    *fOutputMassPt;        //!<! list send on output slot 6
  TList    *fOutputMassY;        //!<! list send on output slot 9
  TList    *fDistr;               //!<! list send on output slot 2
  TH1F     *fNentries;            //!<! histogram with number of events on output slot 3
  THnSparseF *fMCAccPrompt;       //!<!histo for StepMCAcc for D0 prompt (pt,y,ptB)
  THnSparseF *fMCAccBFeed;        //!<!histo for StepMCAcc for D0 FD (pt,y,ptB)
  Bool_t fStepMCAcc;              // flag to activate histos for StepMCAcc
  AliRDHFCutsD0toKpi *fCuts;      //  Cuts - sent to output slot 4
  Bool_t    fEnableCentralityCorrCuts; /// flag to enable centrality correlation event cuts
  AliEventCuts  fEventCuts;       // Event cut object for centrality correlation event cuts
  THnSparseF *fHistMassPtImpParTC[5];   //!<! histograms for impact paramter studies
  Int_t     fArray;               ///  can be D0 or Like Sign candidates
  Bool_t    fReadMC;              ///  flag for MC array: kTRUE = read it, kFALSE = do not read it
  Bool_t    fCutOnDistr;          ///  flag to decide if apply cut also on distributions: 0 no cuts, 1 looser cuts, 2 tighter/ cuts
  Bool_t    fUsePid4Distr;        //  flag to use the particle identification to fill the signal histograms of distributions. It has effect only with fReadMC=kFALSE
  AliNormalizationCounter *fCounter;//!<! AliNormalizationCounter on output slot 5
  Int_t     fNPtBins;             ///  number of pt bins
  Double_t  fLsNormalization;     ///  normalization
  Int_t     fFillOnlyD0D0bar;     /// flag to fill mass histogram with D0/D0bar only (0 = fill with both, 1 = fill with D0 only, 2 = fill with D0bar only)
  TObjArray fDaughterTracks;      /// keeps the daughter tracks
  Int_t     fIsSelectedCandidate; /// selection outcome
  Bool_t    fFillVarHists;        /// flag to enable filling variable histos
  Int_t     fSys;                 /// fSys=0 -> p-p; fSys=1 ->PbPb (in this case fFillVarHists=kFALSE by default: set it to kTRUE *after* if needed)
  Bool_t    fIsRejectSDDClusters; /// flag to reject events with SDD clusters
  Bool_t    fFillPtHist;          /// flag to fill Pt and Impact Parameter Histograms
  Bool_t    fFillYHist;          /// flag to fill Y Histograms
  Bool_t    fFillImpParHist;      /// flag to fill Pt and Impact Parameter Histograms
  Bool_t    fFillSubSampleHist;    /// flag to fill SubSample histogram
  Int_t     fEventCounter; /// event counter used for sub sample test
  Bool_t    fUseSelectionBit;     /// flag to check or not the selection bit
  Int_t     fAODProtection;       /// flag to activate protection against AOD-dAOD mismatch.
                                  /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names

  Bool_t    fWriteVariableTree;       /// flag to decide whether to write the candidate variables on a tree variables
  TTree    *fVariablesTree;           //!<! tree of the candidate variables after track selection on output slot 7
  Double_t *fCandidateVariables;      //!<!  variables to be written to the tree
  Bool_t    fWriteProtosgnVar;        /// flag to decide whether to write the selected candidates variables on a tree for cut optimization
  Bool_t    fSelectTrueD0;            /// flag to decide whether to write only true D0/D0bar
  Bool_t    fUsedMassWindow;          /// flag to activate the mass window selection for output size reduction
  Bool_t	  fPIDCheck;			/// flag to decide whether to fill "PID = x" bins in fNentrie
  Bool_t    fDrawDetSignal;		/// flag to decide whether to draw the TPC dE/dx and TOF signal before/after PID
  Bool_t    fUseQuarkTagInKine;            // flag for quark/hadron level identification of prompt and feeddown
  Bool_t    fFillSparses;                  // flag to activate THnSparse
  Bool_t    fUseRejectionMethod;           // flag to activate the Rejection method
  Float_t   fRejectionFactor;              // rejection factor to be used in the rejection method
  THnSparseF *fhStudyImpParSingleTrackSign; //!<! sparse with imp par residual cuts for MC
  THnSparseF *fhStudyImpParSingleTrackCand;  //!<! sparse with imp par residual cuts for Data
  THnSparseF *fhStudyImpParSingleTrackFd;   //!<! sparse with imp par residual cuts for MC
  TList	   *fDetSignal;		//!<!Detector signal histograms (on output slot 8)
  TH2F *fhMultVZEROTPCoutTrackCorrNoCut;  //!<!
  TH2F *fhMultVZEROTPCoutTrackCorr;  //!<!
  Bool_t    fEnablePileupRejVZEROTPCout;
  TH2F *fhMultVZEROTPCclustersCorrNoCut;  //!<!
  TH2F *fhMultVZEROTPCclustersCorr;  //!<!
  Bool_t    fEnablePileupRejVZEROTPCcls;
  Bool_t    fRejectOutOfBunchPileUp;
  TNtuple *fNtupleD0C;//!<!
  TNtuple *fNtupleD0B;//!<!
  TNtuple *fNtupleD0Data;//!<!
  TNtuple *fNtupleRefl;//!<!



  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSED0MassTMVA,24); /// AliAnalysisTaskSE for D0->Kpi
  /// \endcond
};

#endif
