#ifndef AliAnalysisTaskSED0Correlations_H
#define AliAnalysisTaskSED0Correlations_H

/* Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskSED0Correlations.h 63031 2013-06-17 13:30:19Z arossi $ */

//*************************************************************************
// Class AliAnalysisTaskSED0Correlations
// AliAnalysisTaskSE for D0 candidates (2Prongs) and hadrons correlations
// Authors:
// C.Bianchin, chiara.bianchin@pd.infn.it
// F.Colamaria, fabio.colamaria@ba.infn.it
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1F.h>
#include <THnSparse.h>

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliHFCorrelator.h"
#include "AliNormalizationCounter.h"
#include "AliHFOfflineCorrelator.h"
#include "AliD0hCutOptim.h"

using std::vector;

class AliAODEvent;

class AliAnalysisTaskSED0Correlations : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSED0Correlations();
  AliAnalysisTaskSED0Correlations(const char *name,AliRDHFCutsD0toKpi* cutsD0, AliHFAssociatedTrackCuts* cutsTrk);
  virtual ~AliAnalysisTaskSED0Correlations();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  enum PartType {kTrack,kKCharg,kK0};
  enum FillType {kSE, kME}; //for single event or event mixing histos fill
  enum TreeFill {kNoTrees, kFillTrees, kFillCutOptTree};

  void SetReadMC(Bool_t readMC=kFALSE){fReadMC=readMC;}
  void SetMCReconstructedTracks(Bool_t recoTrMC=kTRUE){fRecoTr=recoTrMC;}
  void SetMCReconstructedD0(Bool_t recoD0MC=kTRUE){fRecoD0=recoD0MC;}
  void SetMCSelEventType(Bool_t sel=kFALSE){fSelEvType=sel;}
  void SetFillOnlyD0D0bar(Int_t flagfill){fFillOnlyD0D0bar=flagfill;}
  void SetSystem(Int_t sys){fSys=sys;}
  void SetRejectSDDClusters(Bool_t flag) {fIsRejectSDDClusters=flag; }
  void SetFillGlobalPlots(Bool_t fill=kTRUE){fFillGlobal=fill;}
  void SetSoftPiFlag(Bool_t piflag) {fSoftPiCut=piflag;}
  void SetMEAxisThresh(Bool_t methresh) {fMEAxisThresh=methresh;}
  void SetKaonCorrelations(Bool_t kaonCorr) {fKaonCorr=kaonCorr;}
  void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}
  void SetPurityStudies(Bool_t puritystudies=kFALSE) {fPurityStudies=puritystudies;}

  Int_t  GetReadMC() const {return fReadMC;}
  Int_t  GetMCReconstructedTracks() const {return fRecoTr;}
  Int_t  GetMCReconstructedD0() const {return fRecoD0;}
  Int_t  GetMCSelEventType() const {return fSelEvType;}
  Int_t  GetFillOnlyD0D0bar() const {return fFillOnlyD0D0bar;}
  Int_t  GetSystem() const {return fSys;}
  Bool_t GetRejectSDDClusters() const {return fIsRejectSDDClusters;}
  Bool_t GetFillGlobalPlots() const {return fFillGlobal;}
  Double_t GetEtaForCorrel() {return fEtaForCorrel;}
  Double_t GetMultEv() {return fMultEv;}
  Bool_t GetSoftPiFlag() const {return fSoftPiCut;}
  Bool_t GetMEAxisThresh() const {return fMEAxisThresh;}
  Bool_t GetKaonCorrelations() const {return fKaonCorr;}
  Bool_t GetFillTrees() const {return fFillTrees;}

  //correlations setters/printers
  void SetNPtBinsCorr(Int_t nbins) {fNPtBinsCorr = nbins;}
  void SetPtBinsLimsCorr(Double_t* ptlims) {for(int i=0;i<=fNPtBinsCorr;i++) {fBinLimsCorr.push_back(ptlims[i]);}}
  void SetPtBinsLimsCorr(Float_t* ptlims) {for(int i=0;i<=fNPtBinsCorr;i++) {fBinLimsCorr.push_back((Double_t)ptlims[i]);}}
  void SetPtTreshLow(Double_t* pttreshlow) {for(int i=0;i<fNPtBinsCorr;i++) {fPtThreshLow.push_back(pttreshlow[i]);}}
  void SetPtTreshUp(Double_t* pttreshup) {for(int i=0;i<fNPtBinsCorr;i++) {fPtThreshUp.push_back(pttreshup[i]);}}
  void SetLSBLowLim(Double_t* LSBLowLim) {for(int i=0;i<fNPtBinsCorr;i++) {fLSBLowLim.push_back(LSBLowLim[i]);}}
  void SetLSBHighLim(Double_t* LSBUppLim) {for(int i=0;i<fNPtBinsCorr;i++) {fLSBUppLim.push_back(LSBUppLim[i]);}}  
  void SetRSBLowLim(Double_t* RSBLowLim) {for(int i=0;i<fNPtBinsCorr;i++) {fRSBLowLim.push_back(RSBLowLim[i]);}}
  void SetRSBHighLim(Double_t* RSBUppLim) {for(int i=0;i<fNPtBinsCorr;i++) {fRSBUppLim.push_back(RSBUppLim[i]);}}
  void SetLeftSignReg_LowPt(Double_t leftlow) {fSignLeft_LowPt=leftlow;}
  void SetRightSignReg_LowPt(Double_t rightlow) {fSignRight_LowPt=rightlow;}
  void SetLeftSignReg_HighPt(Double_t lefthigh) {fSignLeft_HighPt=lefthigh;}
  void SetRightSignReg_HighPt(Double_t righthigh) {fSignRight_HighPt=righthigh;}
  void SetPtAssocLim(Double_t pTlim) {fPtAssocLimit=pTlim;}
  
  void PrintBinsAndLimits();
  Int_t PtBinCorr(Double_t pt) const;
  void SetEvMixing(Bool_t mix) {fMixing=mix;}
  void SetEtaForCorrel(Double_t etacorr) {fEtaForCorrel=etacorr;}
  void SetSpeed(Bool_t speed) {fSpeed=speed;}
  void SetMergePools(Bool_t mergepools) {fMergePools=mergepools;}
  void SetUseDeff(Bool_t useDeff) {fUseDeff=useDeff;}
  void SetUseTrackeff(Bool_t useTrackeff) {fUseTrackeff=useTrackeff;}
  void SetMinDPt(Double_t minDPt) {fMinDPt=minDPt;}
  void SetFillTrees(TreeFill fillTrees, Double_t fractAccME) {fFillTrees=fillTrees; fFractAccME=fractAccME;}
 
 private:

  AliAnalysisTaskSED0Correlations(const AliAnalysisTaskSED0Correlations &source);
  AliAnalysisTaskSED0Correlations& operator=(const AliAnalysisTaskSED0Correlations& source); 
  void FillMassHists(AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi *cuts, TList *listout, AliAODEvent *aod);
  Int_t CheckD0Origin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const;
  //correlation methods
  void CreateCorrelationsObjs();
  void CalculateCorrelations(AliAODRecoDecayHF2Prong* d, Int_t labD0=-1, TClonesArray* mcArray=0x0);
  void CalculateCorrelationsMCKine(AliAODMCParticle* d, TClonesArray* mcArray=0x0);
  void FillSparsePlots(TClonesArray* arrayMC, Double_t mInv[], Int_t origD0, Int_t PdgD0, AliReducedParticle* track, Int_t ptbin, Int_t type, Int_t softpiME, Double_t wg=1.);
  Int_t CheckTrackOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const;
  Bool_t IsDDaughter(AliAODMCParticle* d, AliAODMCParticle* track, TClonesArray* mcArray) const;
  Bool_t SelectV0(AliAODv0* v0, AliAODVertex *vtx, Int_t option, Int_t idArrayV0[][2]) const;
  Bool_t IsSoftPion_MCKine(AliAODMCParticle* d, AliAODMCParticle* track, TClonesArray* arrayMC) const;
  void FillTreeD0(AliAODRecoDecayHF2Prong* d, AliAODEvent* aod);  
  void FillTreeTracks(AliAODEvent* aod);  
  void FillTreeD0ForCutOptim(AliAODRecoDecayHF2Prong* d, AliAODEvent* aod);  
  void ResetBranchD();
  void ResetBranchTracks();
  void ResetBranchDForCutOptim();
  Bool_t AcceptTrackForMEOffline(Double_t pt);
  void FillPurityPlots(TClonesArray* mcArray, AliReducedParticle* track, Int_t ptbin, Double_t deltaphi);
  
  Int_t             	 fNPtBinsCorr;        // number of pt bins per correlations
  std::vector<Double_t>  fBinLimsCorr;        // limits of pt bins per correlations
  std::vector<Double_t>  fPtThreshLow;        // pT threshold of hadrons - low
  std::vector<Double_t>  fPtThreshUp;         // pT threshold of hadrons - up
  std::vector<Double_t>  fLSBLowLim;          // Left SB lower lim
  std::vector<Double_t>  fLSBUppLim;          // Left SB lower lim
  std::vector<Double_t>  fRSBLowLim;          // Right SB upper lim
  std::vector<Double_t>  fRSBUppLim;          // Right SB upper lim
  std::vector<Int_t>     fDaughTrackID;       // ID of tagged daughters
  std::vector<Int_t>     fDaughTrigNum;	      // ID of D-trigger for daughters	

  Int_t     fEvents;		  	// EventCounter
  Bool_t    fAlreadyFilled;	  	// D0 in an event already analyzed (for track distribution plots)
  Int_t	    fNtrigD;			// counter on number of D triggers filled (for association with daughter tracks in TTrees)
  TList    *fOutputMass;          	//!list send on output slot 1
  TList    *fOutputCorr;	  	//!list of correlation histos, output slot 5
  TList    *fOutputStudy;	  	//!list of histos with MC distributions, output slot 6
  TH1F     *fNentries;            	//!histogram with number of events on output slot 2
  AliRDHFCutsD0toKpi *fCutsD0;    	// Cuts for D0, output 3
  AliHFAssociatedTrackCuts *fCutsTracks;// Cuts for tracks and K0, output 7
  AliHFCorrelator* fCorrelatorTr;	// Correlator for tracks
  AliHFCorrelator* fCorrelatorKc;	// Correlator for charged K
  AliHFCorrelator* fCorrelatorK0;	// Correlator for K0
  Bool_t    fReadMC;              	// flag for MC array: kTRUE = read it, kFALSE = do not read it
  Bool_t    fRecoTr;   		       	// flag for using MC reconstructed (kTRUE) or pure kinematic MC (kFALSE) - Associated tracks
  Bool_t    fRecoD0;   		       	// flag for using MC reconstructed (kTRUE) or pure kinematic MC (kFALSE) - D0
  Bool_t    fSelEvType;		       	// flag for enabling selection of event tpye (PP, GS, FE, ...) on MC analysis
  Bool_t    fMixing;			// flag to enable also event mixing
  AliNormalizationCounter *fCounter;	//!AliNormalizationCounter on output slot 4
  Int_t     fNPtBins;             	// Number of pt bins
  Int_t     fFillOnlyD0D0bar;     	// flag to fill mass histogram with D0/D0bar only (0 = fill with both, 1 = fill with D0 only, 2 = fill with D0bar only)
  Int_t     fIsSelectedCandidate; 	// selection outcome
  Int_t     fSys;                 	// fSys=0 -> p-p; fSys=1 ->PbPb
  Double_t  fEtaForCorrel;		// cut for D0 eta to enable correlation with associated particles
  Bool_t    fIsRejectSDDClusters; 	// flag to reject events with SDD clusters
  Bool_t    fFillGlobal;          	// flag to fill global plots (in loops on tracks and V0 for each event)
  Double_t  fMultEv;			// event multiplicity (for trigger eff)
  Double_t  fzVtx;				// event zVtx position (for track eff)
  Bool_t    fSoftPiCut;			// flag to activate soft pion cut on Data
  Bool_t    fMEAxisThresh;		// flag to fill threshold axis in ME plots
  Bool_t    fKaonCorr;			// enables correlations of D0-Kcharg and D0-K0
  Double_t  fSignLeft_LowPt;		// Left bound of "signal region" range - up to 8 GeV/c
  Double_t  fSignRight_LowPt;		// Right bound of "signal region" range - up to 8 GeV/c
  Double_t  fSignLeft_HighPt;		// Left bound of "signal region" range - from 8 GeV/c
  Double_t  fSignRight_HighPt;		// Right bound of "signal region" range - from 8 GeV/c
  Int_t     fPoolNum;			// Number of the pool for the analyzed event
  Bool_t    fSpeed;			// Speed up the execution removing bins and histos
  Bool_t    fMergePools;		// Put all entries from various pools in _pool0 THnSparses (as old approach) - for testing & low stat!
  Bool_t    fUseDeff;			// Use D meson efficiency as weight
  Bool_t    fUseTrackeff;   		// Use track efficiency as weight
  Double_t  fPtAssocLimit;   		// Maximum value for associated pT
  Double_t  fMinDPt;			// Minimum pT of the D0 to allow selection

  TreeFill  fFillTrees;			// Flag to fill ME offline trees
  Double_t  fFractAccME;		// Fraction of tracks to be accepted in the ME offline
  Int_t     fAODProtection;  	        // flag to activate protection against AOD-dAOD mismatch.
  Bool_t    fPurityStudies;		// flag to activate purity studies (primaries, secondaries, charm and beauth tracks rejected by DCA cut, vs pT and deltaPhi)

  AliHFCorrelationBranchD   *fBranchD;
  AliHFCorrelationBranchTr  *fBranchTr;
  AliD0hCutOptim	    *fBranchDCutVars; //for cut optimization!

  TTree	    *fTreeD;			// TTree for ME offline - D0 mesons
  TTree	    *fTreeTr;			// TTree for ME offline - Assoc tracks
  TObjArray *fTrackArray;		// Array with selected tracks for association
  Bool_t    fTrackArrayFilled;		// Flag to fill fTrackArray or not (if already filled)

  ClassDef(AliAnalysisTaskSED0Correlations,14); // AliAnalysisTaskSE for D0->Kpi - h correlations
};

#endif

