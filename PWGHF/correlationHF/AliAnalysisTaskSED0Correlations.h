#ifndef ALIANALYSISTASKSED0CORRELATIONS_H
#define ALIANALYSISTASKSED0CORRELATIONS_H

/* Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

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

  void SetReadMC(Bool_t readMC=kFALSE){fReadMC=readMC;}
  void SetMCReconstructedTracks(Bool_t recoTrMC=kTRUE){fRecoTr=recoTrMC;}
  void SetMCReconstructedD0(Bool_t recoD0MC=kTRUE){fRecoD0=recoD0MC;}
  void SetMCSelEventType(Bool_t sel=kFALSE){fSelEvType=sel;}
  void SetFillOnlyD0D0bar(Int_t flagfill){fFillOnlyD0D0bar=flagfill;}
  void SetSystem(Int_t sys){fSys=sys;}
  void SetRejectSDDClusters(Bool_t flag) {fIsRejectSDDClusters=flag; }
  void SetFillGlobalPlots(Bool_t fill=kTRUE){fFillGlobal=fill;}

  Int_t  GetReadMC() const {return fReadMC;}
  Int_t  GetMCReconstructedTracks() const {return fRecoTr;}
  Int_t  GetMCReconstructedD0() const {return fRecoD0;}
  Int_t  GetMCSelEventType() const {return fSelEvType;}
  Int_t  GetFillOnlyD0D0bar() const {return fFillOnlyD0D0bar;}
  Int_t  GetSystem() const {return fSys;}
  Bool_t GetRejectSDDClusters() const {return fIsRejectSDDClusters;}
  Bool_t GetFillGlobalPlots() const {return fFillGlobal;}
  Double_t GetEtaForCorrel() {return fEtaForCorrel;}

  //correlations setters/printers
  void SetNPtBinsCorr(Int_t nbins) {fNPtBinsCorr = nbins;}
  void SetPtBinsLimsCorr(Double_t* ptlims) {for(int i=0;i<=fNPtBinsCorr;i++) {fBinLimsCorr.push_back(ptlims[i]);}}
  void SetPtBinsLimsCorr(Float_t* ptlims) {for(int i=0;i<=fNPtBinsCorr;i++) {fBinLimsCorr.push_back((Double_t)ptlims[i]);}}
  void SetPtTreshLow(Double_t* pttreshlow) {for(int i=0;i<fNPtBinsCorr;i++) {fPtThreshLow.push_back(pttreshlow[i]);}}
  void SetPtTreshUp(Double_t* pttreshup) {for(int i=0;i<fNPtBinsCorr;i++) {fPtThreshUp.push_back(pttreshup[i]);}}
  void PrintBinsAndLimits();
  Int_t PtBinCorr(Double_t pt) const;
  void SetEvMixing(Bool_t mix) {fMixing=mix;}
  void SetEtaForCorrel(Double_t etacorr) {fEtaForCorrel=etacorr;}

  enum PartType {kTrack,kKCharg,kK0};
  enum FillType {kSE, kME}; //for single event or event mixing histos fill

 private:

  AliAnalysisTaskSED0Correlations(const AliAnalysisTaskSED0Correlations &source);
  AliAnalysisTaskSED0Correlations& operator=(const AliAnalysisTaskSED0Correlations& source); 
  void FillMassHists(AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi *cuts, TList *listout);
  Int_t CheckD0Origin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const;
  //correlation methods
  void CreateCorrelationsObjs();
  void CalculateCorrelations(AliAODRecoDecayHF2Prong* d, Int_t labD0=-1, TClonesArray* mcArray=0x0);
  void CalculateCorrelationsMCKine(AliAODMCParticle* d, TClonesArray* mcArray=0x0);
  void FillSparsePlots(TClonesArray* arrayMC, Double_t mInv[], Int_t origD0, Int_t PdgD0, AliReducedParticle* track, Int_t ptbin, Int_t type, Double_t wg=1.);
  Int_t CheckTrackOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPartCandidate) const;
  Bool_t IsDDaughter(AliAODMCParticle* d, AliAODMCParticle* track, TClonesArray* mcArray) const;
  Bool_t SelectV0(AliAODv0* v0, AliAODVertex *vtx, Int_t option, Int_t idArrayV0[][2]) const;

  Int_t             fNPtBinsCorr;        // number of pt bins per correlations
  std::vector<Double_t>  fBinLimsCorr;        // limits of pt bins per correlations
  std::vector<Double_t>  fPtThreshLow;        // pT treshold of hadrons - low
  std::vector<Double_t>  fPtThreshUp;         // pT treshold of hadrons - up

  Int_t     fEvents;		  	// EventCounter
  Bool_t    fAlreadyFilled;	  	// D0 in an event already analyzed (for track distribution plots)
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

  ClassDef(AliAnalysisTaskSED0Correlations,2); // AliAnalysisTaskSE for D0->Kpi
};

#endif

