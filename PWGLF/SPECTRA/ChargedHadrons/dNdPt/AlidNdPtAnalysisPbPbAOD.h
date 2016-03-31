#ifndef AlidNdPtAnalysisPbPbAOD_H
#define AlidNdPtAnalysisPbPbAOD_H


//------------------------------------------------------------------------------
// AlidNdPtAnalysisPbPbAOD class used for dNdPt analysis in PbPb collision
// via AODs
//
// Author: P. Luettig, 15.05.2013
// last modified: 10.06.2014
//------------------------------------------------------------------------------



class iostream;

#include "AliAnalysisTaskSE.h"
#include "TObject.h"
#include "TList.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "THn.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TObjArray.h"
#include "TObjString.h"

#include "TParticlePDG.h"
#include "TDatabasePDG.h"

#include "AliLog.h"
#include "AliCentrality.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"

#include "AliInputEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliExternalTrackParam.h"
#include "AliESDtrack.h"
#include "AliEventplane.h"

#include "AliEPSelectionTask.h"

#include "TSystem.h"
#include "TROOT.h"

class AlidNdPtAnalysisPbPbAOD : public AliAnalysisTaskSE {
  public :
  enum CheckQuantity { cqCrossedRows = 0, cqNcluster = 1, cqChi = 2, cqLength = 3, cqRowsOverFindable = 4 };
  enum KinematicQuantity { kqPt = 0, kqEta = 1, kqPhi = 2 };
  enum MaxCheckQuantity { cqMax = 5 };
  enum MaxKinematicQuantity { kqMax = 3 };
  
  AlidNdPtAnalysisPbPbAOD(const char *name = "dNdPtPbPbAOD");
  ~AlidNdPtAnalysisPbPbAOD();
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  
  // Set binning for Histograms (if not set default binning is used)
  void SetBinsMult(Int_t nbins, Double_t* edges)      { Printf("[Set] MultBins"); fMultNbins = nbins; fBinsMult = GetArrayClone(nbins,edges); }
  void SetBinsMultFine(Int_t nbins, Double_t* edges) 	{ Printf("[Set] FineMultBins"); fMultFineNbins = nbins; fBinsMultFine = GetArrayClone(nbins,edges); }
  void SetBinsPt(Int_t nbins, Double_t* edges)        { Printf("[Set] pTBins"); fPtNbins = nbins; fBinsPt = GetArrayClone(nbins,edges); }
  void SetBinsPtCorr(Int_t nbins, Double_t* edges)    { Printf("[Set] pTcorrBins"); fPtCorrNbins = nbins; fBinsPtCorr = GetArrayClone(nbins,edges); }
  void SetBinsPtCheck(Int_t nbins, Double_t* edges) 	{ Printf("[Set] pTcheckBins"); fPtCheckNbins = nbins; fBinsPtCheck = GetArrayClone(nbins,edges); }
  void SetBinsEta(Int_t nbins, Double_t* edges)       { Printf("[Set] EtaBins"); fEtaNbins = nbins; fBinsEta = GetArrayClone(nbins,edges); }
  void SetBinsEtaCheck(Int_t nbins, Double_t* edges) 	{ Printf("[Set] EtaCheckBins"); fEtaCheckNbins = nbins; fBinsEtaCheck = GetArrayClone(nbins,edges); }
  void SetBinsZv(Int_t nbins, Double_t* edges)        { Printf("[Set] ZvBins"); fZvNbins = nbins; fBinsZv= GetArrayClone(nbins,edges); }
  void SetBinsCentrality(Int_t nbins, Double_t* edges){ Printf("[Set] CentBins"); fCentralityNbins = nbins; fBinsCentrality = GetArrayClone(nbins,edges); }
  void SetBinsPhi(Int_t nbins, Double_t* edges)       { Printf("[Set] PhiBins"); fPhiNbins = nbins; fBinsPhi = GetArrayClone(nbins,edges); }
  void SetBinsPhiCorr(Int_t nbins, Double_t* edges)   { Printf("[Set] PhiCorrBins"); fPhiCorrNbins = nbins; fBinsPhiCorr = GetArrayClone(nbins,edges); }
  void SetBinsDeltaphi(Int_t nbins, Double_t* edges) 	{ Printf("[Set] DeltaphiBins"); fDeltaphiNbins = nbins; fBinsDeltaphi = GetArrayClone(nbins,edges); }
  void SetBinsRunNumber(Int_t nbins, Double_t* edges) { Printf("[Set] RunNumberBins"); fRunNumberNbins = nbins; fBinsRunNumber = GetArrayClone(nbins,edges); }
  
  // set event cut variables
  void SetCutMaxZVertex( Double_t d)					    { fCutMaxZVertex = d; }
  Double_t GetCutMaxZVertex()						        { return fCutMaxZVertex; }
  
  // set Ncontributors to vertex
  void SetNContributorsVertex(Int_t i)					{ fVertexMinContributors = i; }
  Int_t GetNContributorsVertex()							{ return fVertexMinContributors; }
  
  // set track kinematic cut parameters
  void SetCutPtRange(Double_t ptmin, Double_t ptmax)		{ fCutPtMin = ptmin; fCutPtMax = ptmax; }
  Double_t GetCutPtMin()						            { return fCutPtMin; }
  Double_t GetCutPtMax()						            { return fCutPtMax; }
  
  void SetCutEtaRange(Double_t etamin, Double_t etamax)	{ fCutEtaMin = etamin; fCutEtaMax = etamax; }
  Double_t GetCutEtaMin()						            { return fCutEtaMin; }
  Double_t GetCutEtaMax()						            { return fCutEtaMax; }
  
  void EnableRelativeCuts()								{ Printf("[I] Relative Cuts enabled"); fUseRelativeCuts = kTRUE; }
  Bool_t AreRelativeCutsEnabled()							{ return fUseRelativeCuts; }
  
  void FillEventPtSpectraHistogram(Bool_t b)				{ fFillEventPtSpectraHistogram = b; }
  Bool_t GetFillEventPtSpectraHistogram()					{ return fFillEventPtSpectraHistogram; }
  
  // setter and getter track quality cut parameters
  void SetFilterBit(Int_t b)								{ fFilterBit = b; };
  Int_t GetFilterBit()									{ return fFilterBit; }
  
  void SetRequireHybridTracking(Bool_t b)					{ fHybridTracking = b; }
  Bool_t RequireHybridTracking()							{ return fHybridTracking; }
  
  void SetCutRequireTPCRefit(Bool_t *b) 					{ fCutRequireTPCRefit = b; }
  Bool_t IsTPCRefitRequired() 							{ return fCutRequireTPCRefit; }
  
  void SetCutRequireITSRefit(Bool_t *b) 					{ fCutRequireITSRefit = b; }
  Bool_t IsITSRefitRequired() 							{ return fCutRequireITSRefit; }
  
  void SetCutMinNClustersTPC(Double_t d)					{ fCutMinNumberOfClusters = d; }
  Double_t GetCutMinNClustersTPC()						{ return fCutMinNumberOfClusters; }
  
  void SetCutPercMinNClustersTPC(Double_t d)				{ Printf("[I] Take only %.2f%% tracks with most clusters", d*100.); fCutPercMinNumberOfClusters = d; }
  Double_t GetCutPercMinNClustersTPC()					{ return fCutPercMinNumberOfClusters; }
  
  void SetCutMinNCrossedRowsTPC(Double_t d) 				{ fCutMinNumberOfCrossedRows = d; }
  Double_t GetCutMinNCrossedRowsTPC()						{ return fCutMinNumberOfCrossedRows; }
  
  void SetCutPercMinNCrossedRowsTPC(Double_t d) 			{ Printf("[I] Take only %.2f%% tracks with most crossedRows", d*100.); fCutPercMinNumberOfCrossedRows = d; }
  Double_t GetCutPercMinNCrossedRowsTPC()					{ return fCutPercMinNumberOfCrossedRows; }
  
  void SetCutMinRatioCrossedRowsOverFindableClustersTPC(Double_t d) 	{ fCutMinRatioCrossedRowsOverFindableClustersTPC = d; }
  Double_t GetCutMinRatioCrossedRowsOverFindableClustersTPC()			{ return fCutMinRatioCrossedRowsOverFindableClustersTPC; }
  
  void SetCutLengthInTPCPtDependent(Bool_t b)				{ fCutLengthInTPCPtDependent = b; }
  Bool_t DoCutLengthInTPCPtDependent()					{ return fCutLengthInTPCPtDependent; }
  
  void SetPrefactorLengthInTPCPtDependent(Double_t d)		{ fPrefactorLengthInTPCPtDependent = d; }
  Double_t GetPrefactorLengthInTPCPtDependent()			{ return fPrefactorLengthInTPCPtDependent; }
  
  void SetCutMaxChi2PerClusterTPC(Double_t d) 			{ fCutMaxChi2PerClusterTPC = d; }
  void SetCutMaxFractionSharedTPCClusters(Double_t d) 	{ fCutMaxFractionSharedTPCClusters = d; }
  void SetCutMaxDCAToVertexZ(Double_t d) 					{ fCutMaxDCAToVertexZ = d; }
  void SetCutMaxDCAToVertexXY(Double_t d) 				{ fCutMaxDCAToVertexXY = d; }
  void SetCutMaxChi2PerClusterITS(Double_t d) 			{ fCutMaxChi2PerClusterITS = d; }
  void SetCutDCAToVertex2D(Bool_t *b) 					{ fCutDCAToVertex2D = b; }
  void SetCutRequireSigmaToVertex(Bool_t *b) 				{ fCutRequireSigmaToVertex = b; }
  void SetCutMaxDCAToVertexXYPtDep(Double_t d0, Double_t d1, Double_t d2)
  {
    fCutMaxDCAToVertexXYPtDepPar0 = d0;
    fCutMaxDCAToVertexXYPtDepPar1 = d1;
    fCutMaxDCAToVertexXYPtDepPar2 = d2;
  }
  void SetCutAcceptKinkDaughters(Bool_t *b) 				{ fCutAcceptKinkDaughters = b; }
  void SetCutMaxChi2TPCConstrainedGlobal(Double_t d) 		{ fCutMaxChi2TPCConstrainedGlobal = d; }
  
  // fill function for cross check histos
  Bool_t FillDebugHisto(Double_t *dCrossCheckVar, Double_t *dKineVar, Double_t dCentrality, Bool_t bIsAccepted);
  
  // fill function for cut settings
  void StoreCutSettingsToHistogram();
  
  // getter for DCA
  Bool_t GetDCA(const AliAODTrack *track, AliAODEvent *evt, Double_t d0z0[2]);
  
  THnSparseF * GetHistZvPtEtaCent() const { return fZvPtEtaCent; }
  TH1F * GetHistEventStatistics() const { return fEventStatistics; }
  
  const char * GetParticleName(Int_t pdg);
  
  AliGenHijingEventHeader* GetHijingEventHeader(AliAODMCHeader *header);
  AliGenPythiaEventHeader* GetPythiaEventHeader(AliAODMCHeader *header);
  
  Double_t RotatePhi(Double_t phiTrack, Double_t phiEP, Double_t dMaxDeltaPhi);
  // 	Double_t MoveEventplane(Double_t dMCEP);
  
  Bool_t SetRelativeCuts(AliAODEvent *event);
  
  Bool_t IsTrackAccepted(AliAODTrack *tr, Double_t dCentrality, Double_t bMagZ);
  Bool_t IsMCTrackAccepted(AliAODMCParticle *part);
  
  Bool_t IsHijingParticle(const AliAODMCParticle *part, AliGenHijingEventHeader* hijingGenHeader);
  Bool_t IsPythiaParticle(const AliAODMCParticle *part, AliGenPythiaEventHeader* pythiaGenHeader);
  
  static Double_t* GetArrayClone(Int_t n, Double_t* source);
  
  void SetEventplaneSelector(char *c) { fEPselector = c; }
  TString GetEventplaneSelector() { return fEPselector; }
  
  void SetCentralityEstimator(char *c) { fCentEstimator = c; }
  TString GetCentralityEstimator() { return fCentEstimator; }
  
  void SetDoMinBiasAnalysis(Bool_t b) { fDoMinBiasAnalysis = b; }
  Bool_t GetDoMinBiasAnalysis() { return fDoMinBiasAnalysis; }
  
  void EnableCrossCheckCorrelationHistos() { fCrossCheckCorrelHisto = kTRUE; }
  Bool_t AreCrossCheckCorrelationHistosEnabled() { return fCrossCheckCorrelHisto; }
  
  void DisableOnlineTriggerStrings(char *c) { fDisabledTriggerString = c; }
  TString GetDisabledOnlineTrigger() { return fDisabledTriggerString; }
  
  void SetAnchorPointSystematicFactor(Double_t d) { fAnchorPointCorrectionFactor = d; fDoAnchorPointSystStudy = kTRUE; }
  Bool_t DoAnchorPointSystStudy() { return fDoAnchorPointSystStudy; }
  Double_t GetAnchorPointCorrectionFactor() { return fAnchorPointCorrectionFactor; }
  
  private :
  
  // Output List
  TList		*fOutputList;
  
  // Histograms
  TH1F        *fPt; // simple pT histogramm
  TH1F        *fMCPt; // simple pT truth histogramm
  THnSparseF 	*fZvPtEtaCent; //-> Zv:Pt:Eta:Cent
  THnSparseF 	*fDeltaphiPtEtaPhiCent; //-> DeltaPhi:Pt:Eta:Phi:Cent, was fDeltaphiPtEtaCent
  THnSparseF 	*fDeltaphiPtEtaPhiZvCent; //-> DeltaPhi:Pt:Eta:Phi:Zv:Cent, was fDeltaphiPtEtaCent
  THnSparseF 	*fPtResptCent; //-> 1/pt:ResolutionPt:Cent
  THnSparseF 	*fPtResptptCent; //-> 1/pt:ResolutionPt*Pt:Cent
  TH2F        *fPtEvent; // pT per event, for 200 events
  
  THnSparseF 	*fMCRecPrimZvPtEtaCent; //-> MC Zv:Pt:Eta:Cent
  THnSparseF 	*fMCGenZvPtEtaCent; //-> MC Zv:Pt:Eta:Cent
  THnSparseF 	*fMCRecSecZvPtEtaCent; //-> MC Zv:Pt:Eta:Cent, only secondaries
  
  THnF        *fMCPtEtaPhiCent; //-> MC Pt:Eta:Phi:Cent
  THnF        *fMCRecPrimPtEtaPhiCent; //-> MC Pt:Eta:Phi:Cent, was fMCRecPrimDeltaphiPtEtaCent
  THnF        *fMCGenPtEtaPhiCent; //-> MC Pt:Eta:Phi:Cent, was fMCGenDeltaphiPtEtaCent
  THnF        *fMCRecSecPtEtaPhiCent; //-> MC Pt:Eta:Phi:Cent, only secondaries, was fMCRecSecDeltaphiPtEtaCent
  
  THnF        *fMCPtEtaPhiZvCent; //-> MC Pt:Eta:Phi:Zv:Cent
  THnF        *fMCRecPrimPtEtaPhiZvCent; //-> MC Pt:Eta:Phi:Zv:Cent, was fMCRecPrimDeltaphiPtEtaCent
  THnF        *fMCGenPtEtaPhiZvCent; //-> MC Pt:Eta:Phi:Zv:Cent, was fMCGenDeltaphiPtEtaCent
  THnF        *fMCRecSecPtEtaPhiZvCent; //-> MC Pt:Eta:Phi:Zv:Cent, only secondaries, was fMCRecSecDeltaphiPtEtaCent
  
  TH1F        *fEventStatistics; // contains statistics of number of events after each cut
  TH1F        *fEventStatisticsCentrality; // contains number of events vs centrality, events need to have a track in kinematic range
  TH1F        *fMCEventStatisticsCentrality; // contains MC number of events vs centrality, events need to have a track in kinematic range
  TH1F        *fAllEventStatisticsCentrality; // contains number of events vs centrality, events need to be triggered
  TH2F        *fEventStatisticsCentralityTrigger; // contains number of events vs centrality in 1% bins vs trigger
  THnSparseF	*fZvMultCent; // Zv:Mult:Cent
  TH1F        *fTriggerStatistics; // contains number of events per trigger
  TH1F        *fCharge; // charge distribution in data
  TH1F        *fMCCharge; // charge distribution in MC
  THnSparseF	*fDCAPtAll; //control histo: DCAz:DCAxy:pT:eta:phi for all reconstructed tracks
  THnSparseF	*fDCAPtAccepted; //control histo: DCAz:DCAxy:pT:eta:phi for all accepted reco tracks
  THnSparseF	*fMCDCAPtSecondary; //control histo: DCAz:DCAxy:pT:eta:phi for all accepted reco track, which are secondaries (using MC info)
  THnSparseF	*fMCDCAPtPrimary; //control histo: DCAz:DCAxy:pT:eta:phi for all accepted reco track, which are primaries (using MC info)
  THnF        *fCrossCheckAll[5]; //control histo: {CrossedRows,Ncluster,Chi,Length,CrossedRows/Findable} vs pT,eta,phi,Centrality for all tracks
  THnF        *fCrossCheckAcc[5]; //control histo: {CrossedRows,Ncluster,Chi,Length,CrossedRows/Findable} vs pT,eta,phi,Centrality after cuts
  TH1F        *fCutPercClusters; // control histo: number of clusters, where the relative cut has been set e-by-e
  TH1F        *fCutPercCrossed; // control histo: number of crossed rows, where the relative cut has been set e-by-e
  TH3F        *fCrossCheckRowsLength; // control histo: number of crossed rows vs length in TPC vs cent
  TH3F        *fCrossCheckClusterLength; // control histo: number of clusters vs length in TPC vs cent
  TH3F        *fCrossCheckRowsLengthAcc; // control histo: number of crossed rows vs length in TPC for all accepted tracks vs cent
  TH3F        *fCrossCheckClusterLengthAcc; // control histo: number of clusters vs length in TPC for all accepted tracks vs cent
  TH3F        *fCrossCheckPtresLength; // control histo: relative pt resolution in 1/pt vs lenght in TPC vs cent
  TH3F        *fCrossCheckPtresRows; // control histo: relative pt resolution in 1/pt vs number of crossed rows in TPC vs cent
  TH1F        *fCutSettings; // control histo: cut settings
  
  TH1F        *fEventplaneDist; // event plane distribution in phi
  TH2F        *fEventplaneRunDist; // event plane distribution in phi
  TH1F        *fMCEventplaneDist; // MC event plane distribution in phi
  TH2F        *fCorrelEventplaneMCDATA; // correlation between data and MC eventplane
  THnSparseF	*fCorrelEventplaneDefaultCorrected; // correlation between default and corrected (== subtraction of current track) eventplane
  TH2F        *fEventplaneSubtractedPercentage; // percentage of subtracted tracks
  
  TH2F		*fChargeOverPtRuns; // charge/pT vs run
  
  TH2F		*fVZEROMultCentrality; // VZERO Multiplicity vs Centrality
  TH2F		*fVEROMultRefMult; // VZERO Multiplicity vs Reference Multiplicity
  
  // cross check for event plane resolution
  TH2F      *fEPDistCent; // event plane distribution vs centrality
  TH2F		  *fPhiCent;	// particle phi distribution vs centrality
  TProfile	*fPcosEPCent; // < cos 2 psi_ep > vs centrality
  TProfile	*fPsinEPCent; // < sin 2 psi_ep > vs centrality
  TProfile	*fPcosPhiCent; // < cos 2 phi > vs centrality
  TProfile	*fPsinPhiCent; // < sin 2 phi > vs centrality
  TH2F		  *fEPContributionDifference; // difference between contribution of track to EP between own calculation and lookup
  
  // cross check for event plane determination
  TH2F		*fDeltaPhiCent; // DeltaPhi:Cent - DeltaPhi in the range from -pi to pi
  TH2F		*fDeltaPhiSymCent; // DeltaPhi:Cent - DeltaPhi in the range from 0 to pi/2
  
  TH1F		*fMCRecTracksMult; // number of reconstructed tracks vs reference multiplcity
  TH1F		*fMCGenTracksMult; // number of generated tracks vs reference multiplcity
  
  
  
  THnSparseF	*fCrossCheckFilterBitPhiCent; // FilterBit:Phi:Centrality
  
  TH1F		*fTriggerStringsFired; // distribution of fired trigger strings
  TH1F		*fTriggerStringComplete; // complete fired trigger string
  
  // vertex histograms
  TH1F		*fVertexZ; // global vertex Z distribution
  TH1F		*fVertexZSPD; // SPD vertex Z distribution
  TH1F		*fVertexZTPC; // TPC vertex Z distribution
  TH1F		*fDeltaVertexZGlobalSPD; // difference between global and SPD vertex Z position
  TH1F		*fDeltaVertexZGlobalTPC; // difference between global and TPC vertex Z position
  TH1F		*fVertexContributors; // Ncontributors to vertex
  
  TH1F		*fVertexZAfterCuts; // global vertex Z distribution after cuts on VertexZ and NContrib
  TH1F		*fVertexZSPDAfterCuts; // SPD vertex Z distribution after cuts on VertexZ and NContrib
  TH1F		*fVertexZTPCAfterCuts; // TPC vertex Z distribution after cuts on VertexZ and NContrib
  TH1F		*fDeltaVertexZGlobalSPDAfterCuts; // difference between global and SPD vertex Z position after cuts on VertexZ and NContrib
  TH1F		*fDeltaVertexZGlobalTPCAfterCuts; // difference between global and TPC vertex Z position after cuts on VertexZ and NContrib
  TH1F		*fVertexContributorsAfterCuts; // Ncontributors to vertex after cuts on VertexZ and NContrib
  TH1F		*fVertexContributorsAfterCutsCent; // Ncontributors to vertex after cuts on VertexZ and NContrib for central triggers
  TH1F		*fVertexContributorsAfterCutsSemi; // Ncontributors to vertex after cuts on VertexZ and NContrib for semicentral triggers
  TH1F		*fVertexContributorsAfterCutsMB; // Ncontributors to vertex after cuts on VertexZ and NContrib for MB triggers
  
  // global variables
  Bool_t fIsMonteCarlo;
  Int_t  fEventNumberForPtSpectra;
  Bool_t fFillEventPtSpectraHistogram;
  
  TString fEPselector;
  TString fCentEstimator;
  Bool_t  fDoMinBiasAnalysis;
  TString fDisabledTriggerString;
  
  // event cut variables
  Double_t fCutMaxZVertex;
  Int_t	   fVertexMinContributors; // minimum contributors to vertex
  
  // track kinematic cut variables
  Double_t fCutPtMin;
  Double_t fCutPtMax;
  Double_t fCutEtaMin;
  Double_t fCutEtaMax;
  
  Double_t fAnchorPointCorrectionFactor; // factor to check for systematic effect of systematic uncertainty of setting of anchor point
  Bool_t fDoAnchorPointSystStudy; // do systematic study on anchor point
  
  // track quality cut variables
  Int_t	    fFilterBit;
  Bool_t		fHybridTracking;
  Bool_t 	  fUseRelativeCuts;
  Bool_t  	fCutRequireTPCRefit;
  Bool_t 	  fCutRequireITSRefit;
  Double_t	fCutMinNumberOfClusters;
  Double_t	fCutPercMinNumberOfClusters;
  Double_t 	fCutMinNumberOfCrossedRows;
  Double_t 	fCutPercMinNumberOfCrossedRows;
  Double_t 	fCutMinRatioCrossedRowsOverFindableClustersTPC;
  Double_t 	fCutMaxChi2PerClusterTPC;
  Double_t 	fCutMaxFractionSharedTPCClusters;
  Double_t 	fCutMaxDCAToVertexZ;
  Double_t 	fCutMaxDCAToVertexXY;
  Double_t 	fCutMaxChi2PerClusterITS;
  Bool_t  	fCutDCAToVertex2D;
  Bool_t 	  fCutRequireSigmaToVertex;
  Double_t 	fCutMaxDCAToVertexXYPtDepPar0;
  Double_t 	fCutMaxDCAToVertexXYPtDepPar1;
  Double_t 	fCutMaxDCAToVertexXYPtDepPar2;
  Bool_t 	  fCutAcceptKinkDaughters;
  Double_t 	fCutMaxChi2TPCConstrainedGlobal;
  Bool_t	  fCutLengthInTPCPtDependent;
  Double_t	fPrefactorLengthInTPCPtDependent;
  Bool_t 		fCrossCheckCorrelHisto; //
  
  //binning for THNsparse
  Int_t       fMultNbins;
  Int_t       fMultFineNbins;
  Int_t       fPtNbins;
  Int_t       fPtCorrNbins;
  Int_t       fPtCheckNbins;
  Int_t       fEtaNbins;
  Int_t       fEtaCheckNbins;
  Int_t       fZvNbins;
  Int_t       fCentralityNbins;
  Int_t       fPhiNbins;
  Int_t       fPhiCorrNbins;
  Int_t       fDeltaphiNbins;
  Int_t       fRunNumberNbins;
  Double_t*   fBinsMult; //[fMultNbins]
  Double_t*   fBinsMultFine; //[fMultFineNbins]
  Double_t*   fBinsPt; //[fPtNbins]
  Double_t*   fBinsPtCorr; //[fPtCorrNbins]
  Double_t*   fBinsPtCheck; //[fPtCheckNbins]
  Double_t*   fBinsEta; //[fEtaNbins]
  Double_t*   fBinsEtaCheck; //[fEtaCheckNbins]
  Double_t*   fBinsZv; //[fZvNbins]
  Double_t*   fBinsCentrality; //[fCentralityNbins]
  Double_t*   fBinsPhi; //[fPhiNbins]
  Double_t*   fBinsPhiCorr; //[fPhiCorrNbins]
  Double_t*   fBinsDeltaphi; //[fDeltaphiNbins]
  Double_t*	  fBinsRunNumber; //[fRunNumberNbins]
  
  AlidNdPtAnalysisPbPbAOD(const AlidNdPtAnalysisPbPbAOD&); // not implemented
  AlidNdPtAnalysisPbPbAOD& operator=(const AlidNdPtAnalysisPbPbAOD&); // not implemented
  
  ClassDef(AlidNdPtAnalysisPbPbAOD,21); // has to be at least 1, otherwise not streamable...
};

#endif
