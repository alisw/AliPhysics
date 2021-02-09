#ifndef AliAnalysisTaskPhiCorrelations_H
#define AliAnalysisTaskPhiCorrelations_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//
// Analysis class for Underlying Event studies w.r.t. leading track
//
// Look for correlations on the tranverse regions w.r.t
// the leading track in the event
//
// This class needs input AODs.
// The output is a list of analysis-specific containers.
//
// The AOD can be either connected to the InputEventHandler  
// for a chain of AOD files 
// or 
// to the OutputEventHandler
// for a chain of ESD files,
// in this case the class should be in the train after the jet-finder
//
//    Authors:
//    Jan Fiete Grosse-Oetringhaus
// 
////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"
#include "AliUEHist.h"
#include "TString.h"
#include "AliBasicParticle.h"
#include "AliLog.h"
#include "THn.h" // in cxx file causes .../THn.h:257: error: conflicting declaration ‘typedef class THnT<float> THnF’

class AliAODEvent;
class AliAnalyseLeadingTrackUE;
class AliInputEventHandler;
class AliMCEvent;
class AliMCEventHandler;
class AliUEHistograms;
class TH1;
class TObjArray;
class AliEventPoolManager;
class AliESDEvent;
class AliHelperPID;
class AliAnalysisUtils;
class TFormula;
class TMap;
class AliGenEventHeader;
class AliVEvent;

using std::vector;


class  AliAnalysisTaskPhiCorrelations : public AliAnalysisTask
{
public:
  AliAnalysisTaskPhiCorrelations(const char* name="AliAnalysisTaskPhiCorrelations");
  virtual ~AliAnalysisTaskPhiCorrelations();


  // Implementation of interace methods
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *option);
  virtual void FinishTaskOutput();

  // Setters/Getters
  // general configuration
  virtual void SetDebugLevel(Int_t level) { fDebug = level; }
  virtual void SetMode(Int_t mode) { fMode  = mode;  }
  virtual void SetReduceMemoryFootprint(Bool_t flag) { fReduceMemoryFootprint = flag; }
  virtual void SetEventMixing(Bool_t flag) { fFillMixed = flag; }
  virtual void SetMixingTracks(Int_t tracks) { fMixingTracks = tracks; }
  virtual void SetTwoTrackEfficiencyStudy(Bool_t flag) { fTwoTrackEfficiencyStudy = flag; }
  virtual void SetTwoTrackEfficiencyCut(Float_t value = 0.02, Float_t min = 0.8) { fTwoTrackEfficiencyCut = value; fTwoTrackCutMinRadius = min; }
  virtual void SetUseVtxAxis(Int_t flag) { fUseVtxAxis = flag; }
  virtual void SetCourseCentralityBinning(Bool_t flag) { fCourseCentralityBinning = flag; }
  virtual void SetSkipTrigger(Bool_t flag) { fSkipTrigger = flag; }
  virtual void SetInjectedSignals(Bool_t flag) { fInjectedSignals = flag; }
  virtual void SetV0CL1PileUp(Bool_t flag) { fV0CL1PileUp = flag; }
  virtual void SetESDTPCTrackPileUp(Bool_t flag) { fESDTPCTrackPileUp = flag; }
  virtual void SetTPCITSTOFPileUp(Bool_t flag) { fTPCITSTOFPileUp = flag; }
  void SetRandomizeReactionPlane(Bool_t flag) { fRandomizeReactionPlane = flag; }

  // histogram settings
  void SetEfficiencyCorrectionTriggers(THnF* hist) { fEfficiencyCorrectionTriggers = hist; }
  void SetEfficiencyCorrectionAssociated(THnF* hist) { fEfficiencyCorrectionAssociated = hist; }
  void SetCentralityWeights(TH1* hist) { fCentralityWeights = hist; }
  void SetCentralityMCGen_V0M(TH1* hist) { fCentralityMCGen_V0M = hist; }
  void SetCentralityMCGen_CL1(TH1* hist) { fCentralityMCGen_CL1 = hist; }
  // for event QA
  void SetTracksInVertex(Int_t val) { fnTracksVertex = val; }
  void SetZVertex(Double_t val) { fZVertex = val; }
  void SetAcceptOnlyMuEvents(Bool_t val) { fAcceptOnlyMuEvents = val; }

  // track cuts
  void SetTrackEtaCut(Double_t val) { fTrackEtaCut = val; }
  void SetTrackEtaCutMin(Double_t val) { fTrackEtaCutMin = val; }
  void SetOnlyOneEtaSide(Int_t flag) { fOnlyOneEtaSide = flag; }
  void SetOnlyOneAssocEtaSide(Int_t flag) { fOnlyOneAssocEtaSide = flag; }
  void SetTrackPhiCutEvPlMin(Double_t val) { fTrackPhiCutEvPlMin = val; }
  void SetTrackPhiCutEvPlMax(Double_t val) { fTrackPhiCutEvPlMax = val; }
  void SetPtMin(Double_t val) { fPtMin = val; }
  void SetFilterBit(UInt_t val) { fFilterBit = val;  }
  void SetDCAXYCut(TFormula* value) { fDCAXYCut = value; }
  void SetSharedClusterCut(Float_t value) { fSharedClusterCut = value; }
  void SetCrossedRowsCut(Int_t value) { fCrossedRowsCut = value; }
  void SetFoundFractionCut(Double_t value) { fFoundFractionCut = value; }
  void SetTrackStatus(UInt_t status) { fTrackStatus = status; }
  void SetCheckMotherPDG(Bool_t checkpdg) { fCheckMotherPDG = checkpdg; }

  // track cuts
  void SetTrackletDphiCut(Double_t val) { fTrackletDphiCut = val; }

  void SetEventSelectionBit(UInt_t val) { fSelectBit = val;  }
  void SetUseChargeHadrons(Bool_t val) { fUseChargeHadrons = val; }
  void SetSelectParticleSpecies(Int_t trigger, Int_t associated) { fParticleSpeciesTrigger = trigger; fParticleSpeciesAssociated = associated; }
  void SetSelectCharge(Int_t selectCharge) { fSelectCharge = selectCharge; }
  void SetSelectTriggerCharge(Int_t selectCharge) { fTriggerSelectCharge = selectCharge; }
  void SetSelectAssociatedCharge(Int_t selectCharge) { fAssociatedSelectCharge = selectCharge; }
  void SetTriggerRestrictEta(Float_t eta) { fTriggerRestrictEta = eta; }
  void SetEtaOrdering(Bool_t flag) { fEtaOrdering = flag; }
  void SetPairCuts(Float_t conversions = 0.004, Float_t resonances = 0.005) { fCutConversionsV = conversions; fCutResonancesV = resonances; fCutOnLambdaV = resonances; fCutOnK0sV = resonances; }
  void SetCustomCut(Float_t cutOnCustomMass, Float_t cutOnCustomFirst, Float_t cutOnCustomSecond, Float_t cutOnCustomV) { fCutOnCustomMass = cutOnCustomMass; fCutOnCustomFirst = cutOnCustomFirst; fCutOnCustomSecond = cutOnCustomSecond; fCutOnCustomV = cutOnCustomV; }
  void SetCutOnPhi(bool cutOnPhi) { fCutOnPhi = cutOnPhi; }
  void SetCutOnPhi(Float_t cutOnPhi) { fCutOnPhiV = cutOnPhi; }
  void SetCutOnRho(bool cutOnRho) { fCutOnRho = cutOnRho; }
  void SetCutOnRho(Float_t cutOnRho) { fCutOnRhoV = cutOnRho; }
  void SetCutOnLambda(Float_t cutOnLambda) { fCutOnLambdaV = cutOnLambda; }
  void SetCutOnK0s(Float_t cutOnK0s) { fCutOnK0sV = cutOnK0s; }
  void SetRejectResonanceDaughters(Int_t value) { fRejectResonanceDaughters = value; }
  void SetCentralityMethod(const char* method) { fCentralityMethod = method; }
  void SetFillpT(Bool_t flag) { fFillpT = flag; }
  void SetStepsFillSkip(Bool_t step0, Bool_t step6, Bool_t step9 = kTRUE) { fFillOnlyStep0 = step0; fSkipStep6 = step6; fSkipStep9 = step9; }
  void SetRejectCentralityOutliers(Bool_t flag = kTRUE) { fRejectCentralityOutliers = flag; }
  void SetRejectZeroTrackEvents(Bool_t flag) { fRejectZeroTrackEvents = flag; }
  void SetRemoveWeakDecays(Bool_t flag = kTRUE) { fRemoveWeakDecays = flag; }
  void SetRemoveDuplicates(Bool_t flag = kTRUE) { fRemoveDuplicates = flag; }
  void SetSkipFastCluster(Bool_t flag = kTRUE) { fSkipFastCluster = flag; }
  void SetWeightPerEvent(Bool_t flag = kTRUE) { fWeightPerEvent = flag; }
  void SetCustomBinning(const char* binningStr) { fCustomBinning = binningStr; }
  void SetPtOrder(Bool_t flag) { fPtOrder = flag; }
  void SetTriggersFromDetector(Int_t flag) { fTriggersFromDetector = flag; }
  void SetAssociatedFromDetector(Int_t flag) { fAssociatedFromDetector = flag; }
  void SetUseUncheckedCentrality(Bool_t flag) { fUseUncheckedCentrality = flag; }
  void SetCheckCertainSpecies(Int_t species) { fCheckCertainSpecies = species; }
  void SetRemoveWeakDecaysInMC(Bool_t flag) { fRemoveWeakDecaysInMC = flag; }
  void SetFillYieldRapidity(Bool_t flag) { fFillYieldRapidity = flag; }
  void SetFillCorrelationsRapidity(Bool_t flag) { fFillCorrelationsRapidity = flag; }
  void SetUseDoublePrecision(Bool_t flag) { fUseDoublePrecision = flag; }
  void SetUseNewCentralityFramework(Bool_t flag) { fUseNewCentralityFramework = flag; }

  AliHelperPID* GetHelperPID() { return fHelperPID; }
  void SetHelperPID(AliHelperPID* pid) { fHelperPID = pid; }

  AliAnalysisUtils* GetAnalysisUtils() { return fAnalysisUtils; }
  void SetAnalysisUtils(AliAnalysisUtils* utils) { fAnalysisUtils = utils; }

  TMap* GetMap() { return fMap; }
  void SetMap(TMap* map) { fMap = map; }

  // configuration for tracks with jets removed
  void SetJetBranchName(const char* branchName) { fJetBranchName = branchName; }
  const char* GetJetBranchName() const { return fJetBranchName; }
  void SetTrackEtaMax(Float_t eta) { fTrackEtaMax = eta; }
  Float_t GetTrackEtaMax() const { return fTrackEtaMax; }
  void SetJetEtaMax(Float_t eta) { fJetEtaMax = eta; }
  Float_t GetJetEtaMax() const { return fJetEtaMax; }
  void SetJetPtMin(Float_t pt) { fJetPtMin = pt; }
  Float_t GetJetPtMin() const { return fJetPtMin; }
  void SetJetConstMin(Int_t constMin) { fJetConstMin = constMin; }
  Int_t GetJetConstMin() const { return fJetConstMin; }
  void SetExclusionRadius(Float_t r) { fExclusionRadius = r; }
  Float_t GetExclusionRadius() const { return fExclusionRadius; }

  // Custom particle arrays
  void SetCustomParticlesA(const char* customA) { fCustomParticlesA = customA; }
  void SetCustomParticlesB(const char* customB) { fCustomParticlesB = customB; }

  // ##### External event pool configuration
  void SetExternalEventPoolManager(AliEventPoolManager* mgr) {fPoolMgr = mgr;}
  AliEventPoolManager* GetEventPoolManager() {return fPoolMgr;}
  void SetUsePtBinnedEventPool(Bool_t val) {fUsePtBinnedEventPool = val;}
  void SetCheckEventNumberInMixedEvent(Bool_t val) {fCheckEventNumberInMixedEvent = val;}

  // Set which pools will be saved
  void AddEventPoolsToOutput(Double_t minCent, Double_t maxCent,  Double_t minZvtx, Double_t maxZvtx, Double_t minPt, Double_t maxPt);

private:
  AliAnalysisTaskPhiCorrelations(const AliAnalysisTaskPhiCorrelations &det);
  AliAnalysisTaskPhiCorrelations& operator=(const  AliAnalysisTaskPhiCorrelations &det);
  void AddSettingsTree(); // add list of settings to output list
  // Analysis methods
  void AnalyseCorrectionMode(); // main algorithm to get correction maps
  void AnalyseDataMode();       // main algorithm to get raw distributions
  void Initialize();            // initialize some common pointer
  Double_t GetCentrality(AliVEvent* inputEvent, TObject* mc);
  TObjArray* CloneAndReduceTrackList(TObjArray* tracks, Double_t minPt = 0., Double_t maxPt = -1.);
  void RemoveDuplicates(TObjArray* tracks);
  void CleanUp(TObjArray* tracks, TObject* mcObj, Int_t maxLabel);
  void RemoveWeakDecaysInMC(TObjArray* tracks, TObject* mcObj);
  void SelectCharge(TObjArray* tracks);
  AliGenEventHeader* GetFirstHeader();
  Bool_t AcceptEventCentralityWeight(Double_t centrality);
  void ShiftTracks(TObjArray* tracks, Double_t angle);
  TObjArray* GetParticlesFromDetector(AliVEvent* inputEvent, Int_t idet);
  Bool_t IsMuEvent();
  Bool_t InitiateEventPlane(Double_t& evtPlanePhi, AliVEvent* inputEvent);
  Long64_t GetUniqueEventID(AliVEvent* inputEvent);

  // General configuration
  Int_t fDebug; // Debug flag
  Int_t fMode;  // fMode = 0: data-like analysis, fMode = 1: corrections analysis
  Bool_t  fReduceMemoryFootprint;   // reduce memory consumption by writing less debug histograms
  Bool_t  fFillMixed;               // enable event mixing (default: ON)
  Int_t   fMixingTracks;            // size of track buffer for event mixing
  Bool_t  fTwoTrackEfficiencyStudy; // two-track efficiency study on
  Float_t fTwoTrackEfficiencyCut;   // enable two-track efficiency cut
  Float_t fTwoTrackCutMinRadius;    // minimum radius for two-track efficiency cut
  Int_t   fUseVtxAxis;              // use z vtx as axis (needs 7-10 times more memory!)
  Bool_t  fCourseCentralityBinning; // less centrality bins
  Bool_t  fSkipTrigger;             // skip trigger selection
  Bool_t  fInjectedSignals;         // check header to skip injected signals in MC
  Bool_t  fRandomizeReactionPlane;  // change the orientation of the RP by a random value by shifting all tracks
  Bool_t  fV0CL1PileUp;             // Remove pile up events according to V0 CL1 correlation
  Bool_t  fESDTPCTrackPileUp;       // Remove pile up events according ESD tracks and number of TPC only tracks correlation
  Bool_t  fTPCITSTOFPileUp;         // Remove pile up events according to TPC+ITS tracks (FilterBit 32) and number of tracks matched to TOF with BCid=0 correlation

  AliHelperPID* fHelperPID;         // points to class for PID
  AliAnalysisUtils* fAnalysisUtils; // points to class with common analysis utilities
  TMap* fMap;                       // points to TMap class containing scaling factors for VZERO A signal

  // Pointers to external UE classes
  AliAnalyseLeadingTrackUE* fAnalyseUE; //! points to class containing common analysis algorithms
  AliUEHistograms* fHistos;             //! points to class to handle histograms/containers
  AliUEHistograms* fHistosMixed;        //! points to class to handle mixed histograms/containers

  THnF* fEfficiencyCorrectionTriggers;   // if non-0 this efficiency correction is applied on the fly to the filling for trigger particles. The factor is multiplicative, i.e. should contain 1/efficiency. Axes: eta, pT, centrality, z-vtx
  THnF* fEfficiencyCorrectionAssociated; // if non-0 this efficiency correction is applied on the fly to the filling for associated particles. The factor is multiplicative, i.e. should contain 1/efficiency. Axes: eta, pT, centrality, z-vtx
  TH1* fCentralityWeights;               // for centrality flattening
  TH1* fCentralityMCGen_V0M;             // for centrality from generated MCGen_V0M
  TH1* fCentralityMCGen_CL1;             // for centrality from generated MCGen_CL1

  // Handlers and events
  AliAODEvent*             fAOD;             //! AOD Event
  AliESDEvent*             fESD;             //! ESD Event
  TClonesArray*            fArrayMC;         //! Array of MC particles
  AliInputEventHandler*    fInputHandler;    //! Generic InputEventHandler
  AliMCEvent*              fMcEvent;         //! MC event
  AliInputEventHandler*    fMcHandler;       //! MCEventHandler
  AliEventPoolManager*     fPoolMgr;         // event pool manager

  // Histogram settings
  TList*              fListOfHistos;    //  Output list of containers

  // Event QA cuts
  Int_t               fnTracksVertex;        // QA tracks pointing to principal vertex
  Double_t            fZVertex;              // Position of Vertex in Z direction
  Bool_t              fAcceptOnlyMuEvents;   // Only Events with at least one muon are accepted
  TString             fCentralityMethod;     // Method to determine centrality

  // Track cuts
  Double_t            fTrackEtaCut;          // Maximum Eta cut on particles
  Double_t            fTrackEtaCutMin;       // Minimum Eta cut on particles
  Double_t            fTrackPhiCutEvPlMin;   // Minimum Phi cut on particles with respect to the Event Plane (values between 0 and Pi/2)
  Double_t            fTrackPhiCutEvPlMax;   // Maximum Phi cut on particles with respect to the Event Plane (values between 0 and Pi/2), if = 0, then no cut is performed
  Int_t               fOnlyOneEtaSide;       // decides that only trigger particle from one eta side are considered (0 = all; -1 = negative, 1 = positive)
  Int_t               fOnlyOneAssocEtaSide;  // decides that only associated particle from one eta side are considered (0 = all; -1 = negative, 1 = positive)
  Double_t            fPtMin;                // Min pT to start correlations
  TFormula*           fDCAXYCut;             // additional pt dependent cut on DCA XY (only for AOD)
  Double_t            fSharedClusterCut;     // cut on shared clusters (only for AOD)
  Int_t               fCrossedRowsCut;       // cut on crossed rows (only for AOD)
  Double_t            fFoundFractionCut;     // cut on crossed rows/findable clusters (only for AOD)
  UInt_t              fFilterBit;            // Select tracks from an specific track cut
  UInt_t              fTrackStatus;          // if non-0, the bits set in this variable are required for each track
  UInt_t              fSelectBit;            // Select events according to AliAnalysisTaskJetServices bit maps
  Bool_t              fUseChargeHadrons;     // Only use charge hadrons
  Int_t               fParticleSpeciesTrigger; // Select which particle to use for the trigger [ -1 (all, default) 0 (pions) 1 (kaons) 2 (protons) 3 (others) particles ]
  Int_t               fParticleSpeciesAssociated; // Select which particle to use for the associated [ -1 (all, default) 0 (pions) 1 (kaons) 2 (protons) 3 (others) particles ]
  Bool_t              fCheckMotherPDG;       // Check the PDG code of mother for secondaries

  // Tracklets cuts
  Double_t            fTrackletDphiCut;      // maximum Dphi cut on tracklets

  Int_t fSelectCharge;           // (un)like sign selection when building correlations: 0: no selection; 1: unlike sign; 2: like sign
  Int_t fTriggerSelectCharge;    // select charge of trigger particle: 1: positive; -1 negative
  Int_t fAssociatedSelectCharge; // select charge of associated particle: 1: positive; -1 negative
  Float_t fTriggerRestrictEta;   // restrict eta range for trigger particle (default: -1 [off])
  Bool_t fEtaOrdering;           // eta ordering, see AliUEHistograms.h for documentation
  Float_t fCutConversionsV;      // cut on conversions (inv mass)
  Float_t fCutResonancesV;       // cut on resonances (inv mass)
  Float_t fCutOnLambdaV;         // cut on Lambda (inv mass)
  Float_t fCutOnK0sV;            // cut on K0s (inv mass)
  Bool_t fCutOnPhi;              // cut on Phi
  Float_t fCutOnPhiV;            // cut on Phi (inv mass)
  Bool_t fCutOnRho;              // cut on Rho
  Float_t fCutOnRhoV;            // cut on Rho (inv mass)
  Float_t fCutOnCustomMass;      // user-defined inv mass value
  Float_t fCutOnCustomFirst;     // user-defined mass of the 1st particle
  Float_t fCutOnCustomSecond;    // user-defined mass of the 2nd particle
  Float_t fCutOnCustomV;         // cut on user-defined value (inv mass)
  Int_t fRejectResonanceDaughters; // reject all daughters of all resonance candidates (1: test method (cut at m_inv=0.9); 2: k0; 3: lambda)
  Bool_t fFillOnlyStep0;         // fill only step 0
  Bool_t fSkipStep6;             // skip step 6 when filling
  Bool_t fSkipStep9;             // skip step 9 when filling
  Bool_t fRejectCentralityOutliers; // enable rejection of outliers in centrality vs no track correlation. Interferes with the event plane dependence code
  Bool_t fRejectZeroTrackEvents; // reject events which have no tracks (using the eta, pT cuts defined)
  Bool_t fRemoveWeakDecays;      // remove secondaries from weak decays from tracks and particles
  Bool_t fRemoveDuplicates;      // remove particles with the same label (double reconstruction)
  Bool_t fSkipFastCluster;       // skip kFastOnly flagged events (only for data)
  Bool_t fWeightPerEvent;        // weight with the number of trigger particles per event
  TString fCustomBinning;        // supersedes default binning if set, see AliUEHist::GetBinning or AliUEHistograms::AliUEHistograms for syntax and examples
  Bool_t fPtOrder;               // apply pT,a < pt,t condition; default: kTRUE
  Int_t fTriggersFromDetector;   // 0 = tracks (default); 1 = VZERO_A; 2 = VZERO_C; 3 = SPD tracklets; 4 = forward muons; 5 = tracks w/o jets; 6, 7 = arbitrary AliVParticle-TClonesArrays (see SetCustomParticleArrayA())
  Int_t fAssociatedFromDetector; // 0 = tracks (default); 1 = VZERO_A; 2 = VZERO_C; 3 = SPD tracklets; 4 = forward muons; 5 = tracks w/o jets; 6,7 = arbitrary AliVParticle-TClonesArrays (see SetCustomParticleArrayB())
  Bool_t fUseUncheckedCentrality;// use unchecked centrality; default: kFALSE
  Int_t fCheckCertainSpecies;    // make eta,pt distribution of MC particles with the label of fCheckCertainSpecies, by default switched off (value: -1)
  Bool_t fRemoveWeakDecaysInMC;  // remove weak decays which have been included by mistake as primaries in the stack (bug in AMPT)
  Bool_t fFillYieldRapidity;     // fill a control histogram centrality vs pT vs y
  Bool_t fFillCorrelationsRapidity; // fills correlation histograms with rapidity instead of pseudorapidity (default: kFALSE)
  Bool_t fUseDoublePrecision;    // use double precision for AliTHn
  Bool_t fUseNewCentralityFramework; // use the AliMultSelection framework

  Bool_t fFillpT;                // fill sum pT instead of number density

  // configuration for tracks with jet removal
  TString fJetBranchName;        // name of jet branch for exclusion of in-jet tracks
  Float_t fTrackEtaMax;          // maximum eta to accept track
  Float_t fJetEtaMax;            // maximum eta to accept jet
  Float_t fJetPtMin;             // minimum pt to accept jet
  Int_t   fJetConstMin;          // minimum number of constituents
  Float_t fExclusionRadius;      // radius around jet in which tracks are rejected

  // Variables for arbitrary triggers/associates
  TString fCustomParticlesA;   // name of TClonesArray of AliVParticle-derived particles attached to fInputEvent
  TString fCustomParticlesB;   // name of TClonesArray of AliVParticle-derived particles attached to fInputEvent

  // Event pool variables
  vector<vector<Double_t> > fEventPoolOutputList; // vector representing a list of pools (given by value range) that will be saved
  Bool_t fUsePtBinnedEventPool;                   // uses event pool in pt bins
  Bool_t fCheckEventNumberInMixedEvent;           // check event number before correlation in mixed event

  ClassDef(AliAnalysisTaskPhiCorrelations, 63); // Analysis task for delta phi correlations
};

#endif
