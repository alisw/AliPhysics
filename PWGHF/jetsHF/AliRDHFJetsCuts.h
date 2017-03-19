#ifndef ALIRDHFJETSCUTS_H
#define ALIRDHFJETSCUTS_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliRDHFJetsCuts
// base class for cuts on AOD reconstructed heavy-flavour decays
// Author: A. Rossi, andrea.rossi@cern.ch E. Bruna, elena.bruna@to.infn.it
//***********************************************************

#include <TString.h>

#include "AliAnalysisCuts.h"
#include "AliESDtrackCuts.h"
#include "AliAODPidHF.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliEmcalJet.h"

class AliAODTrack;
class AliAODRecoDecayHF;
class AliESDVertex;

class AliRDHFJetsCuts : public AliAnalysisCuts {
public:

  enum ECentrality {kCentOff, kCentV0M, kCentTRK, kCentTKL, kCentCL1, kCentInvalid};
  enum EPileup {kNoPileupSelection, kRejectPileupEvent, kRejectTracksFromPileupVertex};
  enum ERejBits {kNotSelTrigger, kNoVertex, kTooFewVtxContrib, kZVtxOutFid, kPileupSPD, kOutsideCentrality, kPhysicsSelection, kBadSPDVertex, kZVtxSPDOutFid, kCentralityFlattening};
  enum EV0sel  {kAllV0s = 0, kOnlyOfflineV0s = 1, kOnlyOnTheFlyV0s = 2};

  AliRDHFJetsCuts(const Char_t* name = "RDHFJetsCuts", const Char_t* title = "");
  virtual ~AliRDHFJetsCuts();
  AliRDHFJetsCuts(const AliRDHFJetsCuts& source);
  AliRDHFJetsCuts& operator=(const AliRDHFJetsCuts& source);

  void SetMinCentrality(Float_t minCentrality = 0.) {fMinCentrality = minCentrality;}
  void SetMaxCentrality(Float_t maxCentrality = 100.) {fMaxCentrality = maxCentrality;}
  void SetMinVtxType(Int_t type = 3) {fMinVtxType = type;}
  void SetUseEventsWithOnlySPDVertex(Bool_t flag = kTRUE)
  {
    if (flag) fMinVtxType = 1;
    else fMinVtxType = 3;
  }
  void SetMinVtxContr(Int_t contr = 1) {fMinVtxContr = contr;}
  void SetMaxVtxRdChi2(Float_t chi2 = 1e6) {fMaxVtxRedChi2 = chi2;}
  void SetMaxVtxZ(Float_t z = 1e6) {fMaxVtxZ = z;}
  void SetMinSPDMultiplicity(Int_t mult = 0) {fMinSPDMultiplicity = mult;}

  void SetTriggerMask(ULong64_t mask = 0) {fTriggerMask = mask;}
  void SetUseAnyTrigger() {fTriggerMask = AliVEvent::kAny;}
  void EnableMBTrigger()
  {
    fTriggerMask |= AliVEvent::kMB;
    fUseOnlyOneTrigger = kFALSE;
  }
  void ResetMaskAndEnableMBTrigger()
  {
    fTriggerMask = AliVEvent::kMB;
    fUseOnlyOneTrigger = kFALSE;
  }
  void SetUseMBTriggerExclusively()
  {
    fTriggerMask = AliVEvent::kMB;
    fUseOnlyOneTrigger = kTRUE;
  }
  void EnableCentralTrigger()
  {
    fTriggerMask |= AliVEvent::kCentral;
    fUseOnlyOneTrigger = kFALSE;
  }
  void ResetMaskAndEnableCentralTrigger()
  {
    fTriggerMask = AliVEvent::kCentral;
    fUseOnlyOneTrigger = kFALSE;
  }
  void SetUseCentralTriggerExclusively()
  {
    fTriggerMask = AliVEvent::kCentral;
    fUseOnlyOneTrigger = kTRUE;
  }
  void EnableSemiCentralTrigger()
  {
    fTriggerMask |= AliVEvent::kSemiCentral;
    fUseOnlyOneTrigger = kFALSE;
  }
  void ResetMaskAndEnableSemiCentralTrigger()
  {
    fTriggerMask = AliVEvent::kSemiCentral;
    fUseOnlyOneTrigger = kFALSE;
  }
  void SetUseSemiCentralTriggerExclusively()
  {
    fTriggerMask = AliVEvent::kSemiCentral;
    fUseOnlyOneTrigger = kTRUE;
  }
  void EnableEMCALTrigger()
  {
    fTriggerMask |= (AliVEvent::kEMCEJE | AliVEvent::kEMCEGA);
    fUseOnlyOneTrigger = kFALSE;
  }
  void ResetMaskAndEnableEMCALTrigger()
  {
    fTriggerMask = (AliVEvent::kEMCEJE | AliVEvent::kEMCEGA);
    fUseOnlyOneTrigger = kFALSE;
  }
  void SetUseEMCALTriggerExclusively()
  {
    fTriggerMask = (AliVEvent::kEMCEJE | AliVEvent::kEMCEGA);
    fUseOnlyOneTrigger = kTRUE;
  }
  void SetRemoveTrackletOutliers(Bool_t opt) {fRemoveTrackletOutliers = opt;}
  void SetCutOnzVertexSPD(Int_t opt)
  {
    if (opt >= 0 && opt <= 2) fCutOnzVertexSPD = opt;
    else AliError("Wrong option for cut on zVertexSPD");
  }
  void SetTriggerClass(TString trclass0, TString trclass1 = "") {fTriggerClass[0] = trclass0; fTriggerClass[1] = trclass1;}
  void ApplySPDDeadPbPb2011() {fApplySPDDeadPbPb2011 = kTRUE;}
  void ApplySPDMisalignedCutPP2012() {fApplySPDMisalignedPP2012 = kTRUE;}

  void AddTrackCuts(const AliESDtrackCuts* cuts)
  {delete fTrackCuts; fTrackCuts = new AliESDtrackCuts(*cuts); return;}

  void SetUseAOD049(Bool_t flag = kTRUE) {fUseAOD049 = flag; return;}
  void SetKinkRejection(Bool_t flag = kTRUE) {fKinkReject = flag; return;}
  void SetUseTrackSelectionWithFilterBits(Bool_t flag = kTRUE)
  {
    fUseTrackSelectionWithFilterBits = flag; return;
  }
  void SetUseCentrality(Int_t flag = 1);  // see enum below

  void SetRemoveDaughtersFromPrim(Bool_t removeDaughtersPrim) {fRemoveDaughtersFromPrimary = removeDaughtersPrim;}

  void SetRecomputePrimaryVertex(Bool_t opt) {fRecomputePrimVertex = opt;}

  void SetMinPtJet(Double_t ptJet = -1.) {fMinPtJet = ptJet; return;}
  void SetMaxPtJet(Double_t ptJet = 1000.) {fMaxPtJet = ptJet; return;}
  void SetMaxEtaJet(Double_t etaJet = 0.6) {fMaxEtaJet = etaJet; return;}
  void SetJetRadius(Double_t rad = 0.4) {fJetRadius = rad; return;}


  void SetOptPileup(Int_t opt = 0)
  {
    // see enum below
    fOptPileup = opt;
  }
  void SetHistoForCentralityFlattening(TH1F* h, Double_t minCentr, Double_t maxCentr, Double_t centrRef = 0., Int_t switchTRand = 0);
  void ConfigurePileupCuts(Int_t minContrib = 3, Float_t minDz = 0.6)
  {
    fMinContrPileup = minContrib;
    fMinDzPileup = minDz;
  }

  Double_t GetMaxVtxZ() const {return fMaxVtxZ;}
  Float_t GetCentrality(AliAODEvent* aodEvent) {return GetCentrality(aodEvent, (AliRDHFJetsCuts::ECentrality)fUseCentrality);}
  Float_t GetCentrality(AliAODEvent* aodEvent, AliRDHFJetsCuts::ECentrality estimator);

  AliESDtrackCuts* GetTrackCuts() const {return fTrackCuts;}
  Bool_t  GetUseAOD049() const {return fUseAOD049;}
  Bool_t  GetUseKinkRejection() const {return fKinkReject;}
  Bool_t  GetUseEventsWithOnlySPDVertex() const
  {
    if (fMinVtxType == 1 || fMinVtxType == 2) return kTRUE;
    return kFALSE;
  }
  Bool_t  GetUseTrackSelectionWithFilterBits() const {return fUseTrackSelectionWithFilterBits;}
  Bool_t  GetIsPrimaryWithoutDaughters() const {return fRemoveDaughtersFromPrimary;}
  Bool_t GetOptPileUp() const {return fOptPileup;}
  Int_t GetUseCentrality() const {return fUseCentrality;}
  Float_t GetMinCentrality() const {return fMinCentrality;}
  Float_t GetMaxCentrality() const {return fMaxCentrality;}
  Double_t GetMinPtJet() const {return fMinPtJet;}
  Double_t GetMaxPtJet() const {return fMaxPtJet;}
  Double_t GetMaxEtaJet() const {return fMaxEtaJet;}
  Double_t GetJetRadius() const {return fJetRadius;}

  TH1F* GetHistoForCentralityFlattening() {return fHistCentrDistr;}
  void SetUseCentralityFlatteningInMC(Bool_t opt) {fUseCentrFlatteningInMC = opt;}
  Int_t  IsEventSelectedInCentrality(AliVEvent* event);
  Bool_t IsEventSelectedForCentrFlattening(Float_t centvalue);
  Bool_t IsEventSelected(AliVEvent* event);
  Bool_t AreDaughtersSelected(AliAODRecoDecayHF* rd) const;
  Bool_t IsDaughterSelected(AliAODTrack* track, const AliESDVertex* primary, AliESDtrackCuts* cuts) const;
  void IsDaugElectron(AliAODEvent* aod, const AliAODJet* jet, Bool_t& fFlagElec);
  Bool_t IsJetSelected(const AliEmcalJet* jet) const;

  //  void SetupPID(AliVEvent *event);
  //Bool_t IsSelected(TObject *obj) {return IsSelected(obj,0);}
  Bool_t IsSelected(TObject* obj) {if (!obj) return kTRUE; return kFALSE;}
  Bool_t IsSelected(TList* list) {if (!list) return kTRUE; return kFALSE;}

  //virtual  Int_t IsSelected(TObject* obj,Int_t selectionLevel) {return kFALSE;}

  //virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel,AliAODEvent* /*aod*/)
  //{return IsSelected(obj,selectionLevel);}

  //Int_t PtBin(Double_t pt) const;
  void PrintAll()const;
  void PrintTrigger() const;

  //  virtual Bool_t IsInFiducialAcceptance(Double_t /*pt*/,Double_t /*y*/) const {return kTRUE;}

  void SetWhyRejection(Int_t why) {fWhyRejection = why; return;}
  Int_t GetWhyRejection() const {return fWhyRejection;}
  UInt_t GetEventRejectionBitMap() const {return fEvRejectionBits;}
  Bool_t IsEventRejectedDueToTrigger() const
  {
    return fEvRejectionBits & (1 << kNotSelTrigger);
  }
  Bool_t IsEventRejectedDueToNotRecoVertex() const
  {
    return fEvRejectionBits & (1 << kNoVertex);
  }
  Bool_t IsEventRejectedDueToVertexContributors() const
  {
    return fEvRejectionBits & (1 << kTooFewVtxContrib);
  }
  Bool_t IsEventRejectedDueToZVertexOutsideFiducialRegion() const
  {
    return fEvRejectionBits & (1 << kZVtxOutFid);
  }
  Bool_t IsEventRejectedDueToPileupSPD() const
  {
    return fEvRejectionBits & (1 << kPileupSPD);
  }
  Bool_t IsEventRejectedDueToCentrality() const
  {
    return fEvRejectionBits & (1 << kOutsideCentrality);
  }
  Bool_t IsEventRejectedDuePhysicsSelection() const
  {
    return fEvRejectionBits & (1 << kPhysicsSelection);
  }


  void SetFixRefs(Bool_t fix = kTRUE) {fFixRefs = fix; return;}
  void SetUsePhysicsSelection(Bool_t use = kTRUE) {fUsePhysicsSelection = use; return;}
  Bool_t GetUsePhysicsSelection() const { return fUsePhysicsSelection; }



  Bool_t CompareCuts(const AliRDHFJetsCuts* obj) const;


  void SetUseMCVertex() { fUseMCVertex = kTRUE; }
  Bool_t GetUseMCVertex() const { return fUseMCVertex; }

  Bool_t RecalcOwnPrimaryVtx(AliAODRecoDecayHF* d, AliAODEvent* aod) const;
  Bool_t SetMCPrimaryVtx(AliAODRecoDecayHF* d, AliAODEvent* aod) const;
  void   CleanOwnPrimaryVtx(AliAODRecoDecayHF* d, AliAODEvent* aod, AliAODVertex* origownvtx) const;

  Bool_t CountEventForNormalization() const
  { if (fWhyRejection == 0) {return kTRUE;} else {return kFALSE;} }

  void SetKeepSignalMC() {fKeepSignalMC = kTRUE; return;}


protected:


  Bool_t IsSignalMC(AliAODRecoDecay* d, AliAODEvent* aod, Int_t pdg) const;
  Bool_t RecomputePrimaryVertex(AliAODEvent* event) const;

  // cuts on the event
  Int_t fMinVtxType; // 0: not cut; 1: SPDZ; 2: SPD3D; 3: Tracks
  Int_t fMinVtxContr;   // minimum vertex contributors
  Float_t fMaxVtxRedChi2; // maximum chi2/ndf
  Float_t fMaxVtxZ; // maximum |z| of primary vertex
  Int_t fMinSPDMultiplicity; // SPD multiplicity
  ULong64_t fTriggerMask; // trigger mask
  Bool_t fUseOnlyOneTrigger; // flag to select one trigger only
  TString  fTriggerClass[2]; // trigger class
  // quality cuts on the daughter tracks
  AliESDtrackCuts* fTrackCuts; // tracks for daughter tracks (AOD converted to ESD on the flight!)
  // cuts on the candidate
  Bool_t fUseAOD049; // enable AOD049 centrality cleanup
  Int_t fWhyRejection; // used to code the step at which candidate was rejected
  UInt_t fEvRejectionBits; //bit map storing the full info about event rejection
  Bool_t fRemoveDaughtersFromPrimary; // flag to switch on the removal of duaghters from the primary vertex computation
  Bool_t fRecomputePrimVertex; // flag to recompute primary vertex
  Bool_t fUseMCVertex; // use MC primary vertex
  Bool_t fUsePhysicsSelection; // use Physics selection criteria
  Int_t  fOptPileup;      // option for pielup selection
  Int_t  fMinContrPileup; // min. n. of tracklets in pileup vertex
  Float_t fMinDzPileup;   // min deltaz between main and pileup vertices
  Int_t   fUseCentrality; // off =0 (default)
  // 1 = V0
  // 2 = Tracks
  // 3 = Tracklets
  // 4 = SPD clusters outer
  Float_t fMinCentrality; // minimum centrality for selected events
  Float_t fMaxCentrality; // maximum centrality for selected events
  Bool_t  fFixRefs;       // fix the daughter track references
  Double_t fMinPtJet; // minimum pt of the jet
  Double_t fMaxPtJet; // maximum pt of the jet
  Double_t fMaxEtaJet; // max eta jet to have the jet fully included in eta acceptance
  Double_t fJetRadius; //jet cone radius

  Bool_t  fKeepSignalMC; // IsSelected returns always kTRUE for MC signal

  Bool_t fApplySPDDeadPbPb2011;  // flag to apply SPD dead module map of PbPb2011
  Bool_t fApplySPDMisalignedPP2012; // flag to apply cut on tracks crossing SPD misaligned modules for PP2012 data
  Bool_t fRemoveTrackletOutliers; // flag to apply cut on tracklets vs. centrality for 2011 data
  Int_t fCutOnzVertexSPD; // cut on zSPD vertex to remove outliers in centrality vs. tracklets (0=no cut, 1= cut at 12 cm, 2= cut on difference to z of vtx tracks
  Bool_t fKinkReject; // flag to reject kink daughters
  Bool_t fUseTrackSelectionWithFilterBits; // flag to enable/disable the check on filter bits
  Bool_t fUseCentrFlatteningInMC; // flag for enabling/diabling centrality flattening in MC
  TH1F* fHistCentrDistr;   // histogram with reference centrality distribution for centrality distribution flattening

  ClassDef(AliRDHFJetsCuts, 1); // base class for cuts on AOD reconstructed heavy-flavour decays
};

#endif
