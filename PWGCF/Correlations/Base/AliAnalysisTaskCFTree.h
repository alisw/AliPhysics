#ifndef AliAnalysisTaskCFTree_h
#define AliAnalysisTaskCFTree_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis task to produce trees of lightweight events
// evgeny.kryshen@cern.ch
// lmilano@cern.ch
// alice.ohlson@cern.ch

#include "AliAnalysisTaskSE.h"
class TList;
class TH1I;
class TH1F;
class TH2F;
class TTree;
class TClonesArray;
class AliAnalysisFilter;
class AliVTrack;
class AliCFParticle;
class AliAnalysisUtils;
class AliMuonTrackCuts;
class TPythia6Decayer;

class AliAnalysisTaskCFTree : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskCFTree(const char* name="AliAnalysisTaskCFTree");
  virtual ~AliAnalysisTaskCFTree(){};
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  void SetTrackFilter(AliAnalysisFilter* filter) { fTrackFilter = filter; }
  void SetHybridConstrainedMask(UInt_t mask)  { fHybridConstrainedMask = mask; }
  void SetTPConlyConstrainedMask(UInt_t mask) { fTPConlyConstrainedMask = mask; }
  // Event cut setters
  void SetClassBit(UInt_t val)          { fClassBit = val; }
  void SetEventSelectionBit(UInt_t val) { fSelectBit = val; }
  void SetZVertex(Float_t val)          { fZVertexCut = val; }
  // Track cut setters
  void SetTrackFilterBit(UInt_t val)    { fTrackFilterBit = val; }
  void SetTrackEtaCut(Float_t val)      { fTrackEtaCut = val; }
  void SetPtMin(Float_t val)            { fPtMin = val; } 
  void SetSharedClusterCut(Float_t val) { fSharedClusterCut = val;  }
  void SetCrossedRowsCut(Int_t val)     { fCrossedRowsCut = val; }
  void SetFoundFractionCut(Float_t val) { fFoundFractionCut = val;  }
  // Tracklet cut setters
  void SetTrackletEtaCut(Float_t val)   { fTrackletEtaCut = val; }
  void SetDphiCut(Float_t val)          { fDphiCut = val; }
  // Switchers for additional data to be stored
  void SetStoreTracks(Bool_t val=kTRUE)      { fStoreTracks      = val; }
  void SetStoreTracklets(Bool_t val=kTRUE)   { fStoreTracklets   = val; }
  void SetStoreMuons(Bool_t val=kTRUE)       { fStoreMuons       = val; }
  void SetStoreMcTracks(Bool_t val=kTRUE)    { fStoreMcTracks    = val; }
  void SetStoreMcTracklets(Bool_t val=kTRUE) { fStoreMcTracklets = val; }
  void SetStoreMcMuons(Bool_t val=kTRUE)     { fStoreMcMuons     = val; }
  void SetStoreTrackInfo(Bool_t val=kTRUE)   { fStoreTrackInfo   = val; }
  void SetStorePidInfo(Bool_t val=kTRUE)     { fStorePidInfo     = val; }
  void SetStoreMuonOrigin(Bool_t val=kTRUE)  { fStoreMuonOrigin  = val; }
  void SetIs13TeV(Bool_t val=kTRUE)          { fIs13TeV          = val; }

  void SetApplyPhysicsSelectionCut(Bool_t val=kTRUE) { fApplyPhysicsSelectionCut = val; }
  void SetStoreOnlyEventsWithMuons(Bool_t val=kTRUE) { fStoreOnlyEventsWithMuons = val; }
  void SetStoreCutBitsInTrackMask(Bool_t val=kTRUE)  { fStoreCutBitsInTrackMask  = val; }
 protected:
  AliAnalysisTaskCFTree(const  AliAnalysisTaskCFTree &task);
  AliAnalysisTaskCFTree& operator=(const  AliAnalysisTaskCFTree &task);

  UInt_t GetFilterMap(AliVTrack* part);
  AliCFParticle* AddTrack(AliVTrack* track, UInt_t mask, UInt_t flag=0);

  AliAnalysisFilter* fTrackFilter; // track filter used in ESD analysis
  UInt_t fHybridConstrainedMask;       // Filter mask for hybrid constrained tracks (ESD analysis)
  UInt_t fTPConlyConstrainedMask;      // Filter mask for TPConly constrained tracks (ESD analysis)
  AliMuonTrackCuts* fMuonTrackCuts; // muon track cuts used to extract pxDCA decision
  AliAnalysisUtils* fUtils;   //! analysis utils to detect pileup
  TList* fListOfHistos;       //! list of output histograms
  TH1I*  fEventStatistics;    //! cut-by-cut counter of events
  TH1D*  fClassStatistics;          //! statistics on trigger classes
  TH1F* fV0chan;              //! V0 channel signals (test for Michele F.)
  TTree* fTree;               //! output tree
  // Tree variables
  TClonesArray* fTracks;      //! tree var: selected AliCFParticles
  TClonesArray* fTracklets;   //! tree var: selected tracklets (stored if fStoreTracklets=kTRUE)
  TClonesArray* fMuons;       //! tree var: selected muons (stored if fStoreMuons=kTRUE)
  TClonesArray* fMcParticles; //! tree var: MC particles
  TClonesArray* fMcMuons;     //! tree var: MC muons
  TClonesArray* fMuonOrigin;  //! tree var: Muon origin
  Bool_t fIs13TeV;            //  flag for 13 TeV running (kTRUE by default)
  UInt_t fClassesFired;       //  tree var: classes fired (bit mask, see cxx for details) 
  Float_t fField;             //  tree var: magnetic field value
  Int_t fCurrentRunNumber;    //  tree var: run number
  Float_t fCentrality[6];     //  tree var: centrality
  Float_t fVtxX;              //  tree var: x-vertex
  Float_t fVtxY;              //  tree var: y-vertex
  Float_t fVtxZ;              //  tree var: z-vertex
  Bool_t fVtxTPConly;         //  tree var: is vertex TPC only
  UInt_t fVtxContributors;    //  tree var: number of vertex contributors
  Float_t fSPDVtxZ;           //  tree var: z-vertex (spd vertex)
  Float_t fSPDVtxZRes;        //  tree var: z-vertex resolution (spd vertex)
  UInt_t fSPDVtxContributors; //  tree var: number of vertex contributors (spd vertex)
  UInt_t fPeriod;             //  tree var: period
  UInt_t fOrbit;              //  tree var: orbit
  UShort_t fBc;               //  tree var: bunch crossing
  UInt_t fSelectMask;         //  tree var: physics selection mask
  UInt_t fIsDiffractive;      //  tree var: diffractive events
  Bool_t fIsPileupSPD;        //  tree var: is pileup from SPD flag
  Bool_t fIsPileupMV;         //  tree var: is pileup from MV flag
  TBits fIR1;                 //  tree var: IR1 contains V0 information (VIR)
  TBits fIR2;                 //  tree var: IR2 contains T0 information
  UInt_t  fonlineSPD;         //  tree var: onlineSPD
  UInt_t  fofflineSPD;        //  tree var: offlineSPD
  Int_t  fV0ADecision;        //  tree var: V0A decision
  Int_t  fV0CDecision;        //  tree var: V0C decision
  Int_t  fIsIncomplete;       //  tree var: incomplete events
  UInt_t fNofITSClusters[6];  //  tree var: number of ITS clusters per layer
  Float_t fMultV0Aeq;         //  tree var: multiplicity in V0A from AliMultSelection
  Float_t fMultV0Ceq;         //  tree var: multiplicity in V0C from AliMultSelection
  Float_t fMultV0Meq;         //  tree var: multiplicity in V0M from AliMultSelection
  Float_t fMultV0Aring[4];    //  tree var: V0A multiplicity in 4 rings
  Float_t fMultV0Cring[4];    //  tree var: V0A multiplicity in 4 rings
  Float_t fMultTKL;           //  tree var: multiplicity of tracklets
  Float_t fMultPercV0Aeq;     //  tree var: V0Aeq multiplicity percentile from AliMultSelection
  Float_t fMultPercV0Ceq;     //  tree var: V0Ceq multiplicity percentile from AliMultSelection
  Float_t fMultPercV0Meq;     //  tree var: V0Meq multiplicity percentile from AliMultSelection
  Float_t fMultPercTKL;       //  tree var: TKL multiplicity percentile from AliMultSelection
  Float_t fMultPercV0S;       //  tree var: V0S multiplicity percentile from AliPPVsMultUtils
  Float_t fMultMeanV0A;       //  tree var: V0A mean multiplicity percentile from AliMultSelection
  Float_t fMultMeanV0C;       //  tree var: V0C mean multiplicity percentile from AliMultSelection
  Float_t fMultMeanV0M;       //  tree var: V0M mean multiplicity percentile from AliMultSelection
  Float_t fMultMeanTKL;       //  tree var: TKL mean multiplicity percentile from AliMultSelection
  Int_t fIsEventSel;          //  tree var: is event selected by AliMultSelection class
  Int_t fNchTPC;              //  tree var: Nch in TPC for |eta|<1 and pT>0.2
  Int_t fNchTPCmc;            //  tree var: Nch in TPC for |eta|<1 and pT>0.2 - generated MC
  Int_t fNchV0Amc;            //  tree var: Nch in the V0A acceptance - generated MC
  Int_t fNchV0Cmc;            //  tree var: Nch in the V0C acceptance - generated MC
  Int_t fNchCL1mc;            //  tree var: Nch in the CL1 acceptance - generated MC
  // Event cuts
  UInt_t fClassBit;           // class selection mask (see cxx for datails)
  UInt_t fSelectBit;          // event selection bit
  Float_t fZVertexCut;        // Z-vertex cut
  // Track cuts
  UInt_t fTrackFilterBit;     // track filter bits to be stored
  Float_t fTrackletEtaCut;       // maximum abs(eta) cut for tracklets
  Float_t fTrackEtaCut;       // maximum abs(eta) cut
  Float_t fPtMin;             // minimum pt cut
  Float_t fSharedClusterCut;  // cut on shared clusters
  Int_t   fCrossedRowsCut;    // cut on crossed rows
  Float_t fFoundFractionCut;  // cut on crossed rows/findable clusters
  Float_t fDphiCut;           // cut on tracklet dphi
  //
  Bool_t fStoreTracks;        // if kTRUE - Barrel tracks will be stored as AliCFParticles
  Bool_t fStoreTracklets;     // if kTRUE - SPD tracklets will be stored as AliCFParticles
  Bool_t fStoreMuons;         // if kTRUE - muon tracks will be stored as AliCFParticles
  Bool_t fStoreMcTracks;      // if kTRUE - mc particles will be stored as AliCFParticles
  Bool_t fStoreMcTracklets;   // if kTRUE - Store Monte-Carlo info for tracklets
  Bool_t fStoreMcMuons;       // if kTRUE - Store Monte-Carlo info for muons
  Bool_t fStoreTrackInfo;     // if kTRUE - Store additional track info on tracks
  Bool_t fStorePidInfo;       // if kTRUE - Store PID info for tracks
  Bool_t fStoreMuonOrigin;    // if kTRUE - Store muon origin in a TClonesArray of TObjStrings
  Bool_t fApplyPhysicsSelectionCut; // skip events not passing fSelectionBit mask
  Bool_t fStoreOnlyEventsWithMuons; // if kTRUE store only events with at least one muon
  Bool_t fStoreCutBitsInTrackMask;  // if kTRUE modify additional bits in track mask
  TClonesArray* fDecayArray;
  TPythia6Decayer* fDecayer;

  ClassDef(AliAnalysisTaskCFTree,7);
};
#endif

