#ifndef AliAnalysisTaskCFTree_h
#define AliAnalysisTaskCFTree_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis task to produce trees of lightweight events
// evgeny.kryshen@cern.ch

#include "AliAnalysisTaskSE.h"
class TList;
class TH1I;
class TTree;
class TClonesArray;
class AliAnalysisFilter;
class AliVTrack;
class AliCFParticle;
class AliAnalysisUtils;

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
  void SetEventSelectionBit(UInt_t val) { fSelectBit = val; }
  void SetZVertex(Float_t val)          { fZVertexCut = val; }
  // Track cut setters
  void SetTrackFilterBit(UInt_t val)    { fTrackFilterBit = val; }
  void SetTrackEtaCut(Float_t val)      { fTrackEtaCut = val; }
  void SetPtMin(Float_t val)            { fPtMin = val; } 
  void SetSharedClusterCut(Float_t val) { fSharedClusterCut = val;  }
  void SetCrossedRowsCut(Int_t val)     { fCrossedRowsCut = val; }
  void SetFoundFractionCut(Float_t val) { fFoundFractionCut = val;  }
  void SetDphiCut(Float_t val)          { fDphiCut = val; }
  // Switchers for additional data to be stored
  void SetStoreTracks(Bool_t val=kTRUE)      { fStoreTracks      = val; }
  void SetStoreTracklets(Bool_t val=kTRUE)   { fStoreTracklets   = val; }
  void SetStoreMuons(Bool_t val=kTRUE)       { fStoreMuons       = val; }
  void SetStoreMcTracks(Bool_t val=kTRUE)    { fStoreMcTracks    = val; }
  void SetStoreMcTracklets(Bool_t val=kTRUE) { fStoreMcTracklets = val; }
  void SetStoreMcMuons(Bool_t val=kTRUE)     { fStoreMcMuons     = val; }
  void SetStorePidInfo(Bool_t val=kTRUE)     { fStorePidInfo     = val; }

 protected:
  AliAnalysisTaskCFTree(const  AliAnalysisTaskCFTree &task);
  AliAnalysisTaskCFTree& operator=(const  AliAnalysisTaskCFTree &task);

  UInt_t GetFilterMap(AliVTrack* part);
  AliCFParticle* AddTrack(AliVTrack* track, UInt_t mask, UInt_t flag=0);

  AliAnalysisFilter* fTrackFilter; // track filter used in ESD analysis
  UInt_t fHybridConstrainedMask;       // Filter mask for hybrid constrained tracks (ESD analysis)
  UInt_t fTPConlyConstrainedMask;      // Filter mask for TPConly constrained tracks (ESD analysis)
  AliAnalysisUtils* fUtils;   //! analysis utils to detect pileup
  TList* fListOfHistos;       //! list of output histograms
  TH1I*  fEventStatistics;    //! cut-by-cut counter of events
  TTree* fTree;               //! output tree
  // Tree variables
  TClonesArray* fTracks;      //! tree var: selected AliCFParticles
  TClonesArray* fTracklets;   //! tree var: selected tracklets (stored if fStoreTracklets=kTRUE)
  TClonesArray* fMuons;       //! tree var: selected muons (stored if fStoreMuons=kTRUE)
  TClonesArray* fMcParticles; //! tree var: MC particles
  Float_t fField;             //  tree var: magnetic field value
  Float_t fCentrality[6];     //  tree var: centrality
  Float_t fVtxZ;              //  tree var: z-vertex
  Bool_t fVtxTPConly;         //  tree var: is vertex TPC only
  UInt_t fVtxContributors;    //  tree var: number of vertex contributors
  UInt_t fPeriod;             //  tree var: period
  UInt_t fOrbit;              //  tree var: orbit
  UShort_t fBc;               //  tree var: bunch crossing
  UInt_t fSelectMask;         //  tree var: physics selection mask
  Bool_t fIsPileupSPD;        //  tree var: is pileup from SPD flag
  Bool_t fIsPileupMV;         //  tree var: is pileup from MV flag

  // Event cuts
  UInt_t fSelectBit;          // event selection bit
  Float_t fZVertexCut;        // Z-vertex cut
  // Track cuts
  UInt_t fTrackFilterBit;     // track filter bits to be stored
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
  Bool_t fStorePidInfo;       // if kTRUE - Store PID info for tracks
  ClassDef(AliAnalysisTaskCFTree,4);
};
#endif

