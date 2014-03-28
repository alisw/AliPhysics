#ifndef AliAnalysisTaskCFTree_h
#define AliAnalysisTaskCFTree_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis task to produce trees of lightweight events
// evgeny.kryshen@cern.ch

#include "AliAnalysisTask.h"
class AliInputEventHandler;
class TList;
class TH1I;
class TTree;
class TClonesArray;
class AliAnalysisFilter;
class AliPIDResponse;
class AliVParticle;
class AliCFParticle;

class AliAnalysisTaskCFTree : public AliAnalysisTask {
 public:
  AliAnalysisTaskCFTree(const char* name="AliAnalysisTaskCFTree");
  virtual ~AliAnalysisTaskCFTree(){};
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *option);

  void SetTrackFilter(AliAnalysisFilter* filter) { fTrackFilter = filter; }
  void SetHybridConstrainedMask(UInt_t mask)  { fHybridConstrainedMask = mask; }
  void SetTPConlyConstrainedMask(UInt_t mask) { fTPConlyConstrainedMask = mask; }
  void SetDebug(Int_t val)              { fDebug = val; }
  void SetMode(Int_t val)               { fMode = val; }
  void SetAOD(Bool_t val)               { fIsAOD = val; }
  // Event cut setters
  void SetEventSelectionBit(UInt_t val) { fSelectBit = val; }
  void SetZVertex(Float_t val)          { fZVertexCut = val; }
  void SetTracksInVertex(Int_t val)     { fnTracksVertex = val; }
  void SetCentralityMethod(char* val)   { fCentralityMethod = val; }
  // Track cut setters
  void SetTrackFilterBit(UInt_t val)    { fTrackFilterBit = val; }
  void SetTrackEtaCut(Float_t val)      { fTrackEtaCut = val; }
  void SetPtMin(Float_t val)            { fPtMin = val; } 
  void SetSharedClusterCut(Float_t val) { fSharedClusterCut = val;  }
  void SetCrossedRowsCut(Int_t val)     { fCrossedRowsCut = val; }
  void SetFoundFractionCut(Float_t val) { fFoundFractionCut = val;  }
  // Switchers for additional data to be stored
  void SetStoreTracklets(Bool_t val=kTRUE) { fStoreTracklets = val; }
  void SetStoreMuons(Bool_t val=kTRUE)     { fStoreMuons     = val; }

 protected:
  AliAnalysisTaskCFTree(const  AliAnalysisTaskCFTree &task);
  AliAnalysisTaskCFTree& operator=(const  AliAnalysisTaskCFTree &task);

  UInt_t GetFilterMap(AliVParticle* part);
  AliCFParticle* AddTrack(AliVParticle* track, UInt_t mask, UInt_t flag=0);

  Int_t fDebug;               // debug level
  Int_t fMode;                // Analysis mode: 0 - data, 1 - mc
  Bool_t fIsAOD;              // kTRUE - AOD
  AliInputEventHandler* fInputHandler; // AOD or ESD input handler 
  AliInputEventHandler* fMcHandler;    // MC input handler (ESD)
  AliAnalysisFilter* fTrackFilter;     // track filter used in ESD analysis
  UInt_t fHybridConstrainedMask;       // Filter mask for hybrid constrained tracks (ESD analysis)
  UInt_t fTPConlyConstrainedMask;      // Filter mask for TPConly constrained tracks (ESD analysis)
  AliPIDResponse* fPIDResponse;        //!
  TList* fListOfHistos;       //! list of output histograms
  TH1I*  fEventStatistics;    //! cut-by-cut counter of events
  TTree* fTree;               //! output tree
  // Tree variables
  TClonesArray* fParticles;   //! tree var: selected AliCFParticles
  TClonesArray* fTracklets;   //! tree var: selected tracklets (stored if fStoreTracklets=kTRUE)
  TClonesArray* fMuons;       //! tree var: selected muons (stored if fStoreMuons=kTRUE)
  Float_t fField;             //  tree var: magnetic field value
  Float_t fCentrality;        //  tree var: centrality
  Float_t fZvtx;              //  tree var: z-vertex
  Int_t fRunNumber;           //  tree var: run number
  UInt_t fPeriod;             //  tree var: period
  UInt_t fOrbit;              //  tree var: orbit
  UShort_t fBc;               //  tree var: bunch crossing
  UInt_t fSelectMask;         //  tree var: physics selection mask
  // Event cuts
  UInt_t fSelectBit;          // event selection bit
  Float_t fZVertexCut;        // Z-vertex cut
  Int_t fnTracksVertex;       // min number of vertex contributors
  TString  fCentralityMethod; // method to determine centrality
  // Track cuts
  UInt_t fTrackFilterBit;     // track filter bits to be stored
  Float_t fTrackEtaCut;       // maximum abs(eta) cut
  Float_t fPtMin;             // minimum pt cut
  Float_t fSharedClusterCut;  // cut on shared clusters
  Int_t   fCrossedRowsCut;    // cut on crossed rows
  Float_t fFoundFractionCut;  // cut on crossed rows/findable clusters
  //
  Bool_t fStoreTracklets;     // if kTRUE - SPD tracklets will be stored as AliCFParticles
  Bool_t fStoreMuons;         // if kTRUE - muon tracks will be stored as AliCFParticles
  ClassDef(AliAnalysisTaskCFTree,2);
};
#endif

