#ifndef ALINANOAODREPLICATO_H
#define ALINANOAODREPLICATO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id: AliAODMuonReplicator.h 56492 2012-05-15 18:42:47Z pcrochet $


#ifndef ALIDAODBRANCHREPLICATOR_H
#  include "AliAODBranchReplicator.h"
#endif
#ifndef ROOT_TExMap
#  include "TExMap.h"
#endif

#include <iostream>
#include <list>
//
// Implementation of a branch replicator 
// to produce nano AOD.
//
// Author: Michele Floris, michele.floris@cern.ch

class AliAnalysisCuts;
class TClonesArray;
class AliAODMCHeader;
class AliAODVZERO;
class AliAODTZERO;
class AliPIDResponse;
class AliESDv0;
class TArrayI;
class AliAODv0;  
class TRefArray;
class AliAODRecoDecay;
class AliAODRecoDecayHF;
class AliAODRecoDecayHF2Prong;
class AliVertexerTracks;

class AliESDVertex;
class AliESDtrack;
class AliVEvent;
class AliAODVertex;
class AliVertexerTracks;
class AliESDv0; 
class AliAODv0; 
class AliAODConversionPhoton;
class AliAODHeader;
class AliNanoAODHeader;
class AliAnalysisTaskSE;
class AliNanoAODTrack;
class AliAODTrack;
class AliNanoAODCustomSetter;
class AliAODZDC;

class AliNanoAODReplicator : public AliAODBranchReplicator
{
 public:
  
  AliNanoAODReplicator();
  AliNanoAODReplicator(const char* name, const char* title);
  
  virtual ~AliNanoAODReplicator();
  
  virtual TList* GetList() const ; // FIXME: This is declared const in the interface
  
  virtual void ReplicateAndFilter(const AliAODEvent& source);	

  virtual void Terminate();

  const char * GetVarListTrack() { return fVarList; }
  void  SetVarListTrack (const char * var) { fVarList = var;}
  const char * GetVarListHeader() { return fVarListHeader; }
  void  SetVarListHeader (const char * var) { fVarListHeader = var;}
  
  void SetTrackCuts(AliAnalysisCuts* cuts) { fTrackCuts = cuts; }
  void SetV0Cuts(AliAnalysisCuts* cuts) { fV0Cuts = cuts; }
  void SetCascadeCuts(AliAnalysisCuts* cuts) { fCascadeCuts = cuts; }
  void SetConversionPhotonCuts(AliAnalysisCuts* cuts) { fConversionPhotonCuts = cuts; }
  void SetMCParticleCuts(AliAnalysisCuts* cuts) { fMCParticleCuts = cuts; }

  void AddCustomSetter(AliNanoAODCustomSetter * var) { fCustomSetters.push_back(var);  }
    
  void SetSaveVzero(Bool_t b)  { fSaveVzero = b; }
  void SetSaveZDC(Bool_t b)    { fSaveZDC = b; }
  void SetSaveV0s(Bool_t b)    { fSaveV0s = b; }
  void SetSaveCascades(Bool_t b) { fSaveCascades = b; }
  void SetSaveConversionPhotons(Bool_t b) { fSaveConversionPhotons = b; }
  void SetPhotonDeltaBranchName(TString name) {
    fPhotonFromDeltas = true;
    fDeltaAODBranchName = name;
  }
  
  void SetMCMode(Int_t mode)  { fMCMode = mode; }
  
  void SetInputArrayName(TString name) {fInputArrayName=name;}
  void SetOutputArrayName(TString name) {fOutputArrayName=name;}

  void SetVarListHeaderTC(TString var) {fVarListHeader_fTC=var;}
    
 private:

  void SelectParticle(Int_t i);
  Bool_t IsParticleSelected(Int_t i);
  void CreateLabelMap(const AliAODEvent& source);
  Int_t GetNewLabel(Int_t i);
  void RelabelAODPhotonCandidates(AliAODConversionPhoton *PhotonCandidate);
  void FilterMC(const AliAODEvent& source);
  AliAODVertex* CloneAndStoreVertex(AliAODVertex* toClone);
 
  AliAnalysisCuts* fTrackCuts; // decides which tracks to keep
  AliAnalysisCuts* fV0Cuts;    // decides which V0s to keep
  AliAnalysisCuts* fCascadeCuts; // decides which cascades to keep
  AliAnalysisCuts* fConversionPhotonCuts; // decides which conversion photons to keep
  AliAnalysisCuts* fMCParticleCuts;  // decides which particles from the MC stack to keep
                                                      // Here we need the actual object because we retrieve the MC
                                                      // matching of the V0s from here
  
  mutable TClonesArray* fTracks; //! internal array of arrays of NanoAOD tracks
  mutable AliNanoAODHeader* fHeader; //! internal array of headers
 
  mutable TClonesArray* fVertices; //! internal array of vertices
 
  mutable TList* fList; //! internal list of managed objects (fVertices and fTracks)
  
  mutable TClonesArray* fMCParticles; //! internal array of MC particles
  mutable AliAODMCHeader* fMCHeader; //! internal array of MC header
  Int_t fMCMode; // MC filtering switch (0=none=no mc information,1=normal=simple copy,>=2=aggressive=filter out : keep only particles leading to tracks and trheir relatives + all charged primaries)

  TExMap fLabelMap; //! for MC label remapping (in case of aggressive filtering)
  TExMap fParticleSelected; //! List of selected MC particles
			
  TString fVarList; // list of variables to be filterered
  TString fVarListHeader; // list of variables to be filtered (header)
  TString fVarListHeader_fTC;// list of fired Trigger Classes which are used in the NanoAOD generation 

  std::list<AliNanoAODCustomSetter*> fCustomSetters;       // list of custom setters
    
  mutable AliAODVZERO* fVzero; //! internal array of AliAODVZEROs
  mutable AliAODZDC* fAodZDC; //! internal array of AliAODZDCs
  mutable TClonesArray* fV0s;    //! internal array of AliAODv0
  mutable TClonesArray* fCascades;    //! internal array of AliAODcascade
  mutable TClonesArray* fConversionPhotons;    //! internal array of AliAODConversionPhoton
    
  Bool_t fSaveZDC;    // if kTRUE AliAODZDC will be saved in AliAODEvent
  Bool_t fSaveVzero;  // if kTRUE AliAODVZERO will be saved in AliAODEvent
  Bool_t fSaveV0s;    // if kTRUE AliAODv0 will be saved in AliAODEvent
  Bool_t fSaveCascades; // if kTRUE AliAODcascade will be saved in AliAODEvent
  Bool_t fSaveConversionPhotons; // If kTRUE gamme conversions are stored (needs delta AOD)
  Bool_t fPhotonFromDeltas; // If kTRUE gamma conversions will be directly taken from the Delta AOD
  TString fDeltaAODBranchName; // Name of the photon branch in the Delta AOD

  TString fInputArrayName; // name of array if tracks are stored in a TObjectArray
  TString fOutputArrayName; // name of the output array, where the NanoAODTracks are stored
  
  std::map<AliAODVertex*, std::vector<TObject*> > fKeepDaughters; //! Tracks needed as references to V0s and cascades
  std::map<AliAODVertex*, AliAODVertex*> fClonedVertices; //! avoid that vertices are stored several times

  AliNanoAODReplicator(const AliNanoAODReplicator&);
  AliNanoAODReplicator& operator=(const AliNanoAODReplicator&);

  ClassDef(AliNanoAODReplicator, 7) // Branch replicator for ESD to muon AOD.
};

#endif
