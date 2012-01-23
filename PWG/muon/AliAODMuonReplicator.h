#ifndef ALIAODMUONREPLICATOR_H
#define ALIAODMUONREPLICATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$


#ifndef ALIDAODBRANCHREPLICATOR_H
#  include "AliAODBranchReplicator.h"
#endif
#ifndef ROOT_TExMap
#  include "TExMap.h"
#endif

//
// Implementation of a branch replicator 
// to produce slim muon and dimuon aods.
//
// Author: L. Aphecetche (Subatech)

class AliAnalysisCuts;
class TClonesArray;
class AliAODMCHeader;
class AliAODVZERO;

class AliAODMuonReplicator : public AliAODBranchReplicator
{
public:
  
  AliAODMuonReplicator(const char* name="AliAODMuonReplicator", 
                       const char* title="Branch Replicator for muon related branches",
                       AliAnalysisCuts* trackCut=0x0,
                       AliAnalysisCuts* vertexCut=0x0,
                       Int_t mcMode=0);
  virtual ~AliAODMuonReplicator();
  
  virtual TList* GetList() const;
  
  virtual void ReplicateAndFilter(const AliAODEvent& source);
  
private:
  void FilterMC(const AliAODEvent& source);
  void SelectParticle(Int_t i);
  Bool_t IsParticleSelected(Int_t i);
  void CreateLabelMap(const AliAODEvent& source);
  Int_t GetNewLabel(Int_t i);
  
private:
  AliAnalysisCuts* fTrackCut; // decides which tracks to keep
  mutable TClonesArray* fTracks; //! internal array of muon tracks
  AliAnalysisCuts* fVertexCut; // decides which vertices to keep
  mutable TClonesArray* fVertices; //! internal array of vertices
  mutable TClonesArray* fDimuons; //! internal array of dimuons
  mutable AliAODVZERO* fVZERO; //! internal vzero object
  mutable TList* fList; //! internal list of managed objects (fVertices and fTracks)
  
  mutable TClonesArray* fMCParticles; //! internal array of MC particles
  mutable AliAODMCHeader* fMCHeader; //! internal array of MC header
  Int_t fMCMode; // MC filtering switch (0=none=no mc information,1=normal=simple copy,>=2=aggressive=filter out)
  TExMap fLabelMap; //! for MC label remapping (in case of aggressive filtering)
  TExMap fParticleSelected; //! List of selected MC particles

private:
  AliAODMuonReplicator(const AliAODMuonReplicator&);
  AliAODMuonReplicator& operator=(const AliAODMuonReplicator&);
  
  ClassDef(AliAODMuonReplicator,4) // Branch replicator for ESD to muon AOD.
};

#endif
