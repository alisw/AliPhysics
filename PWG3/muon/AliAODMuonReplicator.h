#ifndef ALIAODMUONREPLICATOR_H
#define ALIAODMUONREPLICATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$


#ifndef ALIDAODBRANCHREPLICATOR_H
#  include "AliAODBranchReplicator.h"
#endif

//
// Implementation of a branch replicator 
// to produce slim muon and dimuon aods.
//
// Author: L. Aphecetche (Subatech)

class AliAnalysisCuts;
class TClonesArray;

class AliAODMuonReplicator : public AliAODBranchReplicator
{
public:
  AliAODMuonReplicator(const char* name="AliAODMuonReplicator", 
                       const char* title="Branch Replicator for muon related branches",
                       AliAnalysisCuts* trackCut=0x0,
                       AliAnalysisCuts* vertexCut=0x0);
  virtual ~AliAODMuonReplicator();
  
  virtual TList* GetList() const;
  
  virtual void ReplicateAndFilter(const AliAODEvent& source);
  
private:
  AliAnalysisCuts* fTrackCut; // decides which tracks to keep
  mutable TClonesArray* fTracks; //! internal array of muon tracks
  AliAnalysisCuts* fVertexCut; // decides which vertices to keep
  mutable TClonesArray* fVertices; //! internal array of vertices
  mutable TClonesArray* fDimuons; //! internal array of dimuons
  mutable TList* fList; //! internal list of managed objects (fVertices and fTracks)
  
private:
  AliAODMuonReplicator(const AliAODMuonReplicator&);
  AliAODMuonReplicator& operator=(const AliAODMuonReplicator&);
  
  ClassDef(AliAODMuonReplicator,2) // Branch replicator for ESD to muon AOD.
};

#endif
