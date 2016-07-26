#ifndef ALIAODUPCReplicator_H
#define ALIAODUPCReplicator_H

/* Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$


#ifndef ALIDAODBRANCHREPLICATOR_H
#  include "AliAODBranchReplicator.h"
#endif

//
// Implementation of a branch replicator 
// to produce slim upc aods.
//
// Based on MUON replicator
//
// Author: Michal Broz

class TClonesArray;
class AliAODVZERO;
class AliAODHeader;
class AliAODTracklets;
class AliAODZDC;
class AliAODAD;
class AliAODTrack;

class AliAODUPCReplicator : public AliAODBranchReplicator
{
public:
  
  AliAODUPCReplicator(const char* name="AliAODUPCReplicator", 
                       const char* title="Branch Replicator for muon related branches",
                       Bool_t replicateHeader=kFALSE,
                       Bool_t replicateTracklets=kFALSE);
  virtual ~AliAODUPCReplicator();
  
  virtual TList* GetList() const;
  
  virtual void ReplicateAndFilter(const AliAODEvent& source);
private:

  Bool_t IsGoodGlobalTrack(const AliAODTrack *trk); 
  Bool_t IsGoodITSsaTrack(const AliAODTrack *trk); 

  mutable TClonesArray* fTracks; //! internal array of muon tracks
  mutable TClonesArray* fVertices; //! internal array of vertices
  mutable AliAODVZERO* fVZERO; //! internal vzero object
  mutable AliAODAD* fAD; //!internal ad object
  mutable AliAODHeader* fHeader; //! internal header object
  mutable AliAODTracklets* fTracklets; //! internal tracklets object
  mutable AliAODZDC* fZDC; //! internal zdc object
  mutable TList* fList; //! internal list of managed objects (fVertices and fTracks)

  Bool_t fReplicateHeader; // whether or not the replicate the AOD Header
  Bool_t fReplicateTracklets; // whether or not the replicate the AOD Tracklets
  
private:
  AliAODUPCReplicator(const AliAODUPCReplicator&);
  AliAODUPCReplicator& operator=(const AliAODUPCReplicator&);
  
  ClassDef(AliAODUPCReplicator,1) // Branch replicator for ESD to muon AOD.
};

#endif
