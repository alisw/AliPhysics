#ifndef ALIMUONREFITTER_H
#define ALIMUONREFITTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONRefitter
/// \brief class to refit the ESD clusters/tracks
/// 
//  Author Philippe Pillot

#include <TObject.h>

class AliMUONGeometryTransformer;
class AliMUONVClusterFinder;
class AliMUONVClusterServer;
class AliMUONVTrackReconstructor;
class AliMUONESDInterface;
class AliMUONVClusterStore;
class AliMUONVTrackStore;
class AliMUONTrack;
class AliMUONRecoParam;

class AliMUONRefitter : public TObject
{
public:
  
  AliMUONRefitter(const AliMUONRecoParam* recoParam);
  virtual ~AliMUONRefitter();
  
  /// connect to the ESD interface containing MUON data to refit
  void Connect(const AliMUONESDInterface* esdInterface) {fkESDInterface = esdInterface;}
  
  // re-reconstruct all tracks (clusters) in the ESD event
  AliMUONVTrackStore* ReconstructFromDigits();
  AliMUONVTrackStore* ReconstructFromClusters();
  
  // refit a particular track in the ESD event
  AliMUONTrack* RetrackFromDigits(UInt_t trackId);
  AliMUONTrack* RetrackFromClusters(UInt_t trackId);
  
  // re-clusterize a particular cluster in the ESD event
  AliMUONVClusterStore* ReClusterize(UInt_t trackId, UInt_t clusterId);
  AliMUONVClusterStore* ReClusterize(UInt_t clusterId);
  
  
protected:
  
  AliMUONRefitter (const AliMUONRefitter&); ///< copy constructor
  AliMUONRefitter& operator=(const AliMUONRefitter&); ///< assignment operator
  
  
private:
  
  void CreateGeometryTransformer();
  void CreateClusterServer(AliMUONGeometryTransformer& transformer);
  
  AliMUONTrack* RetrackFromDigits(const AliMUONTrack& track);
  
  void AddClusterToTracks(const AliMUONVClusterStore &localClusterStore, AliMUONVTrackStore &trackStore);
  
private:
    
  const AliMUONRecoParam*     fkRecoParam;          ///< pointer to reco param (not owner)
  const AliMUONESDInterface*  fkESDInterface;       ///< container of MUON tracks/clusters/digits (not owner)
  AliMUONGeometryTransformer* fGeometryTransformer; ///< geometry transformer (owner)
  AliMUONVClusterServer*      fClusterServer;       ///< clusterizer (owner)
  AliMUONVTrackReconstructor* fTracker;             ///< tracker (owner)
  
  
  ClassDef(AliMUONRefitter,0)
};

#endif

