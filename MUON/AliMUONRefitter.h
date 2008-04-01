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

class AliMUONRefitter : public TObject
{
public:
  
  AliMUONRefitter();
  virtual ~AliMUONRefitter();
  
  /// connect to the ESD interface containing MUON data to refit
  void Connect(AliMUONESDInterface* esdInterface) {fESDInterface = esdInterface;}
  
  // re-reconstruct all tracks (clusters) in the ESD event
  AliMUONVTrackStore* ReconstructFromDigits();
  AliMUONVTrackStore* ReconstructFromClusters();
  
  // refit a particular track in the ESD event
  AliMUONTrack* RetrackFromDigits(Int_t iTrack);
  AliMUONTrack* RetrackFromClusters(Int_t iTrack);
  
  // re-clusterize a particular cluster in the ESD event
  AliMUONVClusterStore* ReClusterize(Int_t iTrack, Int_t iCluster);
  AliMUONVClusterStore* ReClusterize(UInt_t clusterId);
  
  
protected:
  
  AliMUONRefitter (const AliMUONRefitter&); ///< copy constructor
  AliMUONRefitter& operator=(const AliMUONRefitter&); ///< assignment operator
  
  
private:
  
  void                   CreateGeometryTransformer();
  void                   CreateClusterServer(AliMUONGeometryTransformer& transformer);
  
  void AddClusterToTracks(const AliMUONVClusterStore &localClusterStore, AliMUONVTrackStore &trackStore);
  
  
private:
    
  AliMUONGeometryTransformer* fGeometryTransformer; ///< geometry transformer (owner)
  AliMUONVClusterServer*      fClusterServer;       ///< clusterizer (owner)
  AliMUONVTrackReconstructor* fTracker;             ///< tracker (owner)
  AliMUONESDInterface*        fESDInterface;        ///< container of MUON tracks/clusters/digits (not owner)
  
  
  ClassDef(AliMUONRefitter,0)
};

#endif

