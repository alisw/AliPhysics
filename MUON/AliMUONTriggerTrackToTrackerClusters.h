#ifndef ALIMUONTRIGGERTRACKTOTRACKERCLUSTERS_H
#define ALIMUONTRIGGERTRACKTOTRACKERCLUSTERS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONTriggerTrackToTrackerClusters
/// \brief Convertor of trigger track to tracker clusters
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONTriggerTrack;
class AliMUONVClusterStore;
class AliMUONVTriggerTrackStore;
class AliMUONGeometryTransformer;

class AliMUONTriggerTrackToTrackerClusters : public TObject
{
public:
  AliMUONTriggerTrackToTrackerClusters(const AliMUONGeometryTransformer& transformer, AliMUONVTriggerTrackStore* trackStore);
  virtual ~AliMUONTriggerTrackToTrackerClusters();

  Int_t GenerateClusters(Int_t iChamber, AliMUONVClusterStore& clusterStore) const;
  
  Int_t GenerateClusters(Int_t iChamber,
                        const AliMUONTriggerTrack& track,
                        AliMUONVClusterStore& clusterStore) const;
    
  Int_t DetElemId(Int_t chamber, Double_t x, Double_t y,
                  Double_t ex, Double_t ey, Double_t& z) const;
    
private:
    /// not defined
    AliMUONTriggerTrackToTrackerClusters(const AliMUONTriggerTrackToTrackerClusters& rhs);
  /// not defined
  AliMUONTriggerTrackToTrackerClusters& operator=(const AliMUONTriggerTrackToTrackerClusters& rhs);
  
private:
    const AliMUONGeometryTransformer& fkTransformer; ///< to go from local to global
    AliMUONVTriggerTrackStore* fTriggerTrackStore; ///< not owner
  
  ClassDef(AliMUONTriggerTrackToTrackerClusters,1) // Convertor of trigger tracks to tracker clusters
};

#endif
