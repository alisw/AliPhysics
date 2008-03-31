#ifndef ALIMUONLEGACYCLUSTERSERVER_H
#define ALIMUONLEGACYCLUSTERSERVER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONLegacyClusterServer
/// \brief Cluster server that always clusterize everything.
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVCLUSTERSERVER_H
#  include "AliMUONVClusterServer.h"
#endif

class AliMUONTriggerTrackToTrackerClusters;
class AliMUONVClusterStore;
class AliMUONGeometryTransformer;

class AliMUONLegacyClusterServer : public AliMUONVClusterServer
{
public:
  AliMUONLegacyClusterServer(const AliMUONGeometryTransformer& transformer, AliMUONVClusterStore* store=0x0);

  virtual ~AliMUONLegacyClusterServer();
  
  virtual Int_t Clusterize(Int_t chamberId, 
                           AliMUONVClusterStore& clusterStore,
                           const AliMpArea& area);
  
  virtual void UseDigits(TIter& next);

  /// Use trigger tracks. Return kFALSE if not used.
  virtual Bool_t UseTriggerTrackStore(AliMUONVTriggerTrackStore* trackStore);

private:
    /// not defined
    AliMUONLegacyClusterServer(const AliMUONLegacyClusterServer& rhs);
  /// not defined
  AliMUONLegacyClusterServer& operator=(const AliMUONLegacyClusterServer& rhs);

  const AliMUONGeometryTransformer& fTransformer; //!< geometry convertor
    AliMUONVClusterStore* fClusterStore; //!< cluster store 
    AliMUONVTriggerTrackStore* fTriggerTrackStore; //!< trigger track store
    AliMUONTriggerTrackToTrackerClusters* fBypass; //!< bypass 
    
  ClassDef(AliMUONLegacyClusterServer,1) // Implementation of AliMUONVClusterServer
};

#endif
