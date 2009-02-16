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
  AliMUONLegacyClusterServer(const AliMUONGeometryTransformer& transformer, 
                             AliMUONVClusterStore* store=0x0,
                             Bool_t bypassSt4=kFALSE,
                             Bool_t bypassSt5=kFALSE);

  virtual ~AliMUONLegacyClusterServer();
  
  virtual Int_t Clusterize(Int_t chamberId, 
                           AliMUONVClusterStore& clusterStore,
                           const AliMpArea& area,
                           const AliMUONRecoParam* recoParam = 0x0);
  
  virtual void UseDigits(TIter& next);

  /// Use trigger tracks. Return kFALSE if not used.
  virtual Bool_t UseTriggerTrackStore(AliMUONVTriggerTrackStore* trackStore);

private:
    /// not defined
    AliMUONLegacyClusterServer(const AliMUONLegacyClusterServer& rhs);
  /// not defined
  AliMUONLegacyClusterServer& operator=(const AliMUONLegacyClusterServer& rhs);

  const AliMUONGeometryTransformer& fkTransformer; //!< geometry convertor
	AliMUONVClusterStore* fClusterStore; //!< cluster store 
	AliMUONVTriggerTrackStore* fTriggerTrackStore; //!< trigger track store
	AliMUONTriggerTrackToTrackerClusters* fBypass; //!< bypass 
	Bool_t fBypassSt4; //!< whether we should bypass station 4
	Bool_t fBypassSt5; //!< whether we should bypass station 5

  ClassDef(AliMUONLegacyClusterServer,2) // Implementation of AliMUONVClusterServer
};

#endif
