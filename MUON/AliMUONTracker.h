#ifndef ALIMUONTRACKER_H
#define ALIMUONTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONTracker
/// \brief MUON base Tracker
///
//  Authors: Christian Finck, Laurent Aphecetche, SUBATECH Nantes

#include "AliTracker.h"

class AliCluster;
class AliESDEvent;
class AliMUONGeometryTransformer;
class AliMUONRecoParam;
class AliMUONTrackHitPattern;
class AliMUONTriggerCircuit;
class AliMUONVClusterServer;
class AliMUONVClusterStore;
class AliMUONVDigitStore;
class AliMUONVTrackReconstructor;
class AliMUONVTrackStore;
class AliMUONVTriggerStore;
class AliMUONVTriggerTrackStore;
class AliMUONTriggerUtilities;

class AliMUONTracker : public AliTracker
{
 public:

  AliMUONTracker(const AliMUONRecoParam* recoParam,
                 AliMUONVClusterServer* clusterServer,
                 AliMUONVDigitStore& digitStore,
                 const AliMUONGeometryTransformer* transformer=0,
                 const AliMUONTriggerCircuit* triggerCircuit=0,
                 const AliMUONTriggerUtilities* triggerUtilities=0);
  virtual ~AliMUONTracker();
  
  virtual Int_t Clusters2Tracks(AliESDEvent* esd);

  virtual Int_t LoadClusters(TTree* clustersTree);

  virtual void  UnloadClusters();

  /// Return reco parameters
  const AliMUONRecoParam* GetRecoParam() const { return fkRecoParam; }
  
  /// Dummy implementation
  virtual Int_t PropagateBack(AliESDEvent* /*event*/) {return 0;}
  /// Dummy implementation
  virtual Int_t RefitInward(AliESDEvent* /*event*/) {return 0;}
  /// Dummy implementation
  virtual AliCluster *GetCluster(Int_t /*index*/) const {return 0;}

  static AliMUONVTrackReconstructor* CreateTrackReconstructor(const AliMUONRecoParam* recoParam, AliMUONVClusterServer* clusterServer);
  
private:
  /// Not implemented
  AliMUONTracker(const AliMUONTracker& rhs);
  /// Not implemented
  AliMUONTracker& operator=(const AliMUONTracker& rhs);
    
  AliMUONVClusterStore* ClusterStore() const;

  AliMUONVTriggerTrackStore* TriggerTrackStore() const;
  
  void FillESD(const AliMUONVTrackStore& trackStore, AliESDEvent* esd) const;

  void SetupClusterServer(AliMUONVClusterServer& clusterServer);
  
private:
  const AliMUONGeometryTransformer* fkTransformer; //!< geometry transformer (not owner)
  const AliMUONTriggerCircuit* fkTriggerCircuit; //!< trigger circuit (not owner)
  AliMUONTrackHitPattern* fTrackHitPatternMaker; //!< trigger hit pattern maker
  AliMUONVTrackReconstructor* fTrackReco; //!< track reconstructor
  mutable AliMUONVClusterStore* fClusterStore; //!< cluster container
  AliMUONVTriggerStore* fTriggerStore; //!< trigger information
  AliMUONVClusterServer* fClusterServer; //!< to get clusters
  Bool_t fIsOwnerOfClusterServer; //!< whether we are owner of the cluster server
  const AliMUONVDigitStore& fkDigitStore; //!< digit info to fill in ESD
  mutable AliMUONVClusterStore* fInputClusterStore; //!< cluster container
  mutable AliMUONVTriggerTrackStore* fTriggerTrackStore; //!< trigger track store
  const AliMUONRecoParam* fkRecoParam; //!< pointer to reco param
  
  ClassDef(AliMUONTracker,0)  //tracker base class for MUON
};
#endif
