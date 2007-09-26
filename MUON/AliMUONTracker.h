#ifndef ALIMUONTRACKER_H
#define ALIMUONTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/// \ingroup rec
/// \class AliMUONTracker
/// \brief MUON base Tracker
///
//  Author: Christian Finck, SUBATECH Nantes

#include "AliTracker.h"

class AliCluster;
class AliESDEvent;
class AliMUONDigitMaker;
class AliMUONGeometryTransformer;
class AliMUONTrackHitPattern;
class AliMUONTriggerChamberEff;
class AliMUONTriggerCircuit;
class AliMUONVClusterStore;
class AliMUONVTrackReconstructor;
class AliMUONVTrackStore;
class AliMUONVTriggerStore;

class AliMUONTracker : public AliTracker
{
 public:

  AliMUONTracker(const AliMUONDigitMaker* digitMaker=0,
                 const AliMUONGeometryTransformer* transformer=0,
                 const AliMUONTriggerCircuit* triggerCircuit=0,
		 AliMUONTriggerChamberEff* chamberEff=0);
  virtual ~AliMUONTracker();
  
  /// Main entry point.
  virtual Int_t Clusters2Tracks(AliESDEvent* esd);

  /// 
  virtual Int_t LoadClusters(TTree* clustersTree);

  ///
  virtual void  UnloadClusters();

  /// Dummy implementation
  virtual Int_t PropagateBack(AliESDEvent* /*event*/) {return 0;}
  /// Dummy implementation
  virtual Int_t RefitInward(AliESDEvent* /*event*/) {return 0;}
  /// Dummy implementation
  virtual AliCluster *GetCluster(Int_t /*index*/) const {return 0;}
  /// Set option
  void SetOption(Option_t* opt);

private:
  /// Not implemented
  AliMUONTracker(const AliMUONTracker& rhs);
  /// Not implemented
  AliMUONTracker& operator=(const AliMUONTracker& rhs);
    
  Int_t Clusters2Tracks(TTree& tracksTree, AliESDEvent* esd);

  void FillESD(AliMUONVTrackStore& trackStore, AliESDEvent* esd) const;

private:
  const AliMUONDigitMaker* fDigitMaker; //!< digit maker (not owner)
  const AliMUONGeometryTransformer* fTransformer; //!< geometry transformer (not owner)
  const AliMUONTriggerCircuit* fTriggerCircuit;                //!< trigger circuit (not owner)
  AliMUONTriggerChamberEff* fTrigChamberEff; //!< trigger efficiency (not owner)
  AliMUONTrackHitPattern* fTrackHitPatternMaker; //!< trigger hit pattern maker
  AliMUONVTrackReconstructor* fTrackReco;       //!< track reconstructor
  AliMUONVClusterStore* fClusterStore; //!< cluster container
  AliMUONVTriggerStore* fTriggerStore; //!< trigger information
  
  ClassDef(AliMUONTracker,0)  //tracker base class for MUON
};
#endif
