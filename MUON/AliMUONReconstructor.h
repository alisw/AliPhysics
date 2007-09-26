#ifndef ALIMUONRECONSTRUCTOR_H
#define ALIMUONRECONSTRUCTOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONReconstructor
/// \brief Implementation of AliReconstructor for MUON (both tracker and trigger)
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIRECONSTRUCTOR_H
#  include "AliReconstructor.h"
#endif

class AliMUONDigitMaker;
class AliMUONVDigitStore;
class AliMUONVTriggerStore;

class AliMUONGeometryTransformer;

class AliMUONTriggerCrateStore;
class AliMUONTriggerCircuit;
class TClonesArray;
class AliMUONVTriggerStore;

class AliMUONDigitCalibrator;
class AliMUONCalibrationData;

class AliMUONClusterReconstructor;
class AliMUONVClusterStore;

class AliMUONTracker;
class AliMUONVTrackStore;

class AliMUONTriggerChamberEff;

class AliMUONReconstructor : public AliReconstructor
{
public:
  AliMUONReconstructor();
  virtual ~AliMUONReconstructor();
  
  virtual Bool_t HasDigitConversion() const;

  virtual void ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const;
  
  virtual void Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const;
  
  virtual void Reconstruct(TTree* digitsTree, TTree* clustersTree) const;
  
  virtual AliTracker* CreateTracker() const;
    
private:
  /// Not implemented
  AliMUONReconstructor(const AliMUONReconstructor&);
  /// Not implemented
  AliMUONReconstructor& operator=(const AliMUONReconstructor&);
  
  void ConvertDigits(AliRawReader* rawReader, 
                     AliMUONVDigitStore* digitStore,
                     AliMUONVTriggerStore* triggerStore) const;
  void Calibrate(AliMUONVDigitStore& digitStore) const;
  void Clusterize(const AliMUONVDigitStore& digitStore, AliMUONVClusterStore& clusterStore) const;
  AliMUONTriggerCrateStore* CrateManager() const;
  void CreateCalibrator() const;
  void CreateDigitMaker() const;
  void CreateTriggerCircuit() const;
  void CreateClusterReconstructor() const;
  void CreateTriggerChamberEff() const;
  void FillTreeR(AliMUONVTriggerStore* triggerStore,
                 AliMUONVClusterStore* clusterStore,
                 TTree& clustersTree) const;
  
  AliMUONVDigitStore* DigitStore() const;
  AliMUONVClusterStore* ClusterStore() const;
  AliMUONVTriggerStore* TriggerStore() const;

private:
  mutable AliMUONTriggerCrateStore* fCrateManager; //!< Trigger Crate manager
  mutable AliMUONDigitMaker* fDigitMaker; //!< Raw to Digits converter
  AliMUONGeometryTransformer* fTransformer; //!< Geometry transformer (local<->global)
  mutable AliMUONVDigitStore* fDigitStore; //!< Digit container
  mutable AliMUONTriggerCircuit* fTriggerCircuit; //!< Trigger Circuit
  mutable AliMUONCalibrationData* fCalibrationData; //!< Calibration data
  mutable AliMUONDigitCalibrator* fDigitCalibrator; //!<  Digit to calibrate digit converter
  mutable AliMUONClusterReconstructor* fClusterReconstructor; //!<  Clusterizer
  mutable AliMUONVClusterStore* fClusterStore; //!< Cluster container
  mutable AliMUONVTriggerStore* fTriggerStore; //!< Trigger container
  mutable AliMUONVTrackStore* fTrackStore; //!< Track container
  mutable AliMUONTriggerChamberEff* fTrigChamberEff; //!< pointer to trigger chamber efficiency class
  
  ClassDef(AliMUONReconstructor,3) // Implementation of AliReconstructor
};

#endif
