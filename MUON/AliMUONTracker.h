#ifndef ALIMUONTRACKER_H
#define ALIMUONTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/// \ingroup rec
/// \class AliMUONTracker
/// \brief MUON base Tracker

#include "AliTracker.h"
class AliESD;
class AliCluster;
class AliMUONRecData;
class TClonesArray;
class AliMUONVTrackReconstructor;

class AliMUONTracker : public AliTracker
{
 public:

  AliMUONTracker();
  virtual ~AliMUONTracker();
    
  virtual Int_t Clusters2Tracks(AliESD* /*esd*/); 
  
  /// Dummy implementation
  virtual Int_t PropagateBack(AliESD* /*event*/) {return 0;}
  /// Dummy implementation
  virtual Int_t RefitInward(AliESD* /*event*/) {return 0;}
  /// Dummy implementation
  virtual Int_t LoadClusters(TTree* /*tree*/) {return 0;}
  /// Dummy implementation
  virtual void  UnloadClusters() {return;}
  /// Dummy implementation
  virtual AliCluster *GetCluster(Int_t /*index*/) const {return 0;}

  /// Set trigger circuit
  void SetTriggerCircuit(TClonesArray* circuit) {fTriggerCircuit = circuit;}
  /// Set pointer to data container
  void SetMUONData(AliMUONRecData* data) {fMUONData = data;}
  /// Set option
  void SetOption(Option_t* opt);

private:
  /// Not implemented
  AliMUONTracker(const AliMUONTracker& rhs);
  /// Not implemented
  AliMUONTracker& operator=(const AliMUONTracker& rhs);
    
  TClonesArray* fTriggerCircuit;                //!< trigger circuit
  AliMUONRecData*  fMUONData;                   //!< pointer to container
  AliMUONVTrackReconstructor* fTrackReco;       //!< track reconstructor

  ClassDef(AliMUONTracker,0)  //tracker base class for MUON
};
#endif
