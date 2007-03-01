#ifndef ALIMUONTRACKER_H
#define ALIMUONTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/// \ingroup base
/// \class AliMUONTracker
/// \brief MUON base Tracker

#include "AliTracker.h"
class AliESD;
class AliCluster;
class AliMUONData;
class TClonesArray;
class AliMUONVTrackReconstructor;

class AliMUONTracker : public AliTracker
{
 public:

  AliMUONTracker();
  virtual ~AliMUONTracker();
    
  virtual Int_t Clusters2Tracks(AliESD* /*esd*/); 
  
  virtual Int_t PropagateBack(AliESD* /*event*/) {return 0;}
  virtual Int_t RefitInward(AliESD* /*event*/) {return 0;}
  virtual Int_t LoadClusters(TTree* /*tree*/) {return 0;}
  virtual void  UnloadClusters() {return;}
  virtual AliCluster *GetCluster(Int_t /*index*/) const {return 0;}

  void SetTriggerCircuit(TClonesArray* circuit) {fTriggerCircuit = circuit;}
  void SetMUONData(AliMUONData* data) {fMUONData = data;}
  void SetOption(Option_t* opt);

private:

  TClonesArray* fTriggerCircuit;                //!< trigger circuit
  AliMUONData*  fMUONData;                      //!< pointer to container
  AliMUONVTrackReconstructor* fTrackReco;       //!< track reconstructor

  AliMUONTracker(const AliMUONTracker& rhs);
  AliMUONTracker& operator=(const AliMUONTracker& rhs);
    
  ClassDef(AliMUONTracker,0)  //tracker base class for MUON
};
#endif
