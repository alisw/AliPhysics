#ifndef ALIMUONTRIGGERTRACKSTOREV1_H
#define ALIMUONTRIGGERTRACKSTOREV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONTriggerTrackStoreV1
/// \brief Implementation of AliMUONVTriggerTrackStore
/// 
// Author Laurent Aphecetche

#ifndef ALIMUONVTRIGERTRACKSTORE_H
#  include "AliMUONVTriggerTrackStore.h"
#endif

class TClonesArray;

class AliMUONTriggerTrackStoreV1 : public AliMUONVTriggerTrackStore
{
public:
  AliMUONTriggerTrackStoreV1();
  AliMUONTriggerTrackStoreV1(TRootIOCtor* dummy);
  virtual ~AliMUONTriggerTrackStoreV1();
  
  using AliMUONVTriggerTrackStore::Add;
  virtual void Add(const AliMUONTriggerTrack& track);

  /// Whether the Connect(TTree&) method is implemented
  virtual Bool_t CanConnect() const { return kTRUE; }
  
  virtual void Clear(Option_t* opt="");

  virtual Bool_t Connect(TTree& tree, Bool_t alone=kTRUE) const;

  virtual AliMUONTriggerTrackStoreV1* Create() const { return new AliMUONTriggerTrackStoreV1; }
  
  virtual TIterator* CreateIterator() const;  

  using AliMUONVTriggerTrackStore::GetSize;
  
  virtual Int_t GetSize() const;
  
private:
    /// Return the tracks array
    TClonesArray* Tracks() const { return fTracks; }
    /// Return the address of the tracks array
    TClonesArray** TracksPtr() const { return const_cast<TClonesArray**>(&fTracks); }
  
    /// Not implemented
    AliMUONTriggerTrackStoreV1(const AliMUONTriggerTrackStoreV1&);
    /// Not implemented
    AliMUONTriggerTrackStoreV1& operator=(const AliMUONTriggerTrackStoreV1&);
  
private:
    TClonesArray* fTracks; ///< internal array
  
  ClassDef(AliMUONTriggerTrackStoreV1,1) // Implementation of AliMUONVTriggerTrackStore
};

#endif
