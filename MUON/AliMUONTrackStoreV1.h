#ifndef ALIMUONTRACKSTOREV1_H
#define ALIMUONTRACKSTOREV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONTrackStoreV1
/// \brief Implementation of AliMUONVTrackStore
/// 
// Author Laurent Aphecetche

#ifndef ALIMUONVTRACKSTORE_H
#  include "AliMUONVTrackStore.h"
#endif

class TClonesArray;

class AliMUONTrackStoreV1 : public AliMUONVTrackStore
{
public:
  AliMUONTrackStoreV1();
  AliMUONTrackStoreV1(TRootIOCtor* dummy);
  virtual ~AliMUONTrackStoreV1();
  
  using AliMUONVTrackStore::Add;
  
  virtual AliMUONTrack* Add(const AliMUONTrack& track);

  virtual AliMUONTrack* Remove(AliMUONTrack& track);
  
  /// Whether the Connect(TTree&) method is implemented
  virtual Bool_t CanConnect() const { return kTRUE; }
  
  virtual Bool_t Connect(TTree& tree, Bool_t alone=kTRUE) const;

  virtual AliMUONTrackStoreV1* Create() const { return new AliMUONTrackStoreV1; }
  
  virtual TIterator* CreateIterator() const;
  
  virtual void Clear(Option_t* opt="");
  
  using AliMUONVTrackStore::GetSize;
  
  virtual Int_t GetSize() const;
  
private:
    /// Return the tracks array
    TClonesArray* Tracks() const { return fTracks; }
  
    void CreateTracks(); 
  
    /// Return the address of the tracks array
    TClonesArray** TracksPtr() const { return const_cast<TClonesArray**>(&fTracks); }

    /// Not implemented
    AliMUONTrackStoreV1(const AliMUONTrackStoreV1&);
    /// Not implemented
    AliMUONTrackStoreV1& operator=(const AliMUONTrackStoreV1&);
  
private:
    TClonesArray* fTracks; ///< Internal array
  
  ClassDef(AliMUONTrackStoreV1,1) // Implementation of AliMUONVTrackStore
};

#endif
