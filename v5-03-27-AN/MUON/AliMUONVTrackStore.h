#ifndef ALIMUONVTRACKSTORE_H
#define ALIMUONVTRACKSTORE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONVTrackStore
/// \brief Base class of a track container
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVSTORE_H
#  include "AliMUONVStore.h"
#endif

class AliMUONTrack;

class AliMUONVTrackStore : public AliMUONVStore
{
public:
  AliMUONVTrackStore();
  virtual ~AliMUONVTrackStore();
  
  /// Add an object, if of the right type
  virtual Bool_t Add(TObject* object);
  
  /// Add a track
  virtual AliMUONTrack* Add(const AliMUONTrack& track) = 0;
  
  /// Remove a track from the store
  virtual AliMUONTrack* Remove(AliMUONTrack& track) = 0;
  
  using AliMUONVStore::Create;
  
  /// Create a store from the tree (if possible).
  static AliMUONVTrackStore* Create(TTree& tree);

  /// Create an iterator to loop over tracks
  virtual TIterator* CreateIterator() const = 0;
  
  using AliMUONVStore::GetSize;
  
  ClassDef(AliMUONVTrackStore,1) // Base class of a track store
};

#endif
