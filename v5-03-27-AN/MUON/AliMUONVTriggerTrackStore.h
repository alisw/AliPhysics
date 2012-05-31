#ifndef ALIMUONVTRIGGERTRACKSTORE_H
#define ALIMUONVTRIGGERTRACKSTORE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONVTriggerTrackStore
/// \brief Base class of a trigger track store
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVSTORE_H
#  include "AliMUONVStore.h"
#endif

class AliMUONTriggerTrack;

class AliMUONVTriggerTrackStore : public AliMUONVStore
{
public:
  AliMUONVTriggerTrackStore();
  virtual ~AliMUONVTriggerTrackStore();

  /// Add 
  virtual Bool_t Add(TObject* object);
  
  /// Add a trigger track
  virtual void Add(const AliMUONTriggerTrack& track) = 0;
  
  using AliMUONVStore::Create;
  
  /// Create a store from the tree (if possible).
  static AliMUONVTriggerTrackStore* Create(TTree& tree);

  /// Iterator to loop over tracks
  virtual TIterator* CreateIterator() const = 0;

  using AliMUONVStore::GetSize;
  
  ClassDef(AliMUONVTriggerTrackStore,1) // Base class of a trigger track container
};

#endif
