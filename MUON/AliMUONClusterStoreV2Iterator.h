#ifndef ALIMUONCLUSTERSTOREV2ITERATOR_H
#define ALIMUONCLUSTERSTOREV2ITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONClusterStoreV2Iterator
/// \brief Base implementation of TIterator for AliMUONClusterStoreV2
/// 
//  Author Philippe Pillot, Subatech

#ifndef ROOT_TIterator
#  include "TIterator.h"
#endif

class AliMUONClusterStoreV2;

class AliMUONClusterStoreV2Iterator : public TIterator
{
public:
  AliMUONClusterStoreV2Iterator(const AliMUONClusterStoreV2* store,
                                Int_t firstChamberId, Int_t lastChamberId);
  
  virtual ~AliMUONClusterStoreV2Iterator();
  
  TObject* Next();
  
  void Reset();
  
  /// Return 0 as we're not dealing with TCollection objects really
  virtual const TCollection* GetCollection() const { return 0x0; }
  
private:
  TObject* NextInCurrentChamber() const;
  
private:
  /// Not implemented
  AliMUONClusterStoreV2Iterator(const AliMUONClusterStoreV2Iterator& rhs);
  /// Not implemented
  AliMUONClusterStoreV2Iterator& operator=(const AliMUONClusterStoreV2Iterator& rhs);
  /// Overriden TIterator virtual operator=
  AliMUONClusterStoreV2Iterator& operator=(const TIterator& rhs);

  const AliMUONClusterStoreV2* fkStore; ///< store to iterate upon
  Int_t fFirstChamberId; ///< first chamber
  Int_t fLastChamberId; ///< last chamber
  Int_t fCurrentChamberId; ///< current chamber
  TIterator* fChamberIterator; ///< helper iterator
  
  ClassDef(AliMUONClusterStoreV2Iterator,0) // Implementation of TIterator
};

#endif
