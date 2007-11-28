#ifndef ALIMUON1DMAPITERATOR_H
#define ALIMUON1DMAPITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup core
/// \class AliMUON1DMapIterator
/// \brief Implementation of TIterator for 1D maps
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TExMap
#  include "TExMap.h"
#endif
#ifndef ROOT_TIterator
#  include "TIterator.h"
#endif

class AliMpExMap;

//_____________________________________________________________________________
class AliMUON1DMapIterator : public TIterator
{
public:
  AliMUON1DMapIterator(AliMpExMap& theMap);
  AliMUON1DMapIterator(const AliMUON1DMapIterator&);
  AliMUON1DMapIterator& operator=(const AliMUON1DMapIterator& rhs);
  TIterator& operator=(const TIterator& iterator);  
  virtual ~AliMUON1DMapIterator();
  
  /** The returned object must not be deleted by the user ! */
  virtual TObject* Next();
  
  virtual void Reset(); 
  
  /// Return 0 as we're not really dealing with a TCollection
  virtual const TCollection* GetCollection() const { return 0x0; }
  
private:
  TExMapIter fIter; //!< iterator
  Int_t fCurrentI; //!< current index in direction i 
  
  ClassDef(AliMUON1DMapIterator,0) // VDataIterator for 1D maps
};


#endif
