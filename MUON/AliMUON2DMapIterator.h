#ifndef ALIMUON2DMAPITERATOR_H
#define ALIMUON2DMAPITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUON2DMapIterator
/// \brief Implementation of TIterator for 2D maps
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TExMap
#  include "TExMap.h"
#endif
#ifndef ROOT_TIterator
#  include "TIterator.h"
#endif

class AliMpExMap;

//_____________________________________________________________________________
class AliMUON2DMapIterator : public TIterator
{
public:
  AliMUON2DMapIterator(const AliMpExMap& theMap);
  AliMUON2DMapIterator(const AliMUON2DMapIterator& rhs);
  AliMUON2DMapIterator& operator=(const AliMUON2DMapIterator& rhs);
  TIterator& operator=(const TIterator& rhs);
  
  virtual ~AliMUON2DMapIterator();
  
  ///The returned object must not be deleted by the user.  
  virtual TObject* Next();
  
  virtual void Reset(); 
  
  virtual const TCollection* GetCollection() const;
  
private:
  TObject* GetValue(TExMapIter& iter, Int_t& key) const;
  AliMpExMap* GetMap(TExMapIter& iter, Int_t& key) const;
  
private:
  TExMapIter fIter; //!< first iterator
  TExMapIter* fIter2; //!< second iterator
  Int_t fCurrentI; //!< current index in direction i 
  Int_t fCurrentJ; //!< current index in direction j

  ClassDef(AliMUON2DMapIterator,0) // VDataIterator for 2D maps
};


#endif
