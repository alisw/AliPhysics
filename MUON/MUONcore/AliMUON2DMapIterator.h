#ifndef ALIMUON2DMAPITERATOR_H
#define ALIMUON2DMAPITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup core
/// \class AliMUON2DMapIterator
/// \brief Implementation of TIterator for 2D maps
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TIterator
#  include "TIterator.h"
#endif

class AliMpExMap;

//_____________________________________________________________________________
class AliMUON2DMapIterator : public TIterator
{
public:
  AliMUON2DMapIterator(const AliMpExMap& theMap);
  
  virtual ~AliMUON2DMapIterator();
  
  ///The returned object must not be deleted by the user.  
  virtual TObject* Next();
  
  virtual void Reset(); 
  
  virtual const TCollection* GetCollection() const;
  
private:
  TIterator* NextIterator();
  
private:
  /// Not implemented
  AliMUON2DMapIterator(const AliMUON2DMapIterator& rhs);
  /// Not implemented
  AliMUON2DMapIterator& operator=(const AliMUON2DMapIterator& rhs);
  /// Overriden TIterator virtual operator=
  AliMUON2DMapIterator& operator=(const TIterator& rhs);

  const AliMpExMap* fkMap; ///< Top map we iterate upon
  TIterator* fIter1; ///< first iterator
  TIterator* fIter2; ///< second iterator
  
  ClassDef(AliMUON2DMapIterator,0) // TIterator for AliMUON2D maps
};


#endif
