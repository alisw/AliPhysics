#ifndef ALIMUON2DMAPITERATORBYI_H
#define ALIMUON2DMAPITERATORBYI_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup core
/// \class AliMUON2DMapIteratorByI
/// \brief Implementation of TIterator for 2D maps
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TIterator
#  include "TIterator.h"
#endif

class AliMpExMap;
class AliMpExMapIterator;

//_____________________________________________________________________________
class AliMUON2DMapIteratorByI : public TIterator
{
public:
  AliMUON2DMapIteratorByI(const AliMpExMap& theMap, 
                          Int_t firstI, Int_t lastI);
  
  virtual ~AliMUON2DMapIteratorByI();
  
  ///The returned object must not be deleted by the user.  
  virtual TObject* Next();
  
  virtual void Reset(); 

  virtual const TCollection* GetCollection() const;

private:
    AliMpExMapIterator* NextIterator();
  
private:
  /// Not implemented
  AliMUON2DMapIteratorByI(const AliMUON2DMapIteratorByI& rhs);
  /// Not implemented
  AliMUON2DMapIteratorByI& operator=(const AliMUON2DMapIteratorByI& rhs);
  /// Overriden TIterator virtual operator=
  AliMUON2DMapIteratorByI& operator=(const TIterator& rhs);

  const AliMpExMap* fkMap; ///< Top map we iterate upon
  AliMpExMapIterator* fIter1; ///< first iterator
  TIterator* fIter2; ///< second iterator
  Int_t fFirstI; ///< start of range for I
  Int_t fLastI; ///< end of range for I
  Int_t fCurrentI; ///< current value of I 
  
  ClassDef(AliMUON2DMapIteratorByI,0) // VDataIterator for 2D maps
};


#endif
