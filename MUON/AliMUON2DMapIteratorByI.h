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

#ifndef ROOT_TExMap
#  include "TExMap.h"
#endif
#ifndef ROOT_TIterator
#  include "TIterator.h"
#endif

class AliMpExMap;

//_____________________________________________________________________________
class AliMUON2DMapIteratorByI : public TIterator
{
public:
  AliMUON2DMapIteratorByI(const AliMpExMap& theMap, 
                          Int_t firstI, Int_t lastI);
  AliMUON2DMapIteratorByI(const AliMUON2DMapIteratorByI& rhs);
  TIterator& operator=(const TIterator& rhs);
  AliMUON2DMapIteratorByI& operator=(const AliMUON2DMapIteratorByI& rhs);
  
  virtual ~AliMUON2DMapIteratorByI();
  
  ///The returned object must not be deleted by the user.  
  virtual TObject* Next();
  
  virtual void Reset(); 

  virtual const TCollection* GetCollection() const;
  
private:
  
  TObject* GetValue(TExMapIter& iter, Int_t& key) const;
  
private:
  const AliMpExMap* fkMap; //!< map to iterate upon
  TExMapIter* fIter2; //!< second iterator
  Int_t fCurrentI; //!< current index in direction i 
  Int_t fCurrentJ; //!< current index in direction j
  Int_t fFirstI; //!< first I to iterate upon
  Int_t fLastI; //!< last I to iterate upon
  
  ClassDef(AliMUON2DMapIteratorByI,0) // VDataIterator for 2D maps
};


#endif
