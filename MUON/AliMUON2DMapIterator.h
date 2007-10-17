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

#ifndef ROOT_TIterator
#  include "TIterator.h"
#endif

class AliMpExMap;

//_____________________________________________________________________________
class AliMUON2DMapIterator : public TIterator
{
public:
  AliMUON2DMapIterator(AliMpExMap* theMap);
  AliMUON2DMapIterator(const AliMUON2DMapIterator& rhs);
  AliMUON2DMapIterator& operator=(const AliMUON2DMapIterator& rhs);
  TIterator& operator=(const TIterator& rhs);
  
  virtual ~AliMUON2DMapIterator();
  
  ///The returned object must not be deleted by the user.  
  virtual TObject* Next();
  
  virtual void Reset(); 
  
  virtual const TCollection* GetCollection() const;
  
private:
    
  AliMpExMap* Map(Int_t i) const;
  
private:
  AliMpExMap* fMap;        ///< Top map we iterate upon
  AliMpExMap* fCurrentMap; ///< Current map (inside top map) we are iterating upon
  Int_t fI;                ///< Map(fI) is fCurrentMap  
  Int_t fJ;                ///< Current position in fCurrentMap
  
  ClassDef(AliMUON2DMapIterator,0) // VDataIterator for 2D maps
};


#endif
