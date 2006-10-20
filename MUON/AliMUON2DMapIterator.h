#ifndef ALIMUON2DMAPITERATOR_H
#define ALIMUON2DMAPITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUON2DMapIterator
/// \brief Implementation of AliMUONVDataIterator for 2D maps
/// 
/// \author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TExMap
#  include "TExMap.h"
#endif
#ifndef ALIMUONVDATAITERATOR_H
#  include "AliMUONVDataIterator.h"
#endif

class AliMpExMap;

//_____________________________________________________________________________
class AliMUON2DMapIterator : public AliMUONVDataIterator
{
public:
  AliMUON2DMapIterator(AliMpExMap& theMap);
  
  virtual ~AliMUON2DMapIterator();
  
  /** The object returned by this iterator is an AliMUONObjectPair(TObject* key,TObject* value)
    where key is an AliMpIntPair (detElemId,manuId), and value is 
    an AliMUONVCalibParam.
    The returned object must be deleted by the user.
    */
  virtual TObject* Next();
  
  virtual void Reset(); 
  
  virtual Bool_t Remove();
  
private:
    // copy ctor will not implemented
    AliMUON2DMapIterator(const AliMUON2DMapIterator&);
  // assignement operator will not implemented
  AliMUON2DMapIterator& operator=(const AliMUON2DMapIterator&);
  
    TObject* GetValue(TExMapIter& iter, Int_t& key) const;
  AliMpExMap* GetMap(TExMapIter& iter, Int_t& key);
  
private:
    TExMapIter fIter; //! first iterator
  TExMapIter* fIter2; //! second iterator
  Int_t fCurrentI; //! current index in direction i 
  Int_t fCurrentJ; //! current index in direction j
  
  ClassDef(AliMUON2DMapIterator,0) // VDataIterator for 2D maps
};


#endif
