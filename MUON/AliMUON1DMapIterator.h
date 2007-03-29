#ifndef ALIMUON1DMAPITERATOR_H
#define ALIMUON1DMAPITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUON1DMapIterator
/// \brief Implementation of AliMUONVDataIterator for 1D maps
/// 
//  Author Laurent Aphecetche

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
class AliMUON1DMapIterator : public AliMUONVDataIterator
{
public:
  AliMUON1DMapIterator(AliMpExMap& theMap);
  
  virtual ~AliMUON1DMapIterator();
  
  /** The object returned by this iterator is an AliMUONObjectPair(TObject* key,TObject* value)
    where key is an AliMpIntPair (i,0), and value is 
    an AliMUONVCalibParam.
    The returned object must be deleted by the user.
    */
  virtual TObject* Next();
  
  virtual void Reset(); 
  
  virtual Bool_t Remove();
  
private:
  /// copy ctor will not implemented
  AliMUON1DMapIterator(const AliMUON1DMapIterator&);
  /// assignement operator will not implemented
  AliMUON1DMapIterator& operator=(const AliMUON1DMapIterator&);
    
private:
  TExMapIter fIter; //!< iterator
  Int_t fCurrentI; //!< current index in direction i 
  
  ClassDef(AliMUON1DMapIterator,0) // VDataIterator for 1D maps
};


#endif
