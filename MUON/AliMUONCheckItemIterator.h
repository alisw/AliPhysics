#ifndef ALIMUONCHECKITEMITERATOR_H
#define ALIMUONCHECKITEMITERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONCheckItemIterator
/// \brief Iterator on CheckItem
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONCheckItem;
class TExMapIter;

class AliMUONCheckItemIterator : public TObject
{
public:
  AliMUONCheckItemIterator();
  AliMUONCheckItemIterator(const AliMUONCheckItem& object);
  virtual ~AliMUONCheckItemIterator();

  void First();
  TObject* Next();
  
private:
  /// Not implemented
  AliMUONCheckItemIterator(const AliMUONCheckItemIterator&);
  /// Not implemented
  AliMUONCheckItemIterator& operator=(const AliMUONCheckItemIterator&);
  
private:
  TExMapIter* fIter; //!< the actual iterator doing the job
  
  ClassDef(AliMUONCheckItemIterator,1) // Iterator for AliMUONCheckItem objects
};

#endif
