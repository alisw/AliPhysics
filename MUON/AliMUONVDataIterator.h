/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup
/// \class AliMUONVDataIterator
/// \brief
/// 
/// \author Laurent Aphecetche

#ifndef ALIMUONVDATAITERATOR_H
#define ALIMUONVDATAITERATOR_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONVDataIterator : public TObject
{
public:
  AliMUONVDataIterator();
  virtual ~AliMUONVDataIterator();
  
  virtual TObject* Next() = 0;
  
  virtual void Reset() = 0; 
  
  virtual Bool_t Remove() = 0;
  
  ClassDef(AliMUONVDataIterator,0) // Interface for an iterator on AliMUONData.
};

#endif
