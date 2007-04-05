/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONVDataIterator
/// \brief Defines an interface of an iterator over muon data structure(s)
/// 
//  Author Laurent Aphecetche

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
  
  /// Return next object in the iteration.
  /// To know if the client must delete that object, use
  /// the IsOwner() method
  virtual TObject* Next() = 0;
  
  /// Reset (i.e. Rewind) the iterator
  virtual void Reset() = 0; 
  
  /// Remove the current item (might not be implemented yet)
  virtual Bool_t Remove() = 0;
  
  /// Whether or not the objects returned by this iterator are to be deleted
  virtual Bool_t IsOwner() const = 0;
  
  ClassDef(AliMUONVDataIterator,0) // Interface for an iterator on AliMUONData.
};

#endif
