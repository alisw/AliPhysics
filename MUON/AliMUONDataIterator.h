/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONDataIterator
/// \brief An iterator on MUON data structures (so far only Digits).
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUONDATAITERATOR_H
#define ALIMUONDATAITERATOR_H

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONData;
class AliMUONVDataIterator;

class AliMUONDataIterator : public TObject
{
public:
  /// Iteration style
  enum EIterationStyle { kAllChambers, kTrackingChambers, kTriggerChambers };
  
  AliMUONDataIterator();
  AliMUONDataIterator(AliMUONData* data, const char* onWhatToIterate,
                      EIterationStyle howToIterate);
  virtual ~AliMUONDataIterator();
    
  TObject* Next();
  
  Bool_t Remove();
  
  void Reset();
  
private:
  AliMUONVDataIterator* fIterator; //!< the real worker   

private:
  /// Not implemented
  AliMUONDataIterator(const AliMUONDataIterator& rhs);
  /// Not implemented
  AliMUONDataIterator& operator=(const AliMUONDataIterator& rhs);
  
  ClassDef(AliMUONDataIterator,0) // Iterator on MUON data structures.
};

#endif
