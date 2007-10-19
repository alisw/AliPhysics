#ifndef ALIMUONTRIGGERPREPROCESSOR_H
#define ALIMUONTRIGGERPREPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONTriggerPreprocessor
/// \brief Shuttle preprocessor for MUON trigger
/// 
//  Author Laurent Aphecetche, Subatech

#include "AliMUONPreprocessor.h"

class AliMUONTriggerSubprocessor;

class AliMUONTriggerPreprocessor : public AliMUONPreprocessor
{
public:
  AliMUONTriggerPreprocessor(AliShuttleInterface* shuttle);
  virtual ~AliMUONTriggerPreprocessor();
  
  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);

private:
  /// Not implemented
  AliMUONTriggerPreprocessor(const AliMUONTriggerPreprocessor& rhs);
  /// Not implemented
  AliMUONTriggerPreprocessor& operator=(const AliMUONTriggerPreprocessor& rhs);
  

private:
  AliMUONTriggerSubprocessor* fTriggerSubprocessor; //!< the real worker class
  
  ClassDef(AliMUONTriggerPreprocessor,1) // MUON Trigger Shuttle preprocessor
};

#endif
