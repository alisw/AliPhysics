#ifndef ALIMUONHVSUBPROCESSOR_H
#define ALIMUONHVSUBPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONHVSubprocessor
/// \brief A subprocessor to read HV values for one run
/// \author Laurent Aphecetche

#ifndef ALIMUONVSUBPROCESSOR_H
#  include "AliMUONVSubprocessor.h"
#endif

class AliMUONHVSubprocessor : public AliMUONVSubprocessor
{
public:
  AliMUONHVSubprocessor(AliMUONPreprocessor* master);
  virtual ~AliMUONHVSubprocessor();
  
  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  virtual UInt_t Process(TMap* dcsAliasMap);

private:
  AliMUONHVSubprocessor(const AliMUONHVSubprocessor&);
  AliMUONHVSubprocessor& operator=(const AliMUONHVSubprocessor&);
  
  ClassDef(AliMUONHVSubprocessor,1) // Shuttle Subprocessor for MUON TRK HV
};

#endif
