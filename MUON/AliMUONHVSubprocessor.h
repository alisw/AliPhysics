#ifndef ALIMUONHVSUBPROCESSOR_H
#define ALIMUONHVSUBPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONHVSubprocessor
/// \brief A subprocessor to read HV values for one run
///
//  Author Laurent Aphecetche

#ifndef ALIMUONVSUBPROCESSOR_H
#  include "AliMUONVSubprocessor.h"
#endif

class AliMUONHVSubprocessor : public AliMUONVSubprocessor
{
public:
  AliMUONHVSubprocessor(AliMUONPreprocessor* master);
  virtual ~AliMUONHVSubprocessor();
  
  virtual UInt_t Process(TMap* dcsAliasMap);

private:
  /// Not implemented
  AliMUONHVSubprocessor(const AliMUONHVSubprocessor&);
  /// Not implemented
  AliMUONHVSubprocessor& operator=(const AliMUONHVSubprocessor&);
  
  ClassDef(AliMUONHVSubprocessor,1) // Shuttle Subprocessor for MUON TRK HV
};

#endif
