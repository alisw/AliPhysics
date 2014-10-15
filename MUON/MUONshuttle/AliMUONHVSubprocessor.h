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
  AliMUONHVSubprocessor(AliMUONPreprocessor* master, Bool_t includeHVcurrents=kFALSE);
  virtual ~AliMUONHVSubprocessor();
  
  virtual UInt_t Process(TMap* dcsAliasMap);

  Bool_t IncludeHVCurrent() const { return fIncludeHVCurrents; }
  
private:
  /// Not implemented
  AliMUONHVSubprocessor(const AliMUONHVSubprocessor&);
  /// Not implemented
  AliMUONHVSubprocessor& operator=(const AliMUONHVSubprocessor&);
  
  Bool_t fIncludeHVCurrents; // whether or not to transfer also HV current (in addition to HV voltages)
  
  ClassDef(AliMUONHVSubprocessor,2) // Shuttle Subprocessor for MUON TRK HV
};

#endif
