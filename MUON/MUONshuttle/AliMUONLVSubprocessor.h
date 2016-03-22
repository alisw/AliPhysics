#ifndef ALIMUONLVSUBPROCESSOR_H
#define ALIMUONLVSUBPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONLVSubprocessor
/// \brief A subprocessor to read LV values for one run
///
//  Author Laurent Aphecetche

#ifndef ALIMUONVSUBPROCESSOR_H
#  include "AliMUONVSubprocessor.h"
#endif

class AliMUONLVSubprocessor : public AliMUONVSubprocessor
{
public:
  AliMUONLVSubprocessor(AliMUONPreprocessor* master);
  virtual ~AliMUONLVSubprocessor();

  virtual UInt_t Process(TMap* dcsAliasMap);

private:
  /// Not implemented
  AliMUONLVSubprocessor(const AliMUONLVSubprocessor&);
  /// Not implemented
  AliMUONLVSubprocessor& operator=(const AliMUONLVSubprocessor&);

  ClassDef(AliMUONLVSubprocessor,1) // Shuttle Subprocessor for MUON TRK LV
};

#endif
