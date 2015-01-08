#ifndef ALIMUONTRIGGERDCSSUBPROCESSOR_H
#define ALIMUONTRIGGERDCSSUBPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONTriggerDCSSubprocessor
/// \brief A subprocessor to read TriggerDCS values for one run
///
//  Author Laurent Aphecetche

#ifndef ALIMUONVSUBPROCESSOR_H
#  include "AliMUONVSubprocessor.h"
#endif

class AliMUONTriggerDCSSubprocessor : public AliMUONVSubprocessor
{
public:
  AliMUONTriggerDCSSubprocessor(AliMUONPreprocessor* master);
  virtual ~AliMUONTriggerDCSSubprocessor();
  
  virtual UInt_t Process(TMap* dcsAliasMap);

private:
  /// Not implemented
  AliMUONTriggerDCSSubprocessor(const AliMUONTriggerDCSSubprocessor&);
  /// Not implemented
  AliMUONTriggerDCSSubprocessor& operator=(const AliMUONTriggerDCSSubprocessor&);
  
  ClassDef(AliMUONTriggerDCSSubprocessor,1) // Shuttle Subprocessor for MUON TRG HV and currents
};

#endif
