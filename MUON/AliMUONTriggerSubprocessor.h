#ifndef ALIMUONTRIGGERSUBPROCESSOR_H
#define ALIMUONTRIGGERSUBPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONTriggerSubprocessor
/// \brief Implementation of AliMUONVSubprocessor for MUON TRK masks
/// 
//  Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVSUBPROCESSOR_H
#  include "AliMUONVSubprocessor.h"
#endif

class AliMUONTriggerLut;
class AliMUONVStore;
class AliMUONVCalibParam;
class TString;

class AliMUONTriggerSubprocessor : public AliMUONVSubprocessor
{
public:
  AliMUONTriggerSubprocessor(AliMUONPreprocessor* master);
  virtual ~AliMUONTriggerSubprocessor();
  
  void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  UInt_t Process(TMap* dcsAliasMap);
  
private:

  TString GetFileName(const char* fid) const;
  
  /// Not implemented
  AliMUONTriggerSubprocessor(const AliMUONTriggerSubprocessor&);
  /// Not implemented
  AliMUONTriggerSubprocessor& operator=(const AliMUONTriggerSubprocessor&);
  
private:
  AliMUONVStore* fRegionalMasks; //!< regional masks
  AliMUONVStore* fLocalMasks; //!< local masks
  AliMUONVCalibParam* fGlobalMasks; //!< global masks
  AliMUONTriggerLut* fLUT; //!< look-up table(s)
  
  ClassDef(AliMUONTriggerSubprocessor,1) // A shuttle preprocessor for MUON TRK masks
};

#endif
