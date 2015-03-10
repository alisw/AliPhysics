#ifndef ALIMUONGAINSUBPROCESSOR_H
#define ALIMUONGAINSUBPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONGainSubprocessor
/// \brief Implementation of AliMUONVSubprocessor for MUON TRK Gains
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUONVSUBPROCESSOR_H
#  include "AliMUONVSubprocessor.h"
#endif

#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONVStore;
class TObjArray;

class AliMUONGainSubprocessor : public AliMUONVSubprocessor
{
public:
  AliMUONGainSubprocessor(AliMUONPreprocessor* master);
  virtual ~AliMUONGainSubprocessor();
  
  Bool_t Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  UInt_t Process(TMap* dcsAliasMap);
  
private:
  /// Not implemented
  AliMUONGainSubprocessor(const AliMUONGainSubprocessor&);
  /// Not implemented
  AliMUONGainSubprocessor& operator=(const AliMUONGainSubprocessor&);
  
  Int_t ReadFile(const char* filename);

private:
  AliMUONVStore* fGains; //!<! Gains for the MUON TRK
  Bool_t fSkip; //!<! whether we should skip this run (because it's dummy)
  TString fComment; //!<! comment for OCDB entry
  
  ClassDef(AliMUONGainSubprocessor,2) // A shuttle preprocessor for MUON TRK gains
};

#endif
