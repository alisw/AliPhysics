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

class AliMUONV2DStore;
class TObjArray;

class AliMUONGainSubprocessor : public AliMUONVSubprocessor
{
public:
  AliMUONGainSubprocessor(AliMUONPreprocessor* master);
  virtual ~AliMUONGainSubprocessor();
  
  void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  UInt_t Process(TMap* dcsAliasMap);
  void Print(Option_t* opt="") const;
  
private:
  /// Not implemented
  AliMUONGainSubprocessor(const AliMUONGainSubprocessor&);
  /// Not implemented
  AliMUONGainSubprocessor& operator=(const AliMUONGainSubprocessor&);
  
  Int_t ReadFile(const char* filename);

private:
  AliMUONV2DStore* fGains; //!< Gains for the MUON TRK
  
  ClassDef(AliMUONGainSubprocessor,1) // A shuttle preprocessor for MUON TRK gains
};

#endif
