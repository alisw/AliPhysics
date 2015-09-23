#ifndef ALIMUONCONFIGSUBPROCESSOR_H
#define ALIMUONCONFIGSUBPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONConfigSubprocessor
/// \brief Implementation of AliMUONVSubprocessor for MUON TRK config
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUONVSUBPROCESSOR_H
#  include "AliMUONVSubprocessor.h"
#endif

class AliMUONVStore;

class AliMUONConfigSubprocessor : public AliMUONVSubprocessor
{
public:
  AliMUONConfigSubprocessor(AliMUONPreprocessor* master);
  virtual ~AliMUONConfigSubprocessor();
  
  Bool_t Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  UInt_t Process(TMap* dcsAliasMap);
  void Print(Option_t* opt="") const;
  
private:
  /// Not implemented
  AliMUONConfigSubprocessor(const AliMUONConfigSubprocessor&);
  /// Not implemented
  AliMUONConfigSubprocessor& operator=(const AliMUONConfigSubprocessor&);
  
  Int_t ReadConfigFile(const char* filename);

  Bool_t HasConfigChanged(const AliMUONVStore& newConfig) const;

private:
  AliMUONVStore* fConfig; //!<! Configuration (i.e. list of (buspatch,manu)) for the MUON TRK
  Bool_t fConfigChanged; //!<! flag to trigger the saving of the configuration
  
  ClassDef(AliMUONConfigSubprocessor,1) // A shuttle preprocessor for MUON TRK config
};

#endif
