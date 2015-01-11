#ifndef ALIMUONPEDESTALSUBPROCESSOR_H
#define ALIMUONPEDESTALSUBPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONPedestalSubprocessor
/// \brief Implementation of AliMUONVSubprocessor for MUON TRK pedestals
/// 
//  Author Laurent Aphecetche

#ifndef ALIMUONVSUBPROCESSOR_H
#  include "AliMUONVSubprocessor.h"
#endif

class AliMUONVStore;
class TObjArray;

class AliMUONPedestalSubprocessor : public AliMUONVSubprocessor
{
public:
  AliMUONPedestalSubprocessor(AliMUONPreprocessor* master);
  virtual ~AliMUONPedestalSubprocessor();
  
  Bool_t Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  UInt_t Process(TMap* dcsAliasMap);
  void Print(Option_t* opt="") const;
  
private:
  /// Not implemented
  AliMUONPedestalSubprocessor(const AliMUONPedestalSubprocessor&);
  /// Not implemented
  AliMUONPedestalSubprocessor& operator=(const AliMUONPedestalSubprocessor&);
  
  Int_t ReadPedestalFile(const char* filename);
  Int_t ReadConfigFile(const char* filename);

  Bool_t HasConfigChanged(const AliMUONVStore& newConfig) const;

private:
  AliMUONVStore* fPedestals; //!< Pedestals for the MUON TRK
  AliMUONVStore* fConfig; //!< Configuration (i.e. list of (buspatch,manu)) for the MUON TRK
  Bool_t fConfigChanged; //!< flag to trigger the saving of the configuration
  Bool_t fTooFewEvents; //!< whether the current run was a failed ped run, basically
  
  ClassDef(AliMUONPedestalSubprocessor,3) // A shuttle preprocessor for MUON TRK pedestals
};

#endif
