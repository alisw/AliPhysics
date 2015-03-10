#ifndef ALIMUONOCCUPANCYSUBPROCESSOR_H
#define ALIMUONOCCUPANCYSUBPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONOccupancySubprocessor
/// \brief Implementation of AliMUONVSubprocessor for MUON TRK occupancy
/// 
// Author Laurent Aphecetche

#ifndef ALIMUONVSUBPROCESSOR_H
#  include "AliMUONVSubprocessor.h"
#endif

class AliMUONVStore;

class AliMUONOccupancySubprocessor : public AliMUONVSubprocessor
{
public:
  AliMUONOccupancySubprocessor(AliMUONPreprocessor* master);
  virtual ~AliMUONOccupancySubprocessor();
  
  Bool_t Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  UInt_t Process(TMap* dcsAliasMap);
  void Print(Option_t* opt="") const;

private:
  
  /// Not implemented
  AliMUONOccupancySubprocessor(const AliMUONOccupancySubprocessor&);
  /// Not implemented
  AliMUONOccupancySubprocessor& operator=(const AliMUONOccupancySubprocessor&);
  
  Int_t ReadFile(const char* filename);
  
private:
  AliMUONVStore* fOccupancyMap; //!<! Occupancy map (at the manu level) for the MUON TRK
  
  ClassDef(AliMUONOccupancySubprocessor,1) // A shuttle preprocessor for MUON TRK occupancy  
};

#endif
