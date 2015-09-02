#ifndef ALIMUONBUSPATCHEVOLUTIONSUBPROCESSOR_H
#define ALIMUONBUSPATCHEVOLUTIONSUBPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONBusPatchEvolutionSubprocessor
/// \brief Implementation of AliMUONVSubprocessor for MUON TRK bus patch evolution
/// 
// Author Laurent Aphecetche

#ifndef ALIMUONVSUBPROCESSOR_H
#  include "AliMUONVSubprocessor.h"
#endif

#include <vector>

class AliMergeableCollection;
class TH1;

class AliMUONBusPatchEvolutionSubprocessor : public AliMUONVSubprocessor
{
public:
	AliMUONBusPatchEvolutionSubprocessor(AliMUONPreprocessor* master);
  virtual ~AliMUONBusPatchEvolutionSubprocessor();
  
  Bool_t Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  UInt_t Process(TMap* dcsAliasMap);
  void Print(Option_t* opt="") const;

private:
  
  /// Not implemented
  AliMUONBusPatchEvolutionSubprocessor(const AliMUONBusPatchEvolutionSubprocessor&);
  /// Not implemented
  AliMUONBusPatchEvolutionSubprocessor& operator=(const AliMUONBusPatchEvolutionSubprocessor&);

  Bool_t ReadFile(const char* filename);
  
  void ShrinkTimeAxis(AliMergeableCollection& hc);


private:
  AliMergeableCollection* fBPEVO; //!<! Bus patch evolution
  int fProductionMode; //!<! Whether or not we are using this one in production mode
  
  ClassDef(AliMUONBusPatchEvolutionSubprocessor,1) // A shuttle preprocessor for MUON TRK bus patch  evolution
};

#endif
