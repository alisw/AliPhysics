#ifndef ALIMUONTRACKERPREPROCESSOR_H
#define ALIMUONTRACKERPREPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONTrackerPreprocessor
/// \brief Shuttle preprocessor for MUON tracker
/// 
/// \author Laurent Aphecetche

#include "AliMUONPreprocessor.h"

class AliMUONPedestalSubprocessor;
class AliMUONGMSSubprocessor;
class AliMUONHVSubprocessor;

class TObjArray;

class AliMUONTrackerPreprocessor : public AliMUONPreprocessor
{
public:
  AliMUONTrackerPreprocessor(AliShuttleInterface* shuttle);
  virtual ~AliMUONTrackerPreprocessor();
  
  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);

private:
  AliMUONTrackerPreprocessor(const AliMUONTrackerPreprocessor& rhs);
  AliMUONTrackerPreprocessor& operator=(const AliMUONTrackerPreprocessor& rhs);
  
private:
  AliMUONPedestalSubprocessor* fPedestalSubprocessor; ///< Pedestal subprocessor
  AliMUONGMSSubprocessor*      fGMSSubprocessor;      ///< GMS subprocessor
  AliMUONHVSubprocessor*       fHVSubprocessor;       ///< HV subprocessor

  ClassDef(AliMUONTrackerPreprocessor,1) // MUON Tracker Shuttle preprocessor
};

#endif
