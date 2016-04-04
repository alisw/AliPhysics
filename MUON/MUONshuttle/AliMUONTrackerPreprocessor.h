#ifndef ALIMUONTRACKERPREPROCESSOR_H
#define ALIMUONTRACKERPREPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONTrackerPreprocessor
/// \brief Shuttle preprocessor for MUON tracker
///
//  Author Laurent Aphecetche

#include "AliMUONPreprocessor.h"

class AliMUONVSubprocessor;

class TObjArray;

class AliMUONTrackerPreprocessor : public AliMUONPreprocessor
{
public:
  AliMUONTrackerPreprocessor(AliShuttleInterface* shuttle);
  virtual ~AliMUONTrackerPreprocessor();

  virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);

private:
  /// Not implemented
  AliMUONTrackerPreprocessor(const AliMUONTrackerPreprocessor& rhs);
  /// Not implemented
  AliMUONTrackerPreprocessor& operator=(const AliMUONTrackerPreprocessor& rhs);

private:
  AliMUONVSubprocessor* fPedestalSubprocessor; ///< Pedestal subprocessor
  AliMUONVSubprocessor* fGMSSubprocessor; ///< GMS subprocessor
  AliMUONVSubprocessor* fHVSubprocessor;  ///< HV subprocessor
  AliMUONVSubprocessor* fOccupancySubprocessor; ///< Occupancy subprocessor
  AliMUONVSubprocessor* fBusPatchEvolutionSubprocessor; ///< Buspatch evolution subprocessor
  AliMUONVSubprocessor* fConfigSubprocessor; ///< config subprocessor
  AliMUONVSubprocessor* fLVSubprocessor;  ///< LV subprocessor

  ClassDef(AliMUONTrackerPreprocessor,7) // MUON Tracker Shuttle preprocessor
};

#endif
