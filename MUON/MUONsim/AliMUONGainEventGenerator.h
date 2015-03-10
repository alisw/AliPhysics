#ifndef ALIMUONGAINEVENTGENERATOR_H
#define ALIMUONGAINEVENTGENERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup sim
/// \class AliMUONGainEventGenerator
/// \brief Generate gain-calibration-like files
/// 
// Author Laurent Aphecetche

#ifndef ROOT_TTask
#  include "TTask.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONCalibrationData;
class AliMUONVStore;

class AliMUONGainEventGenerator : public TTask
{
public:
  AliMUONGainEventGenerator(Int_t sourceGainRunNumber,                          
                            Int_t sourcePedRunNumber,
                            Int_t nEventsPerFile, 
                            const char* dateBaseFileName);
  virtual ~AliMUONGainEventGenerator();
  
  virtual void Exec(Option_t* option);
  
private:
    /// not implemented
    AliMUONGainEventGenerator(const AliMUONGainEventGenerator& rhs);
  /// not implemented
  AliMUONGainEventGenerator& operator=(const AliMUONGainEventGenerator& rhs);
  
  void GeneratePedestals(Int_t runNumber, Float_t injection);
  void WriteToCDB(TObject* object, Int_t runNumber);

private:
  Int_t fNofEventsPerFile; //!<! number of events to generate per file
  Int_t fSourcePedestalRunNumber; //!<! run number of pedestal to be used
  TString fDateBaseFileName; //!<! base file name of the output file
  AliMUONVStore* fSourceGains; //!<! the gains used to generate the files
  AliMUONVStore* fSourcePedestals; //!<! the pedestals used to generate the fiels
  
  ClassDef(AliMUONGainEventGenerator,1) // Generate gain-like files
};

#endif
