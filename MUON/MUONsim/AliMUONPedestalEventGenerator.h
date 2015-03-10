#ifndef ALIMUONPEDESTALEVENTGENERATOR_H
#define ALIMUONPEDESTALEVENTGENERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup sim
/// \class AliMUONPedestalEventGenerator
/// \brief Generate pedestal events (only for tracker).
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TTask
#  include "TTask.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONCalibrationData;
class TList;
class AliRunLoader;
class AliMUONVDigitStore;
class AliLoader;
class AliMUONVStore;
class AliMUONRawWriter;

class AliMUONPedestalEventGenerator : public TTask
{
public:
  AliMUONPedestalEventGenerator(Int_t runNumber, Int_t nevents, const char* dateFileName);
  virtual ~AliMUONPedestalEventGenerator();
  
  void Exec(Option_t* option);
  
  /// Set option whether to generate DDL ascii files or not
  void MakeDDL(Bool_t value) { fMakeDDL = value; }
  
private:
  /// Not implemented
  AliMUONPedestalEventGenerator(const AliMUONPedestalEventGenerator&);
  /// Not implemented
  AliMUONPedestalEventGenerator& operator=(const AliMUONPedestalEventGenerator&);

  Bool_t ConvertRawFilesToDate();
  AliMUONVDigitStore* DigitStore();
  void GenerateDigits(AliMUONVDigitStore& digitStore);
  AliRunLoader* LoadRun(const char* mode);
  void Digits2Raw(Int_t event);
  
private:
  AliMUONCalibrationData* fCalibrationData; //!<! access to pedestal CDB
  TString fDateFileName; //!<! basefilename of the DATE output file
  TString fGAliceFileName; //!<! absolute path to galice.root file
  Bool_t fMakeDDL; //!<! whether to generate DDL ascii files or not
  AliLoader* fLoader; //!<! to access trees
  AliMUONVStore* fPedestals; //!<! pedestals
  AliMUONVDigitStore* fDigitStore; //!<! digit container
  AliMUONRawWriter* fRawWriter; //!<! to convert digits to raw data
  static Int_t fgCounter; //!<! counter 
  
  ClassDef(AliMUONPedestalEventGenerator,3) // Random generator of pedestal events for MUON TRK
};

#endif
