#ifndef ALIMUONPEDESTALEVENTGENERATOR_H
#define ALIMUONPEDESTALEVENTGENERATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONPedestalEventGenerator
/// \brief Generate pedestal events (only for tracker).
/// 
/// \author Laurent Aphecetche

#ifndef ROOT_TTask
#  include "TTask.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONCalibrationData;
class AliMUONData;
class TList;
class AliRunLoader;

class AliMUONPedestalEventGenerator : public TTask
{
public:
  AliMUONPedestalEventGenerator(Int_t runNumber, Int_t nevents, const char* dateFileName);
  virtual ~AliMUONPedestalEventGenerator();
  
  void Exec(Option_t* option);
  
  void MakeDDL(Bool_t value) { fMakeDDL = value; }
  
private:
  AliMUONPedestalEventGenerator(const AliMUONPedestalEventGenerator&);
  AliMUONPedestalEventGenerator& operator=(const AliMUONPedestalEventGenerator&);
  Bool_t ConvertRawFilesToDate();
  void GenerateDigits(AliMUONData* data);
  AliMUONData* GetDataAccess(const char* mode);
  AliRunLoader* LoadRun(const char* mode);
  void Digits2Raw();
  
private:
  TList* fManuList; //! list of (de,manu) pairs
  AliMUONCalibrationData* fCalibrationData; //! access to pedestal CDB
  TString fDateFileName; //! basefilename of the DATE output file
  TString fGAliceFileName; //! absolute path to galice.root file
  Bool_t fMakeDDL; //! whether to generate DDL ascii files or not
  static Int_t fgCounter; //! 
  
  ClassDef(AliMUONPedestalEventGenerator,1) // Random generator of pedestal events for MUON TRK
};

#endif
