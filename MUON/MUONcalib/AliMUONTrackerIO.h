#ifndef ALIMUONTRACKERIO_H
#define ALIMUONTRACKERIO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONTrackerIO
/// \brief Converts ASCII calibration files (ped, config, occupancy) into AliMUONVStore object
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONVStore;
class TString;

using std::ofstream;

class AliMUONTrackerIO : public TObject
{
public:
  AliMUONTrackerIO();
  virtual ~AliMUONTrackerIO();

  static Int_t ReadConfig(const char* filename, AliMUONVStore& confStore);
  static Int_t DecodeConfig(const char* data, AliMUONVStore& confStore);
  static Int_t WriteConfig(ofstream& out, const AliMUONVStore& confStore);
  
  static Int_t ReadPedestals(const char* filename, AliMUONVStore& pedStore);
  static Int_t DecodePedestals(const char* data, AliMUONVStore& pedStore);
  
  static Int_t ReadOccupancy(const char* filename, AliMUONVStore& occupancyMap);
  static Int_t DecodeOccupancy(const char* data, AliMUONVStore& occupancyMap);
  
  /// Error code constants
  enum ErrorCode
  {
    kCannotOpenFile = -1, /// cannot open given file
    kDummyFile = -2, /// file is a dummy one (e.g. some intermediate gain files from the DA)
    kFormatError = -3, /// file is not of the expected format
    kNoInfoFile = -4, /// file is "empty", i.e. contains to information but that's normal
    kNoMapping = -99 /// mapping not loaded, cannot work
  };
  
  ClassDef(AliMUONTrackerIO,2) // Calibration ASCII file reader for MUON tracker
};

#endif
