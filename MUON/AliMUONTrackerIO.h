#ifndef ALIMUONTRACKERIO_H
#define ALIMUONTRACKERIO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONTrackerIO
/// \brief Converts ASCII calibration files (ped, gains, capa) into AliMUONVStore object
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONVStore;
class TString;

class AliMUONTrackerIO : public TObject
{
public:
  AliMUONTrackerIO();
  virtual ~AliMUONTrackerIO();

  static Int_t ReadPedestals(const char* filename, AliMUONVStore& pedStore);
  static Int_t DecodePedestals(TString data, AliMUONVStore& pedStore);
  
  static Int_t ReadGains(const char* filename, AliMUONVStore& gainStore, TString& comment);
  static Int_t DecodeGains(TString data, AliMUONVStore& gainStore, TString& comment);
  
  static Int_t ReadCapacitances(const char* filename, AliMUONVStore& capaStore);
  
  /// Error code constants
  enum ErrorCode
  {
    kCannotOpenFile = -1, /// cannot open given file
    kDummyFile = -2, /// file is a dummy one (e.g. some intermediate gain files from the DA)
    kFormatError = -3 /// file is not of the expected format
  };
  
  ClassDef(AliMUONTrackerIO,1) // Calibration ASCII file reader for MUON tracker
};

#endif
