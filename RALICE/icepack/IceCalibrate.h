#ifndef IceCalibrate_h
#define IceCalibrate_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TROOT.h"
#include "TTask.h"
#include "TString.h"
#include "TFile.h"

#include "AliJob.h"
#include "IceEvent.h"
#include "IceGOM.h"

class IceCalibrate : public TTask
{
 public :
  IceCalibrate(const char* name="",const char* title=""); // Constructor
  virtual ~IceCalibrate();                                // Destructor
  virtual void Exec(Option_t* opt);                       // Perform the calibrations
  void SetOMdbase(AliObjMatrix* omdb);                    // Set the OM dbase object
  void SetCalibFile(TString name);                        // Set ROOT calibration input file

 protected :
  TFile* fCalfile;     // The (optional) calibration input file in ROOT format
  AliObjMatrix* fOmdb; // The OM database object

 ClassDef(IceCalibrate,2) // TTask derived class to perform the various calibrations
};
#endif
