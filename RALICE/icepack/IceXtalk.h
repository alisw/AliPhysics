#ifndef IceXtalk_h
#define IceXtalk_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TROOT.h"
#include "TTask.h"
#include "TString.h"
#include "TFile.h"

#include "AliJob.h"
#include "IceEvent.h"
#include "IceAOM.h"

class IceXtalk : public TTask
{
 public :
  IceXtalk(const char* name="",const char* title=""); // Constructor
  virtual ~IceXtalk();                                // Destructor
  virtual void Exec(Option_t* opt);                   // Cross talk hit correction
  void SetOMdbase(AliObjMatrix* omdb);                // Set the OM dbase object
  void SetCalibFile(TString name);                    // Set ROOT calibration input file
  void SetMinProb(Float_t pmin);                      // Set minimal probability to induce cross talk
  void SetXtalkPE(Float_t pe);                        // Set nominal Xtalk signal in photo-electrons

 protected :
  TFile* fCalfile;     // The (optional) calibration input file in ROOT format
  AliObjMatrix* fOmdb; // The OM database object
  Float_t fPmin;       // The minimal probability to induce cross talk 
  Float_t fPe;         // The nominal Xtalk signal in photo-electron equivalent

 ClassDef(IceXtalk,2) // TTask derived class to perform cross talk hit correction
};
#endif
