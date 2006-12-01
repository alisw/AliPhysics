#ifndef IceLinefit_h
#define IceLinefit_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TROOT.h"
#include "TTask.h"
#include "TString.h"

#include "AliJob.h"
#include "IceEvent.h"
#include "IceGOM.h"

class IceLinefit : public TTask
{
 public :
  IceLinefit(const char* name="",const char* title=""); // Constructor
  virtual ~IceLinefit();                                // Destructor
  virtual void Exec(Option_t* opt);                     // Linefit reconstruction
  void SetMaxModA(Int_t nmax);   // Set max. number of good fired Amanda modules for events to be processed
  void SetMinModA(Int_t nmin);   // Set min. number of good fired Amanda modules for events to be processed
  void SetMaxHitsA(Int_t nmax);  // Set max. number of good hits per Amanda module to be processed
  void SetTrackName(TString s);  // Set (alternative) name for the produced first guess tracks
  void SetCharge(Float_t charge);// Set user defined charge for the produced first guess tracks

 protected :
  Int_t fMaxmodA;    // The max. number of good fired Amanda modules for events to be processed
  Int_t fMinmodA;    // The min. number of good fired Amanda modules for events to be processed
  Int_t fMaxhitsA;   // The maximum number of good hits per Amanda module to be processed
  TString fTrackname;// The name identifier for the produced first guess tracks
  Float_t fCharge;   // User defined charge of the produced first guess tracks

 ClassDef(IceLinefit,3) // TTask derived class to perform linefit reconstruction
};
#endif
