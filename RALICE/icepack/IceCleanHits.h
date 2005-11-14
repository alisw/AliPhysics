#ifndef IceCleanHits_h
#define IceCleanHits_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TROOT.h"
#include "TTask.h"
#include "TString.h"

#include "AliJob.h"
#include "IceEvent.h"
#include "IceAOM.h"
#include "IceIDOM.h"
#include "IceTDOM.h"

class IceCleanHits : public TTask
{
 public :
  IceCleanHits(const char* name="",const char* title=""); // Constructor
  virtual ~IceCleanHits();                                // Destructor
  virtual void Exec(Option_t* opt);                       // Hit cleaning
  void SetAdcRangeA(Float_t min,Float_t max=999999);      // Set Amanda ADC range (in PE)
  void SetTotRangeAE(Float_t min,Float_t max=2000);       // Set Amanda electrical TOT range (in ns)
  void SetTotRangeAO(Float_t min,Float_t max=2000);       // Set Amanda optical TOT range (in ns)
  void SetIsolationA(Float_t rmax,Float_t dtmax);         // Set Amanda isolation radius (in m) and dt (in ns)
  void SetTwindowA(Float_t dtmax);                        // Set Amanda maximal trigger window (in TDC counts)
  void SetTtimeA(Float_t t);                              // Set Amanda trigger time (in TDC counts)
  void SetTnameA(TString name);                           // Set Amanda trigger name

 protected :
  IceEvent* fEvt;    // Pointer to the current event structure
  Float_t fAdcminA;  // Minimum Amanda ADC value in PE
  Float_t fAdcmaxA;  // Maximum Amanda ADC value in PE
  Float_t fTotminAE; // Minimum Amanda electrical TOT value in ns
  Float_t fTotmaxAE; // Maximum Amanda electrical TOT value in ns
  Float_t fTotminAO; // Minimum Amanda optical TOT value in ns
  Float_t fTotmaxAO; // Maximum Amanda optical TOT value in ns
  Float_t fRmaxA;    // Maximum Amanda isolation radius in m
  Float_t fDtmaxA;   // Maximum Amanda isolation dt in ns
  Float_t fTwinA;    // Maximum Amanda hit time difference from the trigger time in TDC counts
  Float_t fTtimA;    // The Amanda trigger time in TDC counts
  TString fTnamA;    // The Amanda trigger name
  void Amanda();     // Cleaning of Amanda modules
  void InIce();      // Cleaning of IceCube InIce DOMs
  void IceTop();     // Cleaning of IceTop DOMs

 ClassDef(IceCleanHits,3) // TTask derived class to perform hit cleaning
};
#endif
