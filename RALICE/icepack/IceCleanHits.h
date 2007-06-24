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
  IceCleanHits(const char* name="IceCleanHits",const char* title="Hit cleaning"); // Constructor
  virtual ~IceCleanHits();                                // Destructor
  virtual void Exec(Option_t* opt);                       // Hit cleaning
  void SetAdcRangeA(Float_t min,Float_t max=999999,TString s="MuDaq");// Set Amanda ADC range (in PE)
  void SetTotRangeAE(Float_t min,Float_t max=2000,TString s="MuDaq"); // Set Amanda electrical TOT range (in ns)
  void SetTotRangeAO(Float_t min,Float_t max=2000,TString s="MuDaq"); // Set Amanda optical TOT range (in ns)
  void SetIsolationA(Float_t rmax,Float_t dtmax);         // Set Amanda isolation radius (in m) and dt (in ns)
  void SetTwindowA(Float_t dtmax,TString s="MuDaq");      // Set Amanda maximal trigger window
  void SetTtimeA(Float_t t,TString s="MuDaq");            // Set Amanda trigger time
  void SetTnameA(TString name,TString s="MuDaq");         // Set Amanda trigger name

 protected :
  IceEvent* fEvt;     // Pointer to the current event structure
  Float_t fAdcminAM;  // Minimum Amanda MuDaq ADC value in PE
  Float_t fAdcmaxAM;  // Maximum Amanda MuDaq ADC value in PE
  Float_t fTotminAEM; // Minimum Amanda MuDaq electrical TOT value in ns
  Float_t fTotmaxAEM; // Maximum Amanda MuDaq electrical TOT value in ns
  Float_t fTotminAOM; // Minimum Amanda MuDaq optical TOT value in ns
  Float_t fTotmaxAOM; // Maximum Amanda MuDaq optical TOT value in ns
  Float_t fAdcminAT;  // Minimum Amanda TWRDaq ADC value in PE
  Float_t fAdcmaxAT;  // Maximum Amanda TWRDaq ADC value in PE
  Float_t fTotminAET; // Minimum Amanda TWRDaq electrical TOT value in ns
  Float_t fTotmaxAET; // Maximum Amanda TWRDaq electrical TOT value in ns
  Float_t fTotminAOT; // Minimum Amanda TWRDaq optical TOT value in ns
  Float_t fTotmaxAOT; // Maximum Amanda TWRDaq optical TOT value in ns
  Float_t fRmaxA;     // Maximum Amanda isolation radius in m
  Float_t fDtmaxA;    // Maximum Amanda isolation dt in ns
  Float_t fTwinAM;    // Maximum Amanda MuDaq hit time difference from the trigger time in TDC counts
  Float_t fTwinAT;    // Maximum Amanda TWRDaq hit time difference from the trigger time in ns
  Float_t fTtimAM;    // The Amanda MuDaq trigger time in TDC counts
  TString fTnamAM;    // The Amanda MuDaq trigger name
  Float_t fTtimAT;    // The Amanda TWRDaq trigger time in ns
  TString fTnamAT;    // The Amanda TWRDaq trigger name
  void Amanda();      // Cleaning of Amanda modules
  void MuDaq();       // Cleaning of Amanda modules for MuDaq data
  void TWRDaq();      // Cleaning of Amanda modules for TWRDaq data
  void InIce();       // Cleaning of IceCube InIce DOMs
  void IceTop();      // Cleaning of IceTop DOMs

 ClassDef(IceCleanHits,4) // TTask derived class to perform hit cleaning
};
#endif
