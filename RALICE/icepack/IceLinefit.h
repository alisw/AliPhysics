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
  IceLinefit(const char* name="IceLinefit",const char* title="Linefit reconstruction"); // Constructor
  virtual ~IceLinefit();                                // Destructor
  virtual void Exec(Option_t* opt);                     // Linefit reconstruction
  void SetMaxMod(Int_t nmax,TString s="A");  // Set max. number of good fired (D)OMs for events to be processed
  void SetMinMod(Int_t nmin,TString s="A");  // Set min. number of good fired (D)OMs for events to be processed
  void SetMaxHits(Int_t nmax,TString s="A"); // Set max. number of good hits per (D)OM to be processed
  void SetTrackName(TString s);  // Set (alternative) name for the produced first guess tracks
  void SetCharge(Float_t charge);// Set user defined charge for the produced first guess tracks

 protected :
  IceEvent* fEvt;    // Pointer to the current event
  Int_t fMaxmodA;    // The max. number of good fired Amanda OMs for events to be processed
  Int_t fMinmodA;    // The min. number of good fired Amanda OMs for events to be processed
  Int_t fMaxhitsA;   // The maximum number of good hits per Amanda OM to be processed
  Int_t fMaxmodI;    // The max. number of good fired InIce DOMs for events to be processed
  Int_t fMinmodI;    // The min. number of good fired InIce DOMs for events to be processed
  Int_t fMaxhitsI;   // The maximum number of good hits per InIce DOM to be processed
  TString fTrackname;// The name identifier for the produced first guess tracks
  Float_t fCharge;   // User defined charge of the produced first guess tracks
  virtual void Amanda(); // Linefit reconstruction for Amanda OMs
  virtual void InIce();  // Linefit reconstruction for InIce DOMs

 ClassDef(IceLinefit,5) // TTask derived class to perform linefit reconstruction
};
#endif
