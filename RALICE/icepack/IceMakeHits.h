#ifndef IceMakeHits_h
#define IceMakeHits_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TROOT.h"
#include "TTask.h"
#include "TString.h"
#include "TSpectrum.h"
#include "TMath.h"

#include "AliJob.h"
#include "IceEvent.h"
#include "IceAOM.h"
#include "IceIDOM.h"
#include "IceTDOM.h"

class IceMakeHits : public TTask
{
 public :
  IceMakeHits(const char* name="",const char* title=""); // Constructor
  virtual ~IceMakeHits();                                // Destructor
  virtual void Exec(Option_t* opt);                      // Hit extraction
  void SetBasefracA(Float_t val);                        // Set the Amanda fract. baseline update value
  void SetSigmaA(Float_t val);                           // Set the Amanda clipping window width
  void SetMaxPeaksA(Int_t val);                          // Set the Amanda maximum number of peaks in a waveform
  void SetMinPulseHeightA(Float_t val);                  // Set the Amanda threshold for narrow pulses
  void SetThresholdA(Float_t val);                       // Set the Amanda threshold for narrow pulses

 protected :
  IceEvent* fEvt;            // Pointer to the current event structure
  AliDevice* fDaq;           // The DAQ info for the current event
  Float_t fBasefracA;        // The fractional baseline update for Amanda TWR extraction 
  Float_t fSigmaA;           // The width of the clipping window to be used by TSpectrum::Search() in Amanda TWR extraction
  Int_t fMaxPeaksA;          // The maximum number of peaks in a waveform in Amanda TWR extraction
  Float_t fMinPulseHeightA;  // The minimum pulse height for narrow pulses in Amanda TWR extraction
  Float_t fThresholdA;       // The threshold to be used in analysis of narrow pulses in Amanda TWR extraction
  void Amanda();             // Hit extraction from Amanda TWR data
  void InIce();              // Hit extraction from IceCube InIce ATWD data
  void IceTop();             // Hit extraction from IceTop ATWD data

 ClassDef(IceMakeHits,1) // TTask derived class to perform hit extraction from waveforms
};
#endif
