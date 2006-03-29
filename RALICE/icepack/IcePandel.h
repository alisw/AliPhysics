#ifndef IcePandel_h
#define IcePandel_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TROOT.h"
#include "TTask.h"
#include "TString.h"
#include "TObjString.h"
#include "TFitter.h"

#include "AliJob.h"
#include "IceEvent.h"
#include "IceGOM.h"

class IcePandel : public TTask
{
 public :
  IcePandel(const char* name="",const char* title=""); // Constructor
  virtual ~IcePandel();                                // Destructor
  virtual void Exec(Option_t* opt);                    // Perform the fitting procedure
  void SetPrintLevel(Int_t level);                     // Set the fitter (Minuit) printlevel
  void UseTracks(TString classname,Int_t n=-1);        // Specify first guess tracks to be used
  void SelectHits(Int_t mode=1);                       // Specify which hits to be used
  void SetTrackName(TString s);  // Set (alternative) name for the produced tracks
  void SetCharge(Float_t charge);// Set user defined charge for the produced tracks
  void SetPenalty(Float_t val);  // Set penalty value in dB for the minimiser outside the allowed range
  void FitFCN(Int_t&,Double_t*,Double_t&,Double_t*,Int_t); // The minimisation FCN

 protected :
  Int_t fFirst;         // Flag to denote first invokation of the processor
  Int_t fPrint;         // Flag to denote the fitter (Minuit) printlevel
  Int_t fSelhits;       // Flag to denote which hits to be used
  IceEvent* fEvt;       // Pointer to the current event structure
  TObjArray* fUseNames; // The first guess classnames to be used 
  TArrayI* fUseNtk;     // The max. numbers of the various first guess tracks to be used
  AliTrack* fTrack;     // Pointer to the first guess track being processed
  TObjArray* fHits;     // The various hits to be used in the fitting process 
  TFitter* fFitter;     // Pointer to the minimisation processor
  TString fTrackname;   // The name identifier for the produced tracks
  Float_t fCharge;      // User defined charge of the produced tracks
  Float_t fPenalty;     // User defined penalty value in dB for the minimiser outside the allowed range

 ClassDef(IcePandel,4) // TTask derived class to perform Pandel fitting
};
#endif
