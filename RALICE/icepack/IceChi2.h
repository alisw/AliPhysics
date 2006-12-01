#ifndef IceChi2_h
#define IceChi2_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TROOT.h"
#include "TTask.h"
#include "TString.h"
#include "TObjString.h"
#include "TFitter.h"
#include "Math/SpecFuncMathMore.h"

#include "AliJob.h"
#include "AliSample.h"
#include "IceEvent.h"
#include "IceGOM.h"

class IceChi2 : public TTask
{
 public :
  IceChi2(const char* name="",const char* title=""); // Constructor
  virtual ~IceChi2();                                // Destructor
  virtual void Exec(Option_t* opt);                  // Perform the fitting procedure
  void SetPrintLevel(Int_t level);                   // Set the fitter (Minuit) printlevel
  void UseTracks(TString classname,Int_t n=-1);      // Specify first guess tracks to be used
  void SelectHits(Int_t mode=1);                     // Specify which hits to be used
  void SetVgroupUsage(Int_t flag); // (De)activate usage of distinct phase and group velocities
  void SetTrackName(TString s);  // Set (alternative) name for the produced tracks
  void SetCharge(Float_t charge);// Set user defined charge for the produced tracks
  void SetPenalty(Float_t val);  // Set psi penalty value (dB) for extreme distance and/or time values
  Double_t GetPsi(AliTrack* t);  // Provide Bayesian psi value for a track w.r.t. CPandel PDF
  void FitFCN(Int_t&,Double_t*,Double_t&,Double_t*,Int_t); // The minimisation FCN

 protected :
  Int_t fFirst;         // Flag to denote first invokation of the processor
  Int_t fPrint;         // Flag to denote the fitter (Minuit) printlevel
  Int_t fSelhits;       // Flag to denote which hits to be used
  Int_t fVgroup;        // Flag to indicate usage of distinct phase and group velocities
  IceEvent* fEvt;       // Pointer to the current event structure
  TObjArray* fUseNames; // The first guess classnames to be used 
  TArrayI* fUseNtk;     // The max. numbers of the various first guess tracks to be used
  TObjArray* fHits;     // The various hits to be used in the fitting process 
  TFitter* fFitter;     // Pointer to the minimisation processor
  TString fTrackname;   // The name identifier for the produced tracks
  Float_t fCharge;      // User defined charge of the produced tracks
  Float_t fPenalty;     // User defined psi penalty value (dB) for extreme distance and/or time values
  AliTrack* fTkfit;     // Pointer to the produced fitted track
  AliSignal* fFitstats; // The fit details of the produced fitted track
  AliSample fPsistats;  // Statistics of the Bayesian psi value for the best fitted track

 ClassDef(IceChi2,2) // TTask derived class to perform Chi-squared track fitting
};
#endif
