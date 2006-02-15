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

 ClassDef(IcePandel,1) // TTask derived class to perform Pandel fitting
};
#endif
