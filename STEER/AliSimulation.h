#ifndef ALISIMULATION_H
#define ALISIMULATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TNamed.h>
#include <TString.h>

class AliRunLoader;


class AliSimulation: public TNamed {
public:
  AliSimulation(const char* name = "AliSimulation", 
		const char* title = "generation, simulation and digitization");
  AliSimulation(const AliSimulation& sim);
  AliSimulation& operator = (const AliSimulation& sim);
  virtual ~AliSimulation();

  void           SetNumberOfEvents(Int_t nEvents);
  void           SetConfigFile(const char* fileName);

  void           SetRunGeneration(Bool_t run) {fRunGeneration = run;};
  void           SetRunSimulation(Bool_t run) {fRunSimulation = run;};
  void           SetMakeSDigits(const char* detectors) 
                   {fMakeSDigits = detectors;};
  void           MergeWith(const char* fileName, Int_t nSignalPerBkgrd = 1);
  void           SetRegionOfInterest(Bool_t flag) {fRegionOfInterest = flag;};
  void           SetMakeDigits(const char* detectors)
                   {fMakeDigits = detectors;};
  void           SetMakeDigitsFromHits(const char* detectors)
                   {fMakeDigitsFromHits = detectors;};

  virtual Bool_t Run(Int_t nEvents = 0);

private:
  void           Init();

  Bool_t         IsSelected(TString detName, TString& detectors) const;

  Bool_t         RunSimulation();
  Bool_t         RunSDigitization(const TString& detectors);
  Bool_t         RunDigitization(const TString& detectors,
				 const TString& excludeDetectors);
  Bool_t         RunHitsDigitization(const TString& detectors);

  Bool_t         fRunGeneration;      // generate prim. particles or not
  Bool_t         fRunSimulation;      // simulate detectors (hits) or not
  TString        fMakeSDigits;        // create sdigits for these detectors
  TString        fMakeDigits;         // create digits for these detectors
  TString        fMakeDigitsFromHits; // create digits from hits for these detectors
  Bool_t         fStopOnError;        // stop or continue on errors

  Int_t          fNEvents;            // number of events
  TString        fConfigFileName;     // name of the config file
  TString        fGAliceFileName;     // name of the galice file
  TObjArray*     fBkgrdFileNames;     // names of background files for merging
  Bool_t         fRegionOfInterest;   // digitization in region of interest

  AliRunLoader*  fRunLoader;          //! current run loader object

  ClassDef(AliSimulation, 1)  // class for running generation, simulation and digitization
};

#endif
