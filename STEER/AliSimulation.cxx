/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for running generation, simulation and digitization                 //
//                                                                           //
// Hits, sdigits and digits are created for all detectors by typing:         //
//                                                                           //
//   AliSimulation sim;                                                      //
//   sim.Run();                                                              //
//                                                                           //
// The Run method returns kTRUE in case of successful execution.             //
// The number of events can be given as argument to the Run method or it     //
// can be set by                                                             //
//                                                                           //
//   sim.SetNumberOfEvents(n);                                               //
//                                                                           //
// The name of the configuration file can be specified by                    //
//                                                                           //
//   sim.SetConfigFile("...");                                               //
//                                                                           //
// The generation of particles and the simulation of detector hits can be    //
// switched on or off by                                                     //
//                                                                           //
//   sim.SetRunGeneration(kTRUE);   // generation of primary particles       //
//   sim.SetRunSimulation(kFALSE);  // but no tracking                       //
//                                                                           //
// For which detectors sdigits and digits will be created, can be steered    //
// by                                                                        //
//                                                                           //
//   sim.SetMakeSDigits("ALL");     // make sdigits for all detectors        //
//   sim.SetMakeDigits("ITS TPC");  // make digits only for ITS and TPC      //
//                                                                           //
// The argument is a (case sensitive) string with the names of the           //
// detectors separated by a space. An empty string ("") can be used to       //
// disable the creation of sdigits or digits. The special string "ALL"       //
// selects all available detectors. This is the default.                     //
//                                                                           //
// The creation of digits from hits instead of from sdigits can be selected  //
// by                                                                        //
//                                                                           //
//   sim.SetMakeDigitsFromHits("TRD");                                       //
//                                                                           //
// The argument is again a string with the selected detectors. Be aware that //
// this feature is not available for all detectors and that merging is not   //
// possible, when digits are created directly from hits.                     //
//                                                                           //
// Backgound events can be merged by calling                                 //
//                                                                           //
//   sim.MergeWith("background/galice.root", 2);                             //
//                                                                           //
// The first argument is the file name of the background galice file. The    //
// second argument is the number of signal events per background event.      //
// The default value for this is 1. MergeWith can be called several times    //
// to merge more than two event streams. It is assumed that the sdigits      //
// were already produced for the background events.                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliSimulation.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliModule.h"
#include "AliGenerator.h"
#include "AliRunDigitizer.h"
#include "AliDigitizer.h"
#include <TObjString.h>


ClassImp(AliSimulation)


//_____________________________________________________________________________
AliSimulation::AliSimulation(const char* name, const char* title) :
  TNamed(name, title)
{
// create simulation object with default parameters

  Init();
}

//_____________________________________________________________________________
AliSimulation::AliSimulation(const AliSimulation& sim) :
  TNamed(sim)
{
// copy constructor

  fRunGeneration = sim.fRunGeneration;
  fRunSimulation = sim.fRunSimulation;
  fMakeSDigits = sim.fMakeSDigits;
  fMakeDigits = sim.fMakeDigits;
  fMakeDigitsFromHits = sim.fMakeDigitsFromHits;
  fStopOnError = sim.fStopOnError;

  fNEvents = sim.fNEvents;
  fConfigFileName = sim.fConfigFileName;
  fGAliceFileName = sim.fGAliceFileName;
  fBkgrdFileNames = new TObjArray;
  for (Int_t i = 0; i < sim.fBkgrdFileNames->GetEntriesFast(); i++) {
    if (!sim.fBkgrdFileNames->At(i)) continue;
    fBkgrdFileNames->Add(sim.fBkgrdFileNames->At(i)->Clone());
  }

  fRunLoader = NULL;
}

//_____________________________________________________________________________
AliSimulation& AliSimulation::operator = (const AliSimulation& sim)
{
// assignment operator

  this->~AliSimulation();
  new(this) AliSimulation(sim);
  return *this;
}

//_____________________________________________________________________________
AliSimulation::~AliSimulation()
{
// clean up

  fBkgrdFileNames->Delete();
  delete fBkgrdFileNames;
}

//_____________________________________________________________________________
void AliSimulation::Init()
{
// set default parameters

  fRunGeneration = kTRUE;
  fRunSimulation = kTRUE;
  fMakeSDigits = "ALL";
  fMakeDigits = "ALL";
  fMakeDigitsFromHits = "";
  fStopOnError = kFALSE;

  fNEvents = 1;
  fConfigFileName = "Config.C";
  fGAliceFileName = "galice.root";
  fBkgrdFileNames = new TObjArray;
  fRegionOfInterest = kTRUE;

  fRunLoader = NULL;
}


//_____________________________________________________________________________
void AliSimulation::SetNumberOfEvents(Int_t nEvents)
{
// set the number of events for one run

  fNEvents = nEvents;
}

//_____________________________________________________________________________
void AliSimulation::SetConfigFile(const char* fileName)
{
// set the name of the config file

  fConfigFileName = fileName;
}

//_____________________________________________________________________________
void AliSimulation::MergeWith(const char* fileName, Int_t nSignalPerBkgrd)
{
// add a file with background events for merging

  TObjString* fileNameStr = new TObjString(fileName);
  fileNameStr->SetUniqueID(nSignalPerBkgrd);
  fBkgrdFileNames->Add(fileNameStr);
}


//_____________________________________________________________________________
Bool_t AliSimulation::Run(Int_t nEvents)
{
// run the generation, simulation and digitization

  if (nEvents > 0) fNEvents = nEvents;

  // generation and simulation -> hits
  if (fRunGeneration) {
    if (!gAlice) {
      Error("Run", "no gAlice object. Restart aliroot and try again.");
      return kFALSE;
    }
    if (gAlice->Modules()->GetEntries() > 0) {
      Error("Run", "gAlice was already run. Restart aliroot and try again.");
      return kFALSE;
    }
    if (!RunSimulation()) if (fStopOnError) return kFALSE;
  }

  // reopen the run loader
  if (fRunLoader) delete fRunLoader;
  fRunLoader = AliRunLoader::Open(fGAliceFileName.Data(),AliConfig::fgkDefaultEventFolderName,"UPDATE");
  if (!fRunLoader) {
    Error("Run", "no run loader found in file %s", 
	  fGAliceFileName.Data());
    return kFALSE;
  }
  fRunLoader->LoadgAlice();
  gAlice = fRunLoader->GetAliRun();
  if (!gAlice) {
    Error("Run", "no gAlice object found in file %s", 
	  fGAliceFileName.Data());
    return kFALSE;
  }

  // hits -> summable digits
  if (!fMakeSDigits.IsNull()) {
    if (!RunSDigitization(fMakeSDigits)) if (fStopOnError) return kFALSE;
  }

  // summable digits -> digits
  if (!fMakeDigits.IsNull()) {
    if (!RunDigitization(fMakeDigits, fMakeDigitsFromHits)) {
      if (fStopOnError) return kFALSE;
    }
  }

  // hits -> digits
  if (!fMakeDigitsFromHits.IsNull()) {
    if (fBkgrdFileNames->GetEntriesFast() > 0) {
      Warning("Run", "Merging and direct creation of digits from hits " 
	      "was selected for some detectors. "
	      "No merging will be done for the following detectors: %s",
	      fMakeDigitsFromHits.Data());
    }
    if (!RunHitsDigitization(fMakeDigitsFromHits)) {
      if (fStopOnError) return kFALSE;
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::RunSimulation()
{
// run the generation and simulation

  TStopwatch stopwatch;
  stopwatch.Start();

  Info("RunSimulation", "initializing gAlice with config file %s",
       fConfigFileName.Data());
  gAlice->Init(fConfigFileName.Data());
  fRunLoader = gAlice->GetRunLoader();
  if (!fRunLoader) {
    Error("RunSimulation", "gAlice has no run loader object. "
	  "Check your config file: %s", fConfigFileName.Data());
    return kFALSE;
  }
  fGAliceFileName = fRunLoader->GetFileName();

  if (!fRunSimulation) {
    if (!gAlice->Generator()) {
      Error("RunSimulation", "gAlice has no generator object. "
	    "Check your config file: %s", fConfigFileName.Data());
      return kFALSE;
    }
    gAlice->Generator()->SetTrackingFlag(0);
  }

  Info("RunSimulation", "running gAlice");
  gAlice->Run(fNEvents);

  Info("RunSimulation", "execution time:");
  stopwatch.Print();

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::RunSDigitization(const TString& detectors)
{
// run the digitization and produce summable digits

  TStopwatch stopwatch;
  stopwatch.Start();

  TString detStr = detectors;
  TObjArray* detArray = gAlice->Detectors();
  for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
    AliModule* det = (AliModule*) detArray->At(iDet);
    if (!det || !det->IsActive()) continue;
    if (IsSelected(det->GetName(), detStr)) {
      Info("RunSDigitization", "creating summable digits for %s", 
	   det->GetName());
      det->Hits2SDigits();
    }
  }

  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    Error("RunSDigitization", "the following detectors were not found: %s", 
	  detStr.Data());
    if (fStopOnError) return kFALSE;
  }

  Info("RunSDigitization", "execution time:");
  stopwatch.Print();

  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliSimulation::RunDigitization(const TString& detectors, 
				      const TString& excludeDetectors)
{
// run the digitization and produce digits from sdigits

  TStopwatch stopwatch;
  stopwatch.Start();

  Int_t nStreams = fBkgrdFileNames->GetEntriesFast() + 1;
  Int_t signalPerBkgrd = 1;
  if (nStreams > 1) signalPerBkgrd = fBkgrdFileNames->At(0)->GetUniqueID();
  AliRunDigitizer* manager = new AliRunDigitizer(nStreams, signalPerBkgrd);
  manager->SetInputStream(0, fGAliceFileName.Data());
  for (Int_t iStream = 1; iStream < nStreams; iStream++) {
    const char* fileName = ((TObjString*)
			    (fBkgrdFileNames->At(iStream-1)))->GetName();
    manager->SetInputStream(iStream, fileName);
  }

  TString detStr = detectors;
  TString detExcl = excludeDetectors;
  TObjArray* detArray = gAlice->Detectors();
  for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
    AliModule* det = (AliModule*) detArray->At(iDet);
    if (!det || !det->IsActive()) continue;
    if (IsSelected(det->GetName(), detStr) && 
	!IsSelected(det->GetName(), detExcl)) {
      AliDigitizer* digitizer = det->CreateDigitizer(manager);
      if (!digitizer) {
	Error("RunDigitization", "no digitizer for %s", det->GetName());
	if (fStopOnError) return kFALSE;
      } else {
	digitizer->SetRegionOfInterest(fRegionOfInterest);
      }
    }
  }

  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    Error("RunDigitization", "the following detectors were not found: %s", 
	  detStr.Data());
    if (fStopOnError) return kFALSE;
  }

  if (!manager->GetListOfTasks()->IsEmpty()) {
    Info("RunDigitization", "executing digitization");
    manager->Exec("");
  }
  delete manager;

  Info("RunDigitization", "execution time:");
  stopwatch.Print();

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::RunHitsDigitization(const TString& detectors)
{
// run the digitization and produce digits from hits

  TStopwatch stopwatch;
  stopwatch.Start();

  TString detStr = detectors;
  TObjArray* detArray = gAlice->Detectors();
  for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
    AliModule* det = (AliModule*) detArray->At(iDet);
    if (!det || !det->IsActive()) continue;
    if (IsSelected(det->GetName(), detStr)) {
      Info("RunHitsDigitization", "creating digits from hits for %s", 
	   det->GetName());
      det->Hits2Digits();
    }
  }

  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    Error("RunHitsDigitization", "the following detectors were not found: %s", 
	  detStr.Data());
    if (fStopOnError) return kFALSE;
  }

  Info("RunHitsDigitization", "execution time:");
  stopwatch.Print();

  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliSimulation::IsSelected(TString detName, TString& detectors) const
{
// check whether detName is contained in detectors
// if yes, it is removed from detectors

  // check if all detectors are selected
  if ((detectors.CompareTo("ALL") == 0) ||
      detectors.BeginsWith("ALL ") ||
      detectors.EndsWith(" ALL") ||
      detectors.Contains(" ALL ")) {
    detectors = "ALL";
    return kTRUE;
  }

  // search for the given detector
  Bool_t result = kFALSE;
  if ((detectors.CompareTo(detName) == 0) ||
      detectors.BeginsWith(detName+" ") ||
      detectors.EndsWith(" "+detName) ||
      detectors.Contains(" "+detName+" ")) {
    detectors.ReplaceAll(detName, "");
    result = kTRUE;
  }

  // clean up the detectors string
  while (detectors.Contains("  ")) detectors.ReplaceAll("  ", " ");
  while (detectors.BeginsWith(" ")) detectors.Remove(0, 1);
  while (detectors.EndsWith(" ")) detectors.Remove(detectors.Length()-1, 1);

  return result;
}
