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
// The name of the configuration file can be passed as argument to the       //
// AliSimulation constructor or can be specified by                          //
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
// Background events can be merged by calling                                //
//                                                                           //
//   sim.MergeWith("background/galice.root", 2);                             //
//                                                                           //
// The first argument is the file name of the background galice file. The    //
// second argument is the number of signal events per background event.      //
// By default this number is calculated from the number of available         //
// background events. MergeWith can be called several times to merge more    //
// than two event streams. It is assumed that the sdigits were already       //
// produced for the background events.                                       //
//                                                                           //
// The output of raw data can be switched on by calling                      //
//                                                                           //
//   sim.SetWriteRawData("MUON");   // write raw data for MUON               //
//                                                                           //
// The methods RunSimulation, RunSDigitization, RunDigitization,             //
// RunHitsDigitization and WriteRawData can be used to run only parts of     //
// the full simulation chain.                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObjString.h>
#include <TStopwatch.h>
#include <TSystem.h>

#include "AliDigitizer.h"
#include "AliGenerator.h"
#include "AliModule.h"
#include "AliRun.h"
#include "AliRunDigitizer.h"
#include "AliRunLoader.h"
#include "AliSimulation.h"
#include "AliVertexGenFile.h"

ClassImp(AliSimulation)


//_____________________________________________________________________________
AliSimulation::AliSimulation(const char* configFileName,
			     const char* name, const char* title) :
  TNamed(name, title),

  fRunGeneration(kTRUE),
  fRunSimulation(kTRUE),
  fMakeSDigits("ALL"),
  fMakeDigits("ALL"),
  fMakeDigitsFromHits(""),
  fWriteRawData(""),
  fStopOnError(kFALSE),

  fNEvents(1),
  fConfigFileName(configFileName),
  fGAliceFileName("galice.root"),
  fBkgrdFileNames(NULL),
  fUseBkgrdVertex(kTRUE),
  fRegionOfInterest(kTRUE)
{
// create simulation object with default parameters

  SetGAliceFile("galice.root");
}

//_____________________________________________________________________________
AliSimulation::AliSimulation(const AliSimulation& sim) :
  TNamed(sim),

  fRunGeneration(sim.fRunGeneration),
  fRunSimulation(sim.fRunSimulation),
  fMakeSDigits(sim.fMakeSDigits),
  fMakeDigits(sim.fMakeDigits),
  fMakeDigitsFromHits(sim.fMakeDigitsFromHits),
  fWriteRawData(sim.fWriteRawData),
  fStopOnError(sim.fStopOnError),

  fNEvents(sim.fNEvents),
  fConfigFileName(sim.fConfigFileName),
  fGAliceFileName(sim.fGAliceFileName),
  fBkgrdFileNames(NULL),
  fUseBkgrdVertex(sim.fUseBkgrdVertex),
  fRegionOfInterest(sim.fRegionOfInterest)
{
// copy constructor

  fBkgrdFileNames = new TObjArray;
  for (Int_t i = 0; i < sim.fBkgrdFileNames->GetEntriesFast(); i++) {
    if (!sim.fBkgrdFileNames->At(i)) continue;
    fBkgrdFileNames->Add(sim.fBkgrdFileNames->At(i)->Clone());
  }
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

  if (fBkgrdFileNames) {
    fBkgrdFileNames->Delete();
    delete fBkgrdFileNames;
  }
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
void AliSimulation::SetGAliceFile(const char* fileName)
{
// set the name of the galice file
// the path is converted to an absolute one if it is relative

  fGAliceFileName = fileName;
  if (!gSystem->IsAbsoluteFileName(fGAliceFileName)) {
    char* absFileName = gSystem->ConcatFileName(gSystem->WorkingDirectory(),
						fGAliceFileName);
    fGAliceFileName = absFileName;
    delete[] absFileName;
  }
}

//_____________________________________________________________________________
void AliSimulation::MergeWith(const char* fileName, Int_t nSignalPerBkgrd)
{
// add a file with background events for merging

  TObjString* fileNameStr = new TObjString(fileName);
  fileNameStr->SetUniqueID(nSignalPerBkgrd);
  if (!fBkgrdFileNames) fBkgrdFileNames = new TObjArray;
  fBkgrdFileNames->Add(fileNameStr);
}


//_____________________________________________________________________________
Bool_t AliSimulation::Run(Int_t nEvents)
{
// run the generation, simulation and digitization

  if (nEvents > 0) fNEvents = nEvents;

  // generation and simulation -> hits
  if (fRunGeneration) {
    if (!RunSimulation()) if (fStopOnError) return kFALSE;
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
    if (fBkgrdFileNames && (fBkgrdFileNames->GetEntriesFast() > 0)) {
      Warning("Run", "Merging and direct creation of digits from hits " 
	      "was selected for some detectors. "
	      "No merging will be done for the following detectors: %s",
	      fMakeDigitsFromHits.Data());
    }
    if (!RunHitsDigitization(fMakeDigitsFromHits)) {
      if (fStopOnError) return kFALSE;
    }
  }

  // digits -> raw data
  if (!fWriteRawData.IsNull()) {
    if (!WriteRawData(fWriteRawData)) {
      if (fStopOnError) return kFALSE;
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::RunSimulation(Int_t nEvents)
{
// run the generation and simulation

  TStopwatch stopwatch;
  stopwatch.Start();

  if (!gAlice) {
    Error("RunSimulation", "no gAlice object. Restart aliroot and try again.");
    return kFALSE;
  }
  if (gAlice->Modules()->GetEntries() > 0) {
    Error("RunSimulation", 
	  "gAlice was already run. Restart aliroot and try again.");
    return kFALSE;
  }

  Info("RunSimulation", "initializing gAlice with config file %s",
       fConfigFileName.Data());
  gAlice->Init(fConfigFileName.Data());
  AliRunLoader* runLoader = gAlice->GetRunLoader();
  if (!runLoader) {
    Error("RunSimulation", "gAlice has no run loader object. "
	  "Check your config file: %s", fConfigFileName.Data());
    return kFALSE;
  }
  SetGAliceFile(runLoader->GetFileName());

  if (!gAlice->Generator()) {
    Error("RunSimulation", "gAlice has no generator object. "
	  "Check your config file: %s", fConfigFileName.Data());
    return kFALSE;
  }
  if (nEvents <= 0) nEvents = fNEvents;

  // get vertex from background file in case of merging
  if (fUseBkgrdVertex &&
      fBkgrdFileNames && (fBkgrdFileNames->GetEntriesFast() > 0)) {
    Int_t signalPerBkgrd = GetNSignalPerBkgrd(nEvents);
    const char* fileName = ((TObjString*)
			    (fBkgrdFileNames->At(0)))->GetName();
    Info("RunSimulation", "The vertex will be taken from the background "
	 "file %s with nSignalPerBackground = %d", 
	 fileName, signalPerBkgrd);
    AliVertexGenFile* vtxGen = new AliVertexGenFile(fileName, signalPerBkgrd);
    gAlice->Generator()->SetVertexGenerator(vtxGen);
  }

  if (!fRunSimulation) {
    gAlice->Generator()->SetTrackingFlag(0);
  }

  Info("RunSimulation", "running gAlice");
  gAlice->Run(nEvents);

  delete runLoader;

  Info("RunSimulation", "execution time:");
  stopwatch.Print();

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::RunSDigitization(const char* detectors)
{
// run the digitization and produce summable digits

  TStopwatch stopwatch;
  stopwatch.Start();

  AliRunLoader* runLoader = LoadRun();
  if (!runLoader) return kFALSE;

  TString detStr = detectors;
  TObjArray* detArray = runLoader->GetAliRun()->Detectors();
  for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
    AliModule* det = (AliModule*) detArray->At(iDet);
    if (!det || !det->IsActive()) continue;
    if (IsSelected(det->GetName(), detStr)) {
      Info("RunSDigitization", "creating summable digits for %s", 
	   det->GetName());
      TStopwatch stopwatchDet;
      stopwatchDet.Start();
      det->Hits2SDigits();
      Info("RunSDigitization", "execution time for %s:", det->GetName());
      stopwatchDet.Print();
    }
  }

  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    Error("RunSDigitization", "the following detectors were not found: %s", 
	  detStr.Data());
    if (fStopOnError) return kFALSE;
  }

  delete runLoader;

  Info("RunSDigitization", "execution time:");
  stopwatch.Print();

  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliSimulation::RunDigitization(const char* detectors, 
				      const char* excludeDetectors)
{
// run the digitization and produce digits from sdigits

  TStopwatch stopwatch;
  stopwatch.Start();

  while (AliRunLoader::GetRunLoader()) delete AliRunLoader::GetRunLoader();
  if (gAlice) delete gAlice;
  gAlice = NULL;

  Int_t nStreams = 1;
  if (fBkgrdFileNames) nStreams = fBkgrdFileNames->GetEntriesFast() + 1;
  Int_t signalPerBkgrd = GetNSignalPerBkgrd();
  AliRunDigitizer* manager = new AliRunDigitizer(nStreams, signalPerBkgrd);
  manager->SetInputStream(0, fGAliceFileName.Data());
  for (Int_t iStream = 1; iStream < nStreams; iStream++) {
    const char* fileName = ((TObjString*)
			    (fBkgrdFileNames->At(iStream-1)))->GetName();
    manager->SetInputStream(iStream, fileName);
  }

  TString detStr = detectors;
  TString detExcl = excludeDetectors;
  manager->GetInputStream(0)->ImportgAlice();
  AliRunLoader* runLoader = 
    AliRunLoader::GetRunLoader(manager->GetInputStream(0)->GetFolderName());
  TObjArray* detArray = runLoader->GetAliRun()->Detectors();
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
Bool_t AliSimulation::RunHitsDigitization(const char* detectors)
{
// run the digitization and produce digits from hits

  TStopwatch stopwatch;
  stopwatch.Start();

  AliRunLoader* runLoader = LoadRun();
  if (!runLoader) return kFALSE;

  TString detStr = detectors;
  TObjArray* detArray = runLoader->GetAliRun()->Detectors();
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

  delete runLoader;

  Info("RunHitsDigitization", "execution time:");
  stopwatch.Print();

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliSimulation::WriteRawData(const char* detectors)
{
// convert the digits to raw data

  TStopwatch stopwatch;
  stopwatch.Start();

  AliRunLoader* runLoader = LoadRun();
  if (!runLoader) return kFALSE;

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    Info("WriteRawData", "processing event %d", iEvent);
    runLoader->GetEvent(iEvent);
    TString baseDir = gSystem->WorkingDirectory();
    char dirName[256];
    sprintf(dirName, "raw%d", iEvent);
    gSystem->MakeDirectory(dirName);
    if (!gSystem->ChangeDirectory(dirName)) {
      Error("WriteRawData", "couldn't change to directory %s", dirName);
      if (fStopOnError) return kFALSE; else continue;
    }

    TString detStr = detectors;
    TObjArray* detArray = runLoader->GetAliRun()->Detectors();
    for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
      AliModule* det = (AliModule*) detArray->At(iDet);
      if (!det || !det->IsActive()) continue;
      if (IsSelected(det->GetName(), detStr)) {
	Info("WriteRawData", "creating raw data from digits for %s", 
	     det->GetName());
	det->Digits2Raw();
      }
    }

    gSystem->ChangeDirectory(baseDir);
    if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
      Error("WriteRawData", "the following detectors were not found: %s", 
	    detStr.Data());
      if (fStopOnError) return kFALSE;
    }
  }

  delete runLoader;

  Info("WriteRawData", "execution time:");
  stopwatch.Print();

  return kTRUE;
}


//_____________________________________________________________________________
AliRunLoader* AliSimulation::LoadRun() const
{
// delete existing run loaders, open a new one and load gAlice

  while (AliRunLoader::GetRunLoader()) delete AliRunLoader::GetRunLoader();
  AliRunLoader* runLoader = 
    AliRunLoader::Open(fGAliceFileName.Data(), 
		       AliConfig::GetDefaultEventFolderName(), "UPDATE");
  if (!runLoader) {
    Error("LoadRun", "no run loader found in file %s", 
	  fGAliceFileName.Data());
    return NULL;
  }
  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  if (!gAlice) {
    Error("LoadRun", "no gAlice object found in file %s", 
	  fGAliceFileName.Data());
    return NULL;
  }
  return runLoader;
}

//_____________________________________________________________________________
Int_t AliSimulation::GetNSignalPerBkgrd(Int_t nEvents) const
{
// get or calculate the number of signal events per background event

  if (!fBkgrdFileNames) return 1;
  Int_t nBkgrdFiles = fBkgrdFileNames->GetEntriesFast();
  if (nBkgrdFiles == 0) return 1;

  // get the number of signal events
  if (nEvents <= 0) {
    AliRunLoader* runLoader = 
      AliRunLoader::Open(fGAliceFileName.Data(), "SIGNAL");
    if (!runLoader) return 1;
    nEvents = runLoader->GetNumberOfEvents();
    delete runLoader;
  }

  Int_t result = 0;
  for (Int_t iBkgrdFile = 0; iBkgrdFile < nBkgrdFiles; iBkgrdFile++) {
    // get the number of background events
    const char* fileName = ((TObjString*)
			    (fBkgrdFileNames->At(iBkgrdFile)))->GetName();
    AliRunLoader* runLoader = 
      AliRunLoader::Open(fileName, "BKGRD");
    if (!runLoader) continue;
    Int_t nBkgrdEvents = runLoader->GetNumberOfEvents();
    delete runLoader;

    // get or calculate the number of signal per background events
    Int_t nSignalPerBkgrd = fBkgrdFileNames->At(iBkgrdFile)->GetUniqueID();
    if (nSignalPerBkgrd <= 0) {
      nSignalPerBkgrd = (nEvents-1) / nBkgrdEvents + 1;
    } else if (result && (result != nSignalPerBkgrd)) {
      Info("GetNSignalPerBkgrd", "the number of signal events per "
	   "background event will be changed from %d to %d for stream %d", 
	   nSignalPerBkgrd, result, iBkgrdFile+1);
      nSignalPerBkgrd = result;
    }

    if (!result) result = nSignalPerBkgrd;
    if (nSignalPerBkgrd * nBkgrdEvents < nEvents) {
      Warning("GetNSignalPerBkgrd", "not enough background events (%d) for "
	      "%d signal events using %d signal per background events for "
	      "stream %d", 
	      nBkgrdEvents, nEvents, nSignalPerBkgrd, iBkgrdFile+1);
    }
  }

  return result;
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
