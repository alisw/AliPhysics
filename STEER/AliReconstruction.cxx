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
// class for running the reconstruction                                      //
//                                                                           //
// Clusters and tracks are created for all detectors and all events by       //
// typing:                                                                   //
//                                                                           //
//   AliReconstruction rec;                                                  //
//   rec.Run();                                                              //
//                                                                           //
// The Run method returns kTRUE in case of successful execution.             //
//                                                                           //
// If the input to the reconstruction are not simulated digits but raw data, //
// this can be specified by an argument of the Run method or by the method   //
//                                                                           //
//   rec.SetInput("...");                                                    //
//                                                                           //
// The input formats and the corresponding argument are:                     //
// - DDL raw data files: directory name, ends with "/"                       //
// - raw data root file: root file name, extension ".root"                   //
// - raw data DATE file: DATE file name, any other non-empty string          //
// - MC root files     : empty string, default                               //
//                                                                           //
// The name of the galice file can be changed from the default               //
// "galice.root" by passing it as argument to the AliReconstruction          //
// constructor or by                                                         //
//                                                                           //
//   rec.SetGAliceFile("...");                                               //
//                                                                           //
// The local reconstruction can be switched on or off for individual         //
// detectors by                                                              //
//                                                                           //
//   rec.SetRunLocalReconstruction("...");                                   //
//                                                                           //
// The argument is a (case sensitive) string with the names of the           //
// detectors separated by a space. The special string "ALL" selects all      //
// available detectors. This is the default.                                 //
//                                                                           //
// The reconstruction of the primary vertex position can be switched off by  //
//                                                                           //
//   rec.SetRunVertexFinder(kFALSE);                                         //
//                                                                           //
// The tracking and the creation of ESD tracks can be switched on for        //
// selected detectors by                                                     //
//                                                                           //
//   rec.SetRunTracking("...");                                              //
//                                                                           //
// The filling of additional ESD information can be steered by               //
//                                                                           //
//   rec.SetFillESD("...");                                                  //
//                                                                           //
// Again, for both methods the string specifies the list of detectors.       //
// The default is "ALL".                                                     //
//                                                                           //
// The call of the shortcut method                                           //
//                                                                           //
//   rec.SetRunReconstruction("...");                                        //
//                                                                           //
// is equivalent to calling SetRunLocalReconstruction, SetRunTracking and    //
// SetFillESD with the same detector selecting string as argument.           //
//                                                                           //
// The reconstruction requires digits or raw data as input. For the creation //
// of digits and raw data have a look at the class AliSimulation.            //
//                                                                           //
// For debug purposes the method SetCheckPointLevel can be used. If the      //
// argument is greater than 0, files with ESD events will be written after   //
// selected steps of the reconstruction for each event:                      //
//   level 1: after tracking and after filling of ESD (final)                //
//   level 2: in addition after each tracking step                           //
//   level 3: in addition after the filling of ESD for each detector         //
// If a final check point file exists for an event, this event will be       //
// skipped in the reconstruction. The tracking and the filling of ESD for    //
// a detector will be skipped as well, if the corresponding check point      //
// file exists. The ESD event will then be loaded from the file instead.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TArrayF.h>
#include <TFile.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TPluginManager.h>
#include <TStopwatch.h>

#include "AliReconstruction.h"
#include "AliReconstructor.h"
#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliRawReaderFile.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include "AliTracker.h"
#include "AliESD.h"
#include "AliESDVertex.h"
#include "AliVertexer.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliESDpid.h"
#include "AliMagF.h"

ClassImp(AliReconstruction)


//_____________________________________________________________________________
const char* AliReconstruction::fgkDetectorName[AliReconstruction::fgkNDetectors] = {"ITS", "TPC", "TRD", "TOF", "PHOS", "RICH", "EMCAL", "MUON", "FMD", "ZDC", "PMD", "START", "VZERO", "CRT", "HLT"};

//_____________________________________________________________________________
AliReconstruction::AliReconstruction(const char* gAliceFilename,
				     const char* name, const char* title) :
  TNamed(name, title),

  fRunLocalReconstruction("ALL"),
  fRunVertexFinder(kTRUE),
  fRunTracking("ALL"),
  fFillESD("ALL"),
  fGAliceFileName(gAliceFilename),
  fInput(""),
  fStopOnError(kFALSE),
  fCheckPointLevel(0),
  fOptions(),

  fRunLoader(NULL),
  fRawReader(NULL),

  fVertexer(NULL)
{
// create reconstruction object with default parameters
  
  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    fReconstructor[iDet] = NULL;
    fLoader[iDet] = NULL;
    fTracker[iDet] = NULL;
  }
}

//_____________________________________________________________________________
AliReconstruction::AliReconstruction(const AliReconstruction& rec) :
  TNamed(rec),

  fRunLocalReconstruction(rec.fRunLocalReconstruction),
  fRunVertexFinder(rec.fRunVertexFinder),
  fRunTracking(rec.fRunTracking),
  fFillESD(rec.fFillESD),
  fGAliceFileName(rec.fGAliceFileName),
  fInput(rec.fInput),
  fStopOnError(rec.fStopOnError),
  fCheckPointLevel(0),
  fOptions(),

  fRunLoader(NULL),
  fRawReader(NULL),

  fVertexer(NULL)
{
// copy constructor

  for (Int_t i = 0; i < fOptions.GetEntriesFast(); i++) {
    if (rec.fOptions[i]) fOptions.Add(rec.fOptions[i]->Clone());
  }
  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    fReconstructor[iDet] = NULL;
    fLoader[iDet] = NULL;
    fTracker[iDet] = NULL;
  }
}

//_____________________________________________________________________________
AliReconstruction& AliReconstruction::operator = (const AliReconstruction& rec)
{
// assignment operator

  this->~AliReconstruction();
  new(this) AliReconstruction(rec);
  return *this;
}

//_____________________________________________________________________________
AliReconstruction::~AliReconstruction()
{
// clean up

  CleanUp();
  fOptions.Delete();
}


//_____________________________________________________________________________
void AliReconstruction::SetGAliceFile(const char* fileName)
{
// set the name of the galice file

  fGAliceFileName = fileName;
}

//_____________________________________________________________________________
void AliReconstruction::SetOption(const char* detector, const char* option)
{
// set options for the reconstruction of a detector

  TObject* obj = fOptions.FindObject(detector);
  if (obj) fOptions.Remove(obj);
  fOptions.Add(new TNamed(detector, option));
}


//_____________________________________________________________________________
Bool_t AliReconstruction::Run(const char* input)
{
// run the reconstruction

  // set the input
  if (!input) input = fInput.Data();
  TString fileName(input);
  if (fileName.EndsWith("/")) {
    fRawReader = new AliRawReaderFile(fileName);
  } else if (fileName.EndsWith(".root")) {
    fRawReader = new AliRawReaderRoot(fileName);
  } else if (!fileName.IsNull()) {
    fRawReader = new AliRawReaderDate(fileName);
    fRawReader->SelectEvents(7);
  }

  // get the run loader
  if (!InitRunLoader()) return kFALSE;

  // local reconstruction
  if (!fRunLocalReconstruction.IsNull()) {
    if (!RunLocalReconstruction(fRunLocalReconstruction)) {
      if (fStopOnError) {CleanUp(); return kFALSE;}
    }
  }
  if (!fRunVertexFinder && fRunTracking.IsNull() && 
      fFillESD.IsNull()) return kTRUE;

  // get vertexer
  if (fRunVertexFinder && !CreateVertexer()) {
    if (fStopOnError) {
      CleanUp(); 
      return kFALSE;
    }
  }

  // get trackers
  if (!fRunTracking.IsNull() && !CreateTrackers(fRunTracking)) {
    if (fStopOnError) {
      CleanUp(); 
      return kFALSE;
    }      
  }

  // create the ESD output file and tree
  TFile* file = TFile::Open("AliESDs.root", "RECREATE");
  if (!file->IsOpen()) {
    AliError("opening AliESDs.root failed");
    if (fStopOnError) {CleanUp(file); return kFALSE;}    
  }
  AliESD* esd = new AliESD;
  TTree* tree = new TTree("esdTree", "Tree with ESD objects");
  tree->Branch("ESD", "AliESD", &esd);
  delete esd;
  gROOT->cd();

  // loop over events
  if (fRawReader) fRawReader->RewindEvents();
  
  for (Int_t iEvent = 0; iEvent < fRunLoader->GetNumberOfEvents(); iEvent++) {
    AliInfo(Form("processing event %d", iEvent));
    fRunLoader->GetEvent(iEvent);
    if (fRawReader) fRawReader->NextEvent();

    char fileName[256];
    sprintf(fileName, "ESD_%d.%d_final.root", 
	    fRunLoader->GetHeader()->GetRun(), 
	    fRunLoader->GetHeader()->GetEventNrInRun());
    if (!gSystem->AccessPathName(fileName)) continue;

    esd = new AliESD;
    esd->SetRunNumber(fRunLoader->GetHeader()->GetRun());
    esd->SetEventNumber(fRunLoader->GetHeader()->GetEventNrInRun());
    if (gAlice) {
      esd->SetMagneticField(gAlice->Field()->SolenoidField());
    } else {
      // ???
    }

    // vertex finder
    if (fRunVertexFinder) {
      if (!ReadESD(esd, "vertex")) {
	if (!RunVertexFinder(esd)) {
	  if (fStopOnError) {CleanUp(file); return kFALSE;}
	}
	if (fCheckPointLevel > 0) WriteESD(esd, "vertex");
      }
    }

    // barrel tracking
    if (!fRunTracking.IsNull()) {
      if (!ReadESD(esd, "tracking")) {
	if (!RunTracking(esd)) {
	  if (fStopOnError) {CleanUp(file); return kFALSE;}
	}
	if (fCheckPointLevel > 0) WriteESD(esd, "tracking");
      }
    }

    // fill ESD
    if (!fFillESD.IsNull()) {
      if (!FillESD(esd, fFillESD)) {
	if (fStopOnError) {CleanUp(file); return kFALSE;}
      }
    }

    // combined PID
    AliESDpid::MakePID(esd);
    if (fCheckPointLevel > 1) WriteESD(esd, "PID");

    // write ESD
    tree->Fill();

    if (fCheckPointLevel > 0) WriteESD(esd, "final");
    delete esd;
  }

  file->cd();
  tree->Write();
  CleanUp(file);

  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliReconstruction::RunLocalReconstruction(const TString& detectors)
{
// run the local reconstruction

  TStopwatch stopwatch;
  stopwatch.Start();

  TString detStr = detectors;
  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
    AliReconstructor* reconstructor = GetReconstructor(iDet);
    if (!reconstructor) continue;

    AliInfo(Form("running reconstruction for %s", fgkDetectorName[iDet]));
    TStopwatch stopwatchDet;
    stopwatchDet.Start();
    if (fRawReader) {
      fRawReader->RewindEvents();
      reconstructor->Reconstruct(fRunLoader, fRawReader);
    } else {
      reconstructor->Reconstruct(fRunLoader);
    }
    AliInfo(Form("execution time for %s:", fgkDetectorName[iDet]));
    ToAliInfo(stopwatchDet.Print());
  }

  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    AliError(Form("the following detectors were not found: %s",
                  detStr.Data()));
    if (fStopOnError) return kFALSE;
  }

  AliInfo("execution time:");
  ToAliInfo(stopwatch.Print());

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunVertexFinder(AliESD*& esd)
{
// run the barrel tracking

  TStopwatch stopwatch;
  stopwatch.Start();

  AliESDVertex* vertex = NULL;
  Double_t vtxPos[3] = {0, 0, 0};
  Double_t vtxErr[3] = {0.07, 0.07, 0.1};
  TArrayF mcVertex(3); 
  if (fRunLoader->GetHeader() && fRunLoader->GetHeader()->GenEventHeader()) {
    fRunLoader->GetHeader()->GenEventHeader()->PrimaryVertex(mcVertex);
    for (Int_t i = 0; i < 3; i++) vtxPos[i] = mcVertex[i];
  }

  if (fVertexer) {
    AliInfo("running the ITS vertex finder");
    vertex = fVertexer->FindVertexForCurrentEvent(fRunLoader->GetEventNumber());
    if(!vertex){
      AliWarning("Vertex not found");
      vertex = new AliESDVertex();
    }
    else {
      vertex->SetTruePos(vtxPos);  // store also the vertex from MC
    }

  } else {
    AliInfo("getting the primary vertex from MC");
    vertex = new AliESDVertex(vtxPos, vtxErr);
  }

  if (vertex) {
    vertex->GetXYZ(vtxPos);
    vertex->GetSigmaXYZ(vtxErr);
  } else {
    AliWarning("no vertex reconstructed");
    vertex = new AliESDVertex(vtxPos, vtxErr);
  }
  esd->SetVertex(vertex);
  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    if (fTracker[iDet]) fTracker[iDet]->SetVertex(vtxPos, vtxErr);
  }  
  delete vertex;

  AliInfo("execution time:");
  ToAliInfo(stopwatch.Print());

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunTracking(AliESD*& esd)
{
// run the barrel tracking

  TStopwatch stopwatch;
  stopwatch.Start();

  AliInfo("running tracking");

  // pass 1: TPC + ITS inwards
  for (Int_t iDet = 1; iDet >= 0; iDet--) {
    if (!fTracker[iDet]) continue;
    AliDebug(1, Form("%s tracking", fgkDetectorName[iDet]));

    // load clusters
    fLoader[iDet]->LoadRecPoints("read");
    TTree* tree = fLoader[iDet]->TreeR();
    if (!tree) {
      AliError(Form("Can't get the %s cluster tree", fgkDetectorName[iDet]));
      return kFALSE;
    }
    fTracker[iDet]->LoadClusters(tree);

    // run tracking
    if (fTracker[iDet]->Clusters2Tracks(esd) != 0) {
      AliError(Form("%s Clusters2Tracks failed", fgkDetectorName[iDet]));
      return kFALSE;
    }
    if (fCheckPointLevel > 1) {
      WriteESD(esd, Form("%s.tracking", fgkDetectorName[iDet]));
    }
  }

  // pass 2: ALL backwards
  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    if (!fTracker[iDet]) continue;
    AliDebug(1, Form("%s back propagation", fgkDetectorName[iDet]));

    // load clusters
    if (iDet > 1) {     // all except ITS, TPC
      TTree* tree = NULL;
      if (iDet == 3) {   // TOF
	fLoader[iDet]->LoadDigits("read");
	tree = fLoader[iDet]->TreeD();
      } else {
	fLoader[iDet]->LoadRecPoints("read");
	tree = fLoader[iDet]->TreeR();
      }
      if (!tree) {
	AliError(Form("Can't get the %s cluster tree", fgkDetectorName[iDet]));
	return kFALSE;
      }
      fTracker[iDet]->LoadClusters(tree);
    }

    // run tracking
    if (fTracker[iDet]->PropagateBack(esd) != 0) {
      AliError(Form("%s backward propagation failed", fgkDetectorName[iDet]));
      return kFALSE;
    }
    if (fCheckPointLevel > 1) {
      WriteESD(esd, Form("%s.back", fgkDetectorName[iDet]));
    }

    // unload clusters
    if (iDet > 2) {     // all except ITS, TPC, TRD
      fTracker[iDet]->UnloadClusters();
      if (iDet == 3) {   // TOF
	fLoader[iDet]->UnloadDigits();
      } else {
	fLoader[iDet]->UnloadRecPoints();
      }
    }
  }

  // pass 3: TRD + TPC + ITS refit inwards
  for (Int_t iDet = 2; iDet >= 0; iDet--) {
    if (!fTracker[iDet]) continue;
    AliDebug(1, Form("%s inward refit", fgkDetectorName[iDet]));

    // run tracking
    if (fTracker[iDet]->RefitInward(esd) != 0) {
      AliError(Form("%s inward refit failed", fgkDetectorName[iDet]));
      return kFALSE;
    }
    if (fCheckPointLevel > 1) {
      WriteESD(esd, Form("%s.refit", fgkDetectorName[iDet]));
    }

    // unload clusters
    fTracker[iDet]->UnloadClusters();
    fLoader[iDet]->UnloadRecPoints();
  }

  AliInfo("execution time:");
  ToAliInfo(stopwatch.Print());

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::FillESD(AliESD*& esd, const TString& detectors)
{
// fill the event summary data

  TStopwatch stopwatch;
  stopwatch.Start();
  AliInfo("filling ESD");

  TString detStr = detectors;
  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
    AliReconstructor* reconstructor = GetReconstructor(iDet);
    if (!reconstructor) continue;

    if (!ReadESD(esd, fgkDetectorName[iDet])) {
      AliDebug(1, Form("filling ESD for %s", fgkDetectorName[iDet]));
      if (fRawReader) {
        reconstructor->FillESD(fRunLoader, fRawReader, esd);
      } else {
        reconstructor->FillESD(fRunLoader, esd);
      }
      if (fCheckPointLevel > 2) WriteESD(esd, fgkDetectorName[iDet]);
    }
  }

  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    AliError(Form("the following detectors were not found: %s", 
                  detStr.Data()));
    if (fStopOnError) return kFALSE;
  }

  AliInfo("execution time:");
  ToAliInfo(stopwatch.Print());

  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliReconstruction::IsSelected(TString detName, TString& detectors) const
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

//_____________________________________________________________________________
Bool_t AliReconstruction::InitRunLoader()
{
// get or create the run loader

  if (gAlice) delete gAlice;
  gAlice = NULL;

  if (!gSystem->AccessPathName(fGAliceFileName.Data())) {
    fRunLoader = AliRunLoader::Open(fGAliceFileName.Data());
    if (!fRunLoader) {
      AliError(Form("no run loader found in file %s", fGAliceFileName.Data()));
      CleanUp();
      return kFALSE;
    }
    if (fRunLoader->LoadgAlice() == 0) {
      gAlice = fRunLoader->GetAliRun();
      AliTracker::SetFieldMap(gAlice->Field());
    }
    if (!gAlice && !fRawReader) {
      AliError(Form("no gAlice object found in file %s",
		    fGAliceFileName.Data()));
      CleanUp();
      return kFALSE;
    }

  } else {               // galice.root does not exist
    if (!fRawReader) {
      AliError(Form("the file %s does not exist", fGAliceFileName.Data()));
      CleanUp();
      return kFALSE;
    }
    fRunLoader = AliRunLoader::Open(fGAliceFileName.Data(),
				    AliConfig::GetDefaultEventFolderName(),
				    "recreate");
    if (!fRunLoader) {
      AliError(Form("could not create run loader in file %s", 
		    fGAliceFileName.Data()));
      CleanUp();
      return kFALSE;
    }
    fRunLoader->MakeTree("E");
    Int_t iEvent = 0;
    while (fRawReader->NextEvent()) {
      fRunLoader->SetEventNumber(iEvent);
      fRunLoader->GetHeader()->Reset(fRawReader->GetRunNumber(), 
				     iEvent, iEvent);
      fRunLoader->MakeTree("H");
      fRunLoader->TreeE()->Fill();
      iEvent++;
    }
    fRawReader->RewindEvents();
    fRunLoader->WriteHeader("OVERWRITE");
    fRunLoader->CdGAFile();
    fRunLoader->Write(0, TObject::kOverwrite);
//    AliTracker::SetFieldMap(???);
  }

  return kTRUE;
}

//_____________________________________________________________________________
AliReconstructor* AliReconstruction::GetReconstructor(Int_t iDet)
{
// get the reconstructor object and the loader for a detector

  if (fReconstructor[iDet]) return fReconstructor[iDet];

  // load the reconstructor object
  TPluginManager* pluginManager = gROOT->GetPluginManager();
  TString detName = fgkDetectorName[iDet];
  TString recName = "Ali" + detName + "Reconstructor";
  if (gAlice && !gAlice->GetDetector(detName) && (detName != "HLT")) return NULL;

  if (detName == "HLT") {
    if (!gROOT->GetClass("AliLevel3")) {
      gSystem->Load("libAliL3Src.so");
      gSystem->Load("libAliL3Misc.so");
      gSystem->Load("libAliL3Hough.so");
      gSystem->Load("libAliL3Comp.so");
    }
  }

  AliReconstructor* reconstructor = NULL;
  // first check if a plugin is defined for the reconstructor
  TPluginHandler* pluginHandler = 
    pluginManager->FindHandler("AliReconstructor", detName);
  // if not, add a plugin for it
  if (!pluginHandler) {
    AliDebug(1, Form("defining plugin for %s", recName.Data()));
    if (gSystem->Load("lib" + detName + "base.so") == 0) {
      pluginManager->AddHandler("AliReconstructor", detName, 
				recName, detName + "rec", recName + "()");
    } else {
      pluginManager->AddHandler("AliReconstructor", detName, 
				recName, detName, recName + "()");
    }
    pluginHandler = pluginManager->FindHandler("AliReconstructor", detName);
  }
  if (pluginHandler && (pluginHandler->LoadPlugin() == 0)) {
    reconstructor = (AliReconstructor*) pluginHandler->ExecPlugin(0);
  }
  if (reconstructor) {
    TObject* obj = fOptions.FindObject(detName.Data());
    if (obj) reconstructor->SetOption(obj->GetTitle());
    fReconstructor[iDet] = reconstructor;
  }

  // get or create the loader
  if (detName != "HLT") {
    fLoader[iDet] = fRunLoader->GetLoader(detName + "Loader");
    if (!fLoader[iDet]) {
      AliConfig::Instance()
	->CreateDetectorFolders(fRunLoader->GetEventFolder(), 
				detName, detName);
      // first check if a plugin is defined for the loader
      TPluginHandler* pluginHandler = 
	pluginManager->FindHandler("AliLoader", detName);
      // if not, add a plugin for it
      if (!pluginHandler) {
	TString loaderName = "Ali" + detName + "Loader";
	AliDebug(1, Form("defining plugin for %s", loaderName.Data()));
	pluginManager->AddHandler("AliLoader", detName, 
				  loaderName, detName + "base", 
				  loaderName + "(const char*, TFolder*)");
	pluginHandler = pluginManager->FindHandler("AliLoader", detName);
      }
      if (pluginHandler && (pluginHandler->LoadPlugin() == 0)) {
	fLoader[iDet] = 
	  (AliLoader*) pluginHandler->ExecPlugin(2, detName.Data(), 
						 fRunLoader->GetEventFolder());
      }
      if (!fLoader[iDet]) {   // use default loader
	fLoader[iDet] = new AliLoader(detName, fRunLoader->GetEventFolder());
      }
      if (!fLoader[iDet]) {
	AliWarning(Form("couldn't get loader for %s", detName.Data()));
	if (fStopOnError) return NULL;
      } else {
	fRunLoader->AddLoader(fLoader[iDet]);
	fRunLoader->CdGAFile();
	if (gFile && !gFile->IsWritable()) gFile->ReOpen("UPDATE");
	fRunLoader->Write(0, TObject::kOverwrite);
      }
    }
  }
      
  return reconstructor;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::CreateVertexer()
{
// create the vertexer

  fVertexer = NULL;
  AliReconstructor* itsReconstructor = GetReconstructor(0);
  if (itsReconstructor) {
    fVertexer = itsReconstructor->CreateVertexer(fRunLoader);
  }
  if (!fVertexer) {
    AliWarning("couldn't create a vertexer for ITS");
    if (fStopOnError) return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::CreateTrackers(const TString& detectors)
{
// create the trackers

  TString detStr = detectors;
  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    if (!IsSelected(fgkDetectorName[iDet], detStr)) continue;
    AliReconstructor* reconstructor = GetReconstructor(iDet);
    if (!reconstructor) continue;
    TString detName = fgkDetectorName[iDet];
    if (detName == "HLT") continue;

    fTracker[iDet] = reconstructor->CreateTracker(fRunLoader);
    if (!fTracker[iDet] && (iDet < 7)) {
      AliWarning(Form("couldn't create a tracker for %s", detName.Data()));
      if (fStopOnError) return kFALSE;
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliReconstruction::CleanUp(TFile* file)
{
// delete trackers and the run loader and close and delete the file

  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    delete fReconstructor[iDet];
    fReconstructor[iDet] = NULL;
    fLoader[iDet] = NULL;
    delete fTracker[iDet];
    fTracker[iDet] = NULL;
  }
  delete fVertexer;
  fVertexer = NULL;

  delete fRunLoader;
  fRunLoader = NULL;
  delete fRawReader;
  fRawReader = NULL;

  if (file) {
    file->Close();
    delete file;
  }
}


//_____________________________________________________________________________
Bool_t AliReconstruction::ReadESD(AliESD*& esd, const char* recStep) const
{
// read the ESD event from a file

  if (!esd) return kFALSE;
  char fileName[256];
  sprintf(fileName, "ESD_%d.%d_%s.root", 
	  esd->GetRunNumber(), esd->GetEventNumber(), recStep);
  if (gSystem->AccessPathName(fileName)) return kFALSE;

  AliDebug(1, Form("reading ESD from file %s", fileName));
  TFile* file = TFile::Open(fileName);
  if (!file || !file->IsOpen()) {
    AliError(Form("opening %s failed", fileName));
    delete file;
    return kFALSE;
  }

  gROOT->cd();
  delete esd;
  esd = (AliESD*) file->Get("ESD");
  file->Close();
  delete file;
  return kTRUE;
}

//_____________________________________________________________________________
void AliReconstruction::WriteESD(AliESD* esd, const char* recStep) const
{
// write the ESD event to a file

  if (!esd) return;
  char fileName[256];
  sprintf(fileName, "ESD_%d.%d_%s.root", 
	  esd->GetRunNumber(), esd->GetEventNumber(), recStep);

  AliDebug(1, Form("writing ESD to file %s", fileName));
  TFile* file = TFile::Open(fileName, "recreate");
  if (!file || !file->IsOpen()) {
    AliError(Form("opening %s failed", fileName));
  } else {
    esd->Write("ESD");
    file->Close();
  }
  delete file;
}
