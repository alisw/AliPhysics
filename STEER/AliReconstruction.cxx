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
// The tracking in ITS, TPC and TRD and the creation of ESD tracks can be    //
// switched off by                                                           //
//                                                                           //
//   rec.SetRunTracking(kFALSE);                                             //
//                                                                           //
// The filling of additional ESD information can be steered by               //
//                                                                           //
//   rec.SetFillESD("...");                                                  //
//                                                                           //
// Again, the string specifies the list of detectors. The default is "ALL".  //
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
  fRunTracking(kTRUE),
  fFillESD("ALL"),
  fGAliceFileName(gAliceFilename),
  fInput(""),
  fStopOnError(kFALSE),
  fCheckPointLevel(0),

  fRunLoader(NULL),
  fRawReader(NULL),
  fITSLoader(NULL),
  fITSVertexer(NULL),
  fITSTracker(NULL),
  fTPCLoader(NULL),
  fTPCTracker(NULL),
  fTRDLoader(NULL),
  fTRDTracker(NULL),
  fTOFLoader(NULL),
  fTOFTracker(NULL),

  fReconstructors(),
  fOptions()
{
// create reconstruction object with default parameters

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

  fRunLoader(NULL),
  fRawReader(NULL),
  fITSLoader(NULL),
  fITSVertexer(NULL),
  fITSTracker(NULL),
  fTPCLoader(NULL),
  fTPCTracker(NULL),
  fTRDLoader(NULL),
  fTRDTracker(NULL),
  fTOFLoader(NULL),
  fTOFTracker(NULL),

  fReconstructors(),
  fOptions()
{
// copy constructor

  for (Int_t i = 0; i < fOptions.GetEntriesFast(); i++) {
    if (rec.fOptions[i]) fOptions.Add(rec.fOptions[i]->Clone());
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

  // open the run loader
  fRunLoader = AliRunLoader::Open(fGAliceFileName.Data());
  if (!fRunLoader) {
    AliError(Form("no run loader found in file %s", fGAliceFileName.Data()));
    CleanUp();
    return kFALSE;
  }
  fRunLoader->LoadgAlice();
  AliRun* aliRun = fRunLoader->GetAliRun();
  if (!aliRun) {
    AliError(Form("no gAlice object found in file %s",
                  fGAliceFileName.Data()));
    CleanUp();
    return kFALSE;
  }
  gAlice = aliRun;
  AliTracker::SetFieldMap(gAlice->Field());

  // load the reconstructor objects
  TPluginManager* pluginManager = gROOT->GetPluginManager();
  for (Int_t iDet = 0; iDet < fgkNDetectors; iDet++) {
    TString detName = fgkDetectorName[iDet];
    TString recName = "Ali" + detName + "Reconstructor";
    if (!gAlice->GetDetector(detName) && detName != "HLT") continue;

    if(detName == "HLT") {
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
    // if not, but the reconstructor class is implemented, add a plugin for it
    if (!pluginHandler && gROOT->GetClass(recName.Data())) {
      AliDebug(1, Form("defining plugin for %s", recName.Data()));
      pluginManager->AddHandler("AliReconstructor", detName, 
				recName, detName, recName + "()");
      pluginHandler = pluginManager->FindHandler("AliReconstructor", detName);
    }
    if (pluginHandler && (pluginHandler->LoadPlugin() == 0)) {
      reconstructor = (AliReconstructor*) pluginHandler->ExecPlugin(0);
    }
    // if there is no reconstructor class for the detector use the dummy one
    if (!reconstructor && gAlice->GetDetector(detName)) {
      AliDebug(1, Form("using dummy reconstructor for %s", detName.Data()));
      reconstructor = new AliDummyReconstructor(gAlice->GetDetector(detName));
    }
    if (reconstructor) {
      TObject* obj = fOptions.FindObject(detName.Data());
      if (obj) reconstructor->SetOption(obj->GetTitle());
      fReconstructors.Add(reconstructor);
    }
  }

  // local reconstruction
  if (!fRunLocalReconstruction.IsNull()) {
    if (!RunLocalReconstruction(fRunLocalReconstruction)) {
      if (fStopOnError) {CleanUp(); return kFALSE;}
    }
  }
  if (!fRunVertexFinder && !fRunTracking && fFillESD.IsNull()) return kTRUE;

  // get vertexer
  if (fRunVertexFinder && !CreateVertexer()) {
    if (fStopOnError) {
      CleanUp(); 
      return kFALSE;
    }
  }

  // get loaders and trackers
  if (fRunTracking && !CreateTrackers()) {
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
	    aliRun->GetRunNumber(), aliRun->GetEvNumber());
    if (!gSystem->AccessPathName(fileName)) continue;

    esd = new AliESD;
    esd->SetRunNumber(aliRun->GetRunNumber());
    esd->SetEventNumber(aliRun->GetEvNumber());
    esd->SetMagneticField(aliRun->Field()->SolenoidField());

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
    if (fRunTracking) {
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
  for (Int_t iDet = 0; iDet < fReconstructors.GetEntriesFast(); iDet++) {
    AliReconstructor* reconstructor = 
      (AliReconstructor*) fReconstructors[iDet];
    TString detName = reconstructor->GetDetectorName();
    if (IsSelected(detName, detStr)) {
      AliInfo(Form("running reconstruction for %s", detName.Data()));
      TStopwatch stopwatchDet;
      stopwatchDet.Start();
      if (fRawReader) {
	fRawReader->RewindEvents();
	reconstructor->Reconstruct(fRunLoader, fRawReader);
      } else {
	reconstructor->Reconstruct(fRunLoader);
      }
      AliInfo(Form("execution time for %s:", detName.Data()));
      ToAliInfo(stopwatchDet.Print());
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

  if (fITSVertexer) {
    AliInfo("running the ITS vertex finder");
    vertex = fITSVertexer->FindVertexForCurrentEvent(fRunLoader->GetEventNumber());
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
  if (fITSTracker) fITSTracker->SetVertex(vtxPos, vtxErr);
  if (fTPCTracker) fTPCTracker->SetVertex(vtxPos, vtxErr);
  if (fTRDTracker) fTRDTracker->SetVertex(vtxPos, vtxErr);
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

  if (!fTPCTracker) {
    AliError("no TPC tracker");
    return kFALSE;
  }
  AliInfo("running tracking");

  // TPC tracking
  AliDebug(1, "TPC tracking");
  fTPCLoader->LoadRecPoints("read");
  TTree* tpcTree = fTPCLoader->TreeR();
  if (!tpcTree) {
    AliError("Can't get the TPC cluster tree");
    return kFALSE;
  }
  fTPCTracker->LoadClusters(tpcTree);
  if (fTPCTracker->Clusters2Tracks(esd) != 0) {
    AliError("TPC Clusters2Tracks failed");
    return kFALSE;
  }
  if (fCheckPointLevel > 1) WriteESD(esd, "TPC.tracking");

  if (!fITSTracker) {
    AliWarning("no ITS tracker");
  } else {

    GetReconstructor("TPC")->FillESD(fRunLoader, esd); // preliminary
    AliESDpid::MakePID(esd);                  // PID for the ITS tracker

    // ITS tracking
    AliDebug(1, "ITS tracking");
    fITSLoader->LoadRecPoints("read");
    TTree* itsTree = fITSLoader->TreeR();
    if (!itsTree) {
      Error("RunTracking", "Can't get the ITS cluster tree");
      return kFALSE;
    }
    fITSTracker->LoadClusters(itsTree);
    if (fITSTracker->Clusters2Tracks(esd) != 0) {
      AliError("ITS Clusters2Tracks failed");
      return kFALSE;
    }
    if (fCheckPointLevel > 1) WriteESD(esd, "ITS.tracking");

    if (!fTRDTracker) {
      AliWarning("no TRD tracker");
    } else {
      // ITS back propagation
      AliDebug(1, "ITS back propagation");
      if (fITSTracker->PropagateBack(esd) != 0) {
	AliError("ITS backward propagation failed");
	return kFALSE;
      }
      if (fCheckPointLevel > 1) WriteESD(esd, "ITS.back");

      // TPC back propagation
      AliDebug(1, "TPC back propagation");
      if (fTPCTracker->PropagateBack(esd) != 0) {
	AliError("TPC backward propagation failed");
	return kFALSE;
      }
      if (fCheckPointLevel > 1) WriteESD(esd, "TPC.back");

      // TRD back propagation
      AliDebug(1, "TRD back propagation");
      fTRDLoader->LoadRecPoints("read");
      TTree* trdTree = fTRDLoader->TreeR();
      if (!trdTree) {
	AliError("Can't get the TRD cluster tree");
	return kFALSE;
      }
      fTRDTracker->LoadClusters(trdTree);
      if (fTRDTracker->PropagateBack(esd) != 0) {
	AliError("TRD backward propagation failed");
	return kFALSE;
      }
      if (fCheckPointLevel > 1) WriteESD(esd, "TRD.back");

      if (!fTOFTracker) {
	AliWarning("no TOF tracker");
      } else {
	// TOF back propagation
	AliDebug(1, "TOF back propagation");
	fTOFLoader->LoadDigits("read");
	TTree* tofTree = fTOFLoader->TreeD();
	if (!tofTree) {
	  AliError("Can't get the TOF digits tree");
	  return kFALSE;
	}
	fTOFTracker->LoadClusters(tofTree);
	if (fTOFTracker->PropagateBack(esd) != 0) {
	  AliError("TOF backward propagation failed");
	  return kFALSE;
	}
	if (fCheckPointLevel > 1) WriteESD(esd, "TOF.back");
	fTOFTracker->UnloadClusters();
	fTOFLoader->UnloadDigits();
      }

      // TRD inward refit
      AliDebug(1, "TRD inward refit");
      if (fTRDTracker->RefitInward(esd) != 0) {
	AliError("TRD inward refit failed");
	return kFALSE;
      }
      if (fCheckPointLevel > 1) WriteESD(esd, "TRD.refit");
      fTRDTracker->UnloadClusters();
      fTRDLoader->UnloadRecPoints();
    
      // TPC inward refit
      AliInfo("TPC inward refit");
      if (fTPCTracker->RefitInward(esd) != 0) {
	AliError("TPC inward refit failed");
	return kFALSE;
      }
      if (fCheckPointLevel > 1) WriteESD(esd, "TPC.refit");
    
      // ITS inward refit
      AliInfo("ITS inward refit");
      if (fITSTracker->RefitInward(esd) != 0) {
	AliError("ITS inward refit failed");
	return kFALSE;
      }
      if (fCheckPointLevel > 1) WriteESD(esd, "ITS.refit");

    }  // if TRD tracker
    fITSTracker->UnloadClusters();
    fITSLoader->UnloadRecPoints();

  }  // if ITS tracker
  fTPCTracker->UnloadClusters();
  fTPCLoader->UnloadRecPoints();

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
  for (Int_t iDet = 0; iDet < fReconstructors.GetEntriesFast(); iDet++) {
    AliReconstructor* reconstructor = 
      (AliReconstructor*) fReconstructors[iDet];
    TString detName = reconstructor->GetDetectorName();
    if (IsSelected(detName, detStr)) {
      if (!ReadESD(esd, detName.Data())) {
	AliDebug(1, Form("filling ESD for %s", detName.Data()));
	if (fRawReader) {
	  reconstructor->FillESD(fRunLoader, fRawReader, esd);
	} else {
	  reconstructor->FillESD(fRunLoader, esd);
	}
	if (fCheckPointLevel > 2) WriteESD(esd, detName.Data());
      }
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
AliReconstructor* AliReconstruction::GetReconstructor(const char* detName) const
{
// get the reconstructor object for a detector

  for (Int_t iDet = 0; iDet < fReconstructors.GetEntriesFast(); iDet++) {
    AliReconstructor* reconstructor = 
      (AliReconstructor*) fReconstructors[iDet];
    if (strcmp(reconstructor->GetDetectorName(), detName) == 0) {
      return reconstructor;
    }
  }
  return NULL;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::CreateVertexer()
{
// create the vertexer

  fITSVertexer = NULL;
  AliReconstructor* itsReconstructor = GetReconstructor("ITS");
  if (itsReconstructor) {
    fITSVertexer = itsReconstructor->CreateVertexer(fRunLoader);
  }
  if (!fITSVertexer) {
    AliWarning("couldn't create a vertexer for ITS");
    if (fStopOnError) return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::CreateTrackers()
{
// get the loaders and create the trackers

  fITSTracker = NULL;
  fITSLoader = fRunLoader->GetLoader("ITSLoader");
  if (!fITSLoader) {
    AliWarning("no ITS loader found");
    if (fStopOnError) return kFALSE;
  } else {
    AliReconstructor* itsReconstructor = GetReconstructor("ITS");
    if (itsReconstructor) {
      fITSTracker = itsReconstructor->CreateTracker(fRunLoader);
    }
    if (!fITSTracker) {
      AliWarning("couldn't create a tracker for ITS");
      if (fStopOnError) return kFALSE;
    }
  }
    
  fTPCTracker = NULL;
  fTPCLoader = fRunLoader->GetLoader("TPCLoader");
  if (!fTPCLoader) {
    AliError("no TPC loader found");
    if (fStopOnError) return kFALSE;
  } else {
    AliReconstructor* tpcReconstructor = GetReconstructor("TPC");
    if (tpcReconstructor) {
      fTPCTracker = tpcReconstructor->CreateTracker(fRunLoader);
    }
    if (!fTPCTracker) {
      AliError("couldn't create a tracker for TPC");
      if (fStopOnError) return kFALSE;
    }
  }
    
  fTRDTracker = NULL;
  fTRDLoader = fRunLoader->GetLoader("TRDLoader");
  if (!fTRDLoader) {
    AliWarning("no TRD loader found");
    if (fStopOnError) return kFALSE;
  } else {
    AliReconstructor* trdReconstructor = GetReconstructor("TRD");
    if (trdReconstructor) {
      fTRDTracker = trdReconstructor->CreateTracker(fRunLoader);
    }
    if (!fTRDTracker) {
      AliWarning("couldn't create a tracker for TRD");
      if (fStopOnError) return kFALSE;
    }
  }
    
  fTOFTracker = NULL;
  fTOFLoader = fRunLoader->GetLoader("TOFLoader");
  if (!fTOFLoader) {
    AliWarning("no TOF loader found");
    if (fStopOnError) return kFALSE;
  } else {
    AliReconstructor* tofReconstructor = GetReconstructor("TOF");
    if (tofReconstructor) {
      fTOFTracker = tofReconstructor->CreateTracker(fRunLoader);
    }
    if (!fTOFTracker) {
      AliWarning("couldn't create a tracker for TOF");
      if (fStopOnError) return kFALSE;
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliReconstruction::CleanUp(TFile* file)
{
// delete trackers and the run loader and close and delete the file

  fReconstructors.Delete();

  delete fITSVertexer;
  fITSVertexer = NULL;
  delete fITSTracker;
  fITSTracker = NULL;
  delete fTPCTracker;
  fTPCTracker = NULL;
  delete fTRDTracker;
  fTRDTracker = NULL;
  delete fTOFTracker;
  fTOFTracker = NULL;

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
