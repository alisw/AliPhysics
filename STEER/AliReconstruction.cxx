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
// The name of the galice file can be changed from the default               //
// "galice.root" by                                                          //
//                                                                           //
//   rec.SetGAliceFile("...");                                               //
//                                                                           //
// The reconstruction can be switched on or off for individual detectors by  //
//                                                                           //
//   rec.SetRunReconstruction("...");                                        //
//                                                                           //
// The argument is a (case sensitive) string with the names of the           //
// detectors separated by a space. The special string "ALL" selects all      //
// available detectors. This is the default.                                 //
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
// The reconstruction requires digits as input. For the creation of digits   //
// have a look at the class AliSimulation.                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliReconstruction.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliModule.h"
#include "AliDetector.h"
#include "AliTracker.h"
#include "AliESD.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliESDpid.h"
#include <TArrayF.h>


ClassImp(AliReconstruction)


//_____________________________________________________________________________
AliReconstruction::AliReconstruction(const char* name, const char* title) :
  TNamed(name, title)
{
// create reconstruction object with default parameters

  Init();
}

//_____________________________________________________________________________
AliReconstruction::AliReconstruction(const AliReconstruction& rec) :
  TNamed(rec)
{
// copy constructor

  fRunReconstruction = rec.fRunReconstruction;
  fRunTracking = rec.fRunTracking;
  fStopOnError = rec.fStopOnError;

  fGAliceFileName = rec.fGAliceFileName;

  fRunLoader = NULL;
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

}

//_____________________________________________________________________________
void AliReconstruction::Init()
{
// set default parameters

  fRunReconstruction = "ALL";
  fRunTracking = kTRUE;
  fFillESD = "ALL";
  fStopOnError = kFALSE;

  fGAliceFileName = "galice.root";

  fRunLoader = NULL;
}


//_____________________________________________________________________________
void AliReconstruction::SetGAliceFile(const char* fileName)
{
// set the name of the galice file

  fGAliceFileName = fileName;
}


//_____________________________________________________________________________
Bool_t AliReconstruction::Run()
{
// run the reconstruction

  // open the run loader
  if (fRunLoader) delete fRunLoader;
  fRunLoader = AliRunLoader::Open(fGAliceFileName.Data());
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

  // local reconstruction
  if (!fRunReconstruction.IsNull()) {
    if (!RunReconstruction(fRunReconstruction)) {
      if (fStopOnError) return kFALSE;
    }
  }
  if (!fRunTracking && fFillESD.IsNull()) return kTRUE;

  // get loaders and trackers
  fITSLoader = fRunLoader->GetLoader("ITSLoader");
  if (!fITSLoader) {
    Error("Run", "no ITS loader found");
    if (fStopOnError) return kFALSE;
  }
  fITSTracker = NULL;
  if (gAlice->GetDetector("ITS")) {
    fITSTracker = gAlice->GetDetector("ITS")->CreateTracker();
  }
  if (!fITSTracker) {
    Error("Run", "couldn't create a tracker for ITS");
    if (fStopOnError) return kFALSE;
  }

  fTPCLoader = fRunLoader->GetLoader("TPCLoader");
  if (!fTPCLoader) {
    Error("Run", "no TPC loader found");
    if (fStopOnError) return kFALSE;
  }
  fTPCTracker = NULL;
  if (gAlice->GetDetector("TPC")) {
    fTPCTracker = gAlice->GetDetector("TPC")->CreateTracker();
  }
  if (!fTPCTracker) {
    Error("Run", "couldn't create a tracker for TPC");
    if (fStopOnError) return kFALSE;
  }

  fTRDLoader = fRunLoader->GetLoader("TRDLoader");
  if (!fTRDLoader) {
    Error("Run", "no TRD loader found");
    if (fStopOnError) return kFALSE;
  }
  fTRDTracker = NULL;
  if (gAlice->GetDetector("TRD")) {
    fTRDTracker = gAlice->GetDetector("TRD")->CreateTracker();
  }
  if (!fTRDTracker) {
    Error("Run", "couldn't create a tracker for TRD");
    if (fStopOnError) return kFALSE;
  }

  fTOFLoader = fRunLoader->GetLoader("TOFLoader");
  if (!fTOFLoader) {
    Error("Run", "no TOF loader found");
    if (fStopOnError) return kFALSE;
  }
  fTOFTracker = NULL;
  if (gAlice->GetDetector("TOF")) {
    fTOFTracker = gAlice->GetDetector("TOF")->CreateTracker();
  }
  if (!fTOFTracker) {
    Error("Run", "couldn't create a tracker for TOF");
    if (fStopOnError) return kFALSE;
  }

  // create the ESD output file
  TFile* file = TFile::Open("AliESDs.root", "RECREATE");
  if (!file->IsOpen()) {
    Error("Run", "opening AliESDs.root failed");
    if (fStopOnError) return kFALSE;    
  }

  // loop over events
  for (Int_t iEvent = 0; iEvent < fRunLoader->GetNumberOfEvents(); iEvent++) {
    Info("Run", "processing event %d", iEvent);
    AliESD* esd = new AliESD;
    fRunLoader->GetEvent(iEvent);
    esd->SetRunNumber(gAlice->GetRunNumber());
    esd->SetEventNumber(gAlice->GetEvNumber());

    // barrel tracking
    if (fRunTracking) {
      if (!RunTracking(esd)) {
	if (fStopOnError) return kFALSE;
      }
    }

    // fill ESD
    if (!fFillESD.IsNull()) {
      if (!FillESD(esd, fFillESD)) {
	if (fStopOnError) return kFALSE;
      }
    }

    // combined PID
    AliESDpid::MakePID(esd);

    // write ESD
    char name[100]; 
    sprintf(name, "ESD%d", iEvent);
    file->cd();
    if (!esd->Write(name)) {
      Error("Run", "writing ESD failed");
      if (fStopOnError) return kFALSE;
    }
  }

  file->Close();

  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliReconstruction::RunReconstruction(const TString& detectors)
{
// run the reconstruction

  TStopwatch stopwatch;
  stopwatch.Start();

  TString detStr = detectors;
  TObjArray* detArray = gAlice->Detectors();
  for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
    AliModule* det = (AliModule*) detArray->At(iDet);
    if (!det || !det->IsActive()) continue;
    if (IsSelected(det->GetName(), detStr)) {
      Info("RunReconstruction", "running reconstruction for %s", 
	   det->GetName());
      TStopwatch stopwatchDet;
      stopwatchDet.Start();
      det->Reconstruct();
      Info("RunReconstruction", "execution time for %s:", det->GetName());
      stopwatchDet.Print();
    }
  }

  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    Error("RunReconstruction", "the following detectors were not found: %s", 
	  detStr.Data());
    if (fStopOnError) return kFALSE;
  }

  Info("RunReconstruction", "execution time:");
  stopwatch.Print();

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::RunTracking(AliESD* esd)
{
// run the barrel tracking

  TStopwatch stopwatch;
  stopwatch.Start();

  // get the primary vertex (from MC for the moment)
  TArrayF vertex(3);     
  fRunLoader->GetHeader()->GenEventHeader()->PrimaryVertex(vertex);
  Double_t vtxPos[3] = {vertex[0], vertex[1], vertex[2]};
  Double_t vtxCov[6] = {
    0.005,
    0.000, 0.005,
    0.000, 0.000, 0.010
  };
  Double_t vtxErr[3] = {vtxCov[0], vtxCov[2], vtxCov[5]}; // diag. elements
  esd->SetVertex(vtxPos, vtxCov);
  fITSTracker->SetVertex(vtxPos, vtxErr);
  fTPCTracker->SetVertex(vtxPos, vtxErr);
  fTRDTracker->SetVertex(vtxPos, vtxErr);

  // TPC tracking
  Info("RunTracking", "TPC tracking");
  fTPCLoader->LoadRecPoints("read");
  TTree* tpcTree = fTPCLoader->TreeR();
  if (!tpcTree) {
    Error("RunTracking", "Can't get the TPC cluster tree");
    return kFALSE;
  }     
  fTPCTracker->LoadClusters(tpcTree);
  if (fTPCTracker->Clusters2Tracks(esd) != 0) {
    Error("RunTracking", "TPC Clusters2Tracks failed");
    return kFALSE;
  }

  gAlice->GetDetector("TPC")->FillESD(esd); // preliminary PID
  AliESDpid::MakePID(esd);                  // for the ITS tracker

  // ITS tracking
  Info("RunTracking", "ITS tracking");
  fITSLoader->LoadRecPoints("read");
  TTree* itsTree = fITSLoader->TreeR();
  if (!itsTree) {
    Error("RunTracking", "Can't get the ITS cluster tree");
    return kFALSE;
  }     
  fITSTracker->LoadClusters(itsTree);
  if (fITSTracker->Clusters2Tracks(esd) != 0) {
    Error("RunTracking", "ITS Clusters2Tracks failed");
    return kFALSE;
  }

  // ITS back propagation
  Info("RunTracking", "ITS back propagation");
  if (fITSTracker->PropagateBack(esd) != 0) {
    Error("RunTracking", "ITS backward propagation failed");
    return kFALSE;
  }

  // TPC back propagation
  Info("RunTracking", "TPC back propagation");
  if (fTPCTracker->PropagateBack(esd) != 0) {
    Error("RunTracking", "TPC backward propagation failed");
    return kFALSE;
  }

  // TRD back propagation
  Info("RunTracking", "TRD back propagation");
  fTRDLoader->LoadRecPoints("read");
  TTree* trdTree = fTRDLoader->TreeR();
  if (!trdTree) {
    Error("RunTracking", "Can't get the TRD cluster tree");
    return kFALSE;
  }     
  fTRDTracker->LoadClusters(trdTree);
  if (fTRDTracker->PropagateBack(esd) != 0) {
    Error("RunTracking", "TRD backward propagation failed");
    return kFALSE;
  }

  // TOF back propagation
  Info("RunTracking", "TOF back propagation");
  fTOFLoader->LoadDigits("read");
  TTree* tofTree = fTOFLoader->TreeD();
  if (!tofTree) {
    Error("RunTracking", "Can't get the TOF digits tree");
    return kFALSE;
  }     
  fTOFTracker->LoadClusters(tofTree);
  if (fTOFTracker->PropagateBack(esd) != 0) {
    Error("RunTracking", "TOF backward propagation failed");
    return kFALSE;
  }
  fTOFTracker->UnloadClusters();
  fTOFLoader->UnloadDigits();

  // TRD inward refit
  Info("RunTracking", "TRD inward refit");
  if (fTRDTracker->RefitInward(esd) != 0) {
    Error("RunTracking", "TRD inward refit failed");
    return kFALSE;
  }
  fTRDTracker->UnloadClusters();
  fTRDLoader->UnloadRecPoints();
    
  // TPC inward refit
  Info("RunTracking", "TPC inward refit");
  if (fTPCTracker->RefitInward(esd) != 0) {
    Error("RunTracking", "TPC inward refit failed");
    return kFALSE;
  }
  fTPCTracker->UnloadClusters();
  fTPCLoader->UnloadRecPoints();
    
  // ITS inward refit
  Info("RunTracking", "ITS inward refit");
  if (fITSTracker->RefitInward(esd) != 0) {
    Error("RunTracking", "ITS inward refit failed");
    return kFALSE;
  }
  fITSTracker->UnloadClusters();
  fITSLoader->UnloadRecPoints();

  Info("RunTracking", "execution time:");
  stopwatch.Print();

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliReconstruction::FillESD(AliESD* esd, const TString& detectors)
{
// fill the event summary data

  TStopwatch stopwatch;
  stopwatch.Start();

  TString detStr = detectors;
  TObjArray* detArray = gAlice->Detectors();
  for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
    AliModule* det = (AliModule*) detArray->At(iDet);
    if (!det || !det->IsActive()) continue;
    if (IsSelected(det->GetName(), detStr)) {
      Info("FillESD", "filling ESD for %s", 
	   det->GetName());
      det->FillESD(esd);
    }
  }

  if ((detStr.CompareTo("ALL") != 0) && !detStr.IsNull()) {
    Error("FillESD", "the following detectors were not found: %s", 
	  detStr.Data());
    if (fStopOnError) return kFALSE;
  }

  Info("FillESD", "execution time:");
  stopwatch.Print();

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
