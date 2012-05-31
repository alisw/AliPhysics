////////////////////////////////////////////////////////////////////////
//
// AliTPCTracking.C 
//
// date: 22.08.2002
// author: Jiri Chudoba based on code of Jourij Belikov
// version: 1.0
// description: 
//      reconstructs of tracks in TPC in the following steps:
//         TPC cluster finding
//         TPC track finding
// input parameters: 
//        Int_t nEvents      ... nr of events to process (<0 means all)
//        Int_t firstEventNr ... first event number (starts from 0)
//        const char* fileName ... name of galice file
//        Bool_t makeClusters ... run the cluster finder or not?
//        Bool_t makeTracks ... run the track finder or not?
//        const char* fileNameRaw ... if not NULL, the cluster finder uses
//                                    the given file as raw data input
//
// History:
//
//     21.07.2003 ... NewIO
//
//     18.03.2003 ... Char_t* replaced by const char*
//
//     03.03.2003 ... SetFieldFactor moved to AliTracker class and
//                    LoadTPCParam moved to AliTPC class
//                    TString replaced by Char_t*
//
//     20.11.2002 ... Changes due to a changed interface of AliTPCtracker. 
//                    Use Riostream.h instead of iostream.h
//
//     22.08.2002 ... first version
//
////////////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TBenchmark.h>
#include "AliRunLoader.h"
#include "AliTPC.h"
#include "AliTPCParam.h"
#include "AliTPCclustererMI.h"
#include "AliTPCtrackerMI.h"
#include "AliRawReaderRoot.h"
#endif

Bool_t AliTPCTracking(Int_t nEvents = -1, Int_t firstEvent = 0,
		      const char* fileName = "galice.root",
		      Bool_t makeClusters = kTRUE,
		      Bool_t makeTracks = kTRUE,
		      const char* fileNameRaw = NULL);

/*
Int_t TPCRefitInward(Int_t nEvents=1, Int_t firstEvent=0,
		     const char* fileNameClusters="tpc.clusters.root",
		     const char* fileNameTracks="tpc.tracks.root",
		     const char* fileNameTRDTracks="trd.tracks.root",
		     const char* fileNameRefittedTracks="tpc.refitted.tracks.root");
*/

////////////////////////////////////////////////////////////////////////
Bool_t AliTPCTracking(Int_t nEvents, Int_t firstEvent,
		      const char* fileName,
		      Bool_t makeClusters,
		      Bool_t makeTracks,
		      const char* fileNameRaw) 
{
  // get the loaders
  AliRunLoader* runLoader = AliRunLoader::Open(fileName);
  if (!runLoader) {
    cerr << "AliTPCTracking: no run loader found\n";
    return kFALSE;
  }
  AliLoader* tpcLoader = runLoader->GetLoader("TPCLoader");
  if (!tpcLoader) {
    cerr << "AliTPCTracking: no TPC loader found\n";
    return kFALSE;
  }

  // get the TPC parameters
  runLoader->CdGAFile();
  AliTPCParam* param = AliTPC::LoadTPCParam(gFile);
  if (!param) {
    cerr << "AliTPCTracking: no TPC parameters found\n";
    return kFALSE;
  }

  // create the clusterer object
  AliTPCclustererMI* clusterer = NULL;
  if (makeClusters) {
    clusterer = new AliTPCclustererMI(param);
    if (!fileNameRaw) tpcLoader->LoadDigits();
    tpcLoader->LoadRecPoints("recreate");
  }

  // create the tracker object
  AliTPCtrackerMI* tracker = NULL;
  if (makeTracks) {
//    tracker = new AliTPCtrackerMI(param);
    if (!makeClusters) tpcLoader->LoadRecPoints();
    tpcLoader->LoadTracks("recreate");
  }

  // get the event number range
  Int_t maxEvent = 0;
  if (fileNameRaw) {
    TFile* file = TFile::Open(fileNameRaw);
    if (file && file->IsOpen()) {
      TTree* tree = (TTree*) file->Get("T");
      if (tree) maxEvent = (Int_t) tree->GetEntries();
    }
  } else {
    maxEvent = runLoader->GetNumberOfEvents();
  }
  if (nEvents < 0) nEvents = maxEvent - firstEvent;
  Int_t lastEvent = firstEvent + nEvents;
  if (lastEvent > maxEvent) lastEvent = maxEvent;

  // loop over the events
  for (Int_t iEvent = firstEvent; iEvent < lastEvent; iEvent++) {

    runLoader->GetEvent(iEvent);

    // run the cluster finder
    if (makeClusters) {
      if (!tpcLoader->TreeR()) tpcLoader->MakeRecPointsContainer();
      clusterer->SetOutput(tpcLoader->TreeR());
      if (fileNameRaw) {
	AliRawReaderRoot rawReader(fileNameRaw, iEvent);
	clusterer->Digits2Clusters(&rawReader);
      } else {
	clusterer->SetInput(tpcLoader->TreeD());
	clusterer->Digits2Clusters();
      }
      tpcLoader->WriteRecPoints("OVERWRITE");
    }

    // run the track finder
    if (makeTracks) {
      tracker = new AliTPCtrackerMI(param);
      tracker->Clusters2Tracks();
      delete tracker;
    }
  }

  if (tracker) delete tracker;
  if (clusterer) delete clusterer;

  return kTRUE;
}

/*
////////////////////////////////////////////////////////////////////////
Int_t TPCRefitInward(Int_t nEvents, Int_t firstEvent,
		     const char* fileNameClusters,
		     const char* fileNameTracks,
		     const char* fileNameTRDTracks,
		     const char* fileNameRefittedTracks)
{
  Int_t rc = 0;
  const Char_t *name="TPCRefitInward";
  if (gDEBUG>1) cout<<name<<" starts"<<endl;
  if (gDEBUG>1) gBenchmark->Start(name);
  TFile *fileClusters = TFile::Open(fileNameClusters);
  TFile *fileTracks = TFile::Open(fileNameTracks);
  TFile *fileTRDTracks = TFile::Open(fileNameTRDTracks);
  TFile *fileRefittedTracks = TFile::Open(fileNameRefittedTracks, "recreate");

  AliTPCParam* paramTPC = AliTPC::LoadTPCParam(fileClusters);
  if (!paramTPC) return 1;

  for (Int_t iEvent = firstEvent; iEvent < firstEvent+nEvents; iEvent++){
    if (gDEBUG > 2) cout<<"TPCRefitInward: event "<<iEvent<<endl;
    AliTPCtracker *tracker = new AliTPCtracker(paramTPC);
    tracker->SetEventNumber(iEvent);
    fileClusters->cd();
    rc = tracker->RefitInward(fileTRDTracks, fileTracks, fileRefittedTracks);
    delete tracker;
    if (rc) return rc;
  }

  fileClusters->Close();
  fileTracks->Close();
  fileTRDTracks->Close();
  fileRefittedTracks->Close();
  delete fileClusters;
  delete fileTracks;
  delete fileTRDTracks;
  delete fileRefittedTracks;
  if (gDEBUG>1) gBenchmark->Show(name);
  return rc;
}
*/
