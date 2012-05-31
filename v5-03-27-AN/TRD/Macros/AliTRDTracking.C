////////////////////////////////////////////////////////////////////////
//
// AliTRDTracking.C 
//
// date: 06.02.2003
// author: Thomas Kuhr based on AliTPCTracking.C
// version: 1.0
// description: 
//      reconstructs of tracks in TRD in the following steps:
//         TRD cluster finding
//         TRD track finding
// input parameters: 
//        Int_t nEvents      ... nr of events to process
//        Int_t firstEventNr ... first event number (starts from 0)
//        const char* fileNameParam ... name of file with TPC parameters
//        const char* fileNameHits ... name of file with hits
//        const char* fileNameDigits .. name of file with TRD digits
//        const char* fileNameSeeds .. name of file with TPC track seeds
//        const char* fileNameClusters .. name of file with TRD clusters (output)
//        const char* fileNameTracks .. name of file with TRD tracks (output)
//
//        default file names correspond to pp production (2002-04)
//
////////////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "TFile.h"
#include "TBenchmark.h"
#include "AliTRD.h"
#include "AliTRDparameter.h"
#include "AliTRDtracker.h"
#include "AliTRDclusterizerV1.h"
#include "AliTPC.h"
#include "AliTPCParamSR.h"
#include "AliTPCtracker.h"
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITStrackerV2.h"
#include "AliRun.h"
#endif

Int_t gDEBUG = 2;

Int_t AliTRDTracking(Int_t nEvents=1, Int_t firstEvent=0,
		     const char* fileNameParam="trd.sdigits.root",
		     const char* fileNameDigits="trd.digits.root",
		     const char* fileNameSeed="tpc.tracks.root",
		     const char* fileNameClusters="trd.clusters.root",
		     const char* fileNameTracks="trd.tracks.root");

Bool_t ITSPropagateBack(Int_t nEvents=1, Int_t firstEvent=0,
			const char* fileNameClusters="its.clusters.root",
			const char* fileNameTracks="its.tracks.root",
			const char* fileNameBackTracks="its.back.tracks.root");
Bool_t TPCPropagateBack(Int_t nEvents=1, Int_t firstEvent=0,
			const char* fileNameClusters="tpc.clusters.root",
			const char* fileNameTracks="tpc.tracks.root",
			const char* fileNameITSBackTracks="its.back.tracks.root",
			const char* fileNameBackTracks="tpc.back.tracks.root");
Bool_t TRDFindClusters(Int_t nEvents=1, Int_t firstEvent=0,
		      const char* fileNameParam="trd.sdigits.root", 
		      const char* fileNameDigits="trd.digits.root", 
		      const char* fileNameClusters="trd.clusters.root");
Bool_t TRDFindTracks(Int_t nEvents=1, Int_t firstEvent=0,
		    const char* fileNameSeeds="tpc.tracks.root",
		    const char* fileNameClusters="trd.clusters.root",
		    const char* fileNameTracks="trd.tracks.root",
		     Bool_t inwards = kTRUE, Bool_t addSeeds = kTRUE);


////////////////////////////////////////////////////////////////////////
Int_t AliTRDTracking( Int_t nEvents, Int_t firstEvent,
		      const char* fileNameParam,
		      const char* fileNameDigits,
		      const char* fileNameSeeds,
		      const char* fileNameClusters,
		      const char* fileNameTracks) {

  AliTracker::SetFieldFactor(fileNameParam,kFALSE);

// ********** Find TRD clusters *********** //
  if (fileNameParam && fileNameDigits && fileNameClusters){
    if (!TRDFindClusters(nEvents,firstEvent,fileNameParam,fileNameDigits,fileNameClusters)) {
      cerr<<"Failed to get TRD clusters: !\n";
      return 1;
    }
  }      

// ********** Find TRD tracks *********** //
  if (fileNameSeeds && fileNameClusters && fileNameTracks) {
    if (!TRDFindTracks(nEvents,firstEvent,fileNameSeeds,fileNameClusters,fileNameTracks)) {
      cerr<<"Failed to get TRD tracks !\n";
      return 2;
    }
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
Bool_t ITSPropagateBack(Int_t nEvents, Int_t firstEvent,
			const char* fileNameClusters,
			const char* fileNameTracks,
			const char* fileNameBackTracks)
{
  const char *name="ITSPropagateBack";
  if (gDEBUG>1) cout<<name<<" starts...\n";
  if (gDEBUG>1) gBenchmark->Start(name);

  TFile* clustersFile = TFile::Open(fileNameClusters);
  if (!clustersFile->IsOpen()) {
    cerr<<"Cannot open "<<fileNameClusters<<" !\n"; 
    return kFALSE;
  }
  TFile* tracksFile = NULL;
  TFile* backTracksFile = NULL;
  if (strcmp(fileNameTracks, fileNameBackTracks) == 0) {
    tracksFile = backTracksFile = TFile::Open(fileNameTracks, "UPDATE");
  } else {
    backTracksFile = TFile::Open(fileNameBackTracks, "RECREATE");
    tracksFile = TFile::Open(fileNameTracks);
  }
  if (!tracksFile->IsOpen()) {
    cerr<<"Cannot open "<<fileNameTracks<<" !\n"; 
    return kFALSE;
  }

  gROOT->cd();
  AliITSgeom* geom = (AliITSgeom*) clustersFile->Get("AliITSgeom");
  AliITStrackerV2 tracker(geom);
  clustersFile->cd();

  for (Int_t iEvent = firstEvent; iEvent < firstEvent+nEvents; iEvent++) {
    tracker.SetEventNumber(iEvent);
    tracker.PropagateBack(tracksFile, backTracksFile);
  }

  if (tracksFile != backTracksFile) {
    backTracksFile->Close();
    delete backTracksFile;
  }
  tracksFile->Close();
  delete tracksFile;
  clustersFile->Close();
  delete clustersFile;

  if (gDEBUG>1) gBenchmark->Show(name);
  return kTRUE;
}
////////////////////////////////////////////////////////////////////////
Bool_t TPCPropagateBack(Int_t nEvents, Int_t firstEvent,
			const char* fileNameClusters,
			const char* fileNameTracks,
			const char* fileNameITSBackTracks,
			const char* fileNameBackTracks)
{
  const char *name="TPCPropagateBack";
  if (gDEBUG>1) cout<<name<<" starts...\n";
  if (gDEBUG>1) gBenchmark->Start(name);

  TFile* clustersFile = TFile::Open(fileNameClusters);
  if (!clustersFile->IsOpen()) {
    cerr<<"Cannot open "<<fileNameClusters<<" !\n"; 
    return kFALSE;
  }
  TFile* itsBackFile = NULL;
  if (fileNameITSBackTracks) {
    itsBackFile = TFile::Open(fileNameITSBackTracks);
    if (!itsBackFile->IsOpen()) {
      cerr<<"Cannot open "<<fileNameITSBackTracks<<" !\n"; 
      return kFALSE;
    }
  }
  TFile* tracksFile = NULL;
  TFile* backTracksFile = NULL;
  if (strcmp(fileNameTracks, fileNameBackTracks) == 0) {
    tracksFile = backTracksFile = TFile::Open(fileNameTracks, "UPDATE");
  } else {
    backTracksFile = TFile::Open(fileNameBackTracks, "RECREATE");
    tracksFile = TFile::Open(fileNameTracks);
  }
  if (!tracksFile->IsOpen()) {
    cerr<<"Cannot open "<<fileNameTracks<<" !\n"; 
    return kFALSE;
  }

  gROOT->cd();
  AliTPCParam* tpcParam = (AliTPCParam*) 
    clustersFile->Get("75x40_100x60_150x60");
  if (!tpcParam) tpcParam = new AliTPCParamSR;
  AliTPCtracker tracker(tpcParam);
  clustersFile->cd();

  for (Int_t iEvent = firstEvent; iEvent < firstEvent+nEvents; iEvent++) {
    tracker.SetEventNumber(iEvent);
    tracker.PropagateBack(tracksFile, itsBackFile, backTracksFile);
  }

  if (tracksFile != backTracksFile) {
    backTracksFile->Close();
    delete backTracksFile;
  }
  tracksFile->Close();
  delete tracksFile;
  clustersFile->Close();
  delete clustersFile;

  if (gDEBUG>1) gBenchmark->Show(name);
  return kTRUE;
}
////////////////////////////////////////////////////////////////////////
Bool_t TRDFindClusters(Int_t nEvents, Int_t firstEvent,
		      const char* fileNameParam, 
		      const char* fileNameDigits, const char* fileNameClusters) {
  
  const char *name="TRDFindClusters";
  if (gDEBUG>1) cout<<name<<" starts...\n";
  if (gDEBUG>1) gBenchmark->Start(name);

  TFile* paramFile = TFile::Open(fileNameParam);
  if (!paramFile->IsOpen()) {
    cerr<<"Cannot open "<<fileNameParam<<" !\n"; 
    return kFALSE;
  }
  AliTRDparameter* trdParam = 
    (AliTRDparameter*) paramFile->Get("TRDparameter"); 
  if (!trdParam) {
    cerr << "TRD parameters have not been found !\n"; 
    return kFALSE;
  }
  trdParam->ReInit();
  AliTRDclusterizerV1 clusterizer("clusterizer", "Clusterizer class"); 
  clusterizer.SetParameter(trdParam);

  for (Int_t iEvent = firstEvent; iEvent < firstEvent+nEvents; iEvent++) {
    clusterizer.Open(fileNameDigits, fileNameClusters, iEvent);
    clusterizer.ReadDigits();
    clusterizer.MakeClusters();
    clusterizer.WriteClusters(-1);
  }

  paramFile->Close();
  delete paramFile;

  if (gDEBUG>1) gBenchmark->Show(name);
  return kTRUE;
}
////////////////////////////////////////////////////////////////////////
Bool_t TRDFindTracks(Int_t nEvents, Int_t firstEvent,
		    const char* fileNameSeeds,
		     const char* fileNameClusters, const char* fileNameTracks,
		     Bool_t inwards, Bool_t addSeeds) {

  const char *name="TRDFindTracks";
  if (gDEBUG>1) cout<<name<<" starts...\n";
  if (gDEBUG>1) gBenchmark->Start(name);

  TFile* seedsFile = TFile::Open(fileNameSeeds);
  if (!seedsFile->IsOpen()) {
    cerr<<"Cannot open "<<fileNameSeeds<<" !\n"; 
    return kFALSE;
  }
  TFile* clustersFile = TFile::Open(fileNameClusters);
  if (!clustersFile->IsOpen()) {
    cerr<<"Cannot open "<<fileNameClusters<<" !\n"; 
    return kFALSE;
  }
  TFile* tracksFile = TFile::Open(fileNameTracks, "recreate");
  if (!tracksFile->IsOpen()) {
    cerr<<"Cannot open "<<fileNameTracks<<" !\n"; 
    return kFALSE;
  }

  AliTRDtracker tracker(clustersFile);
  if (addSeeds) tracker.SetAddTRDseeds();
  clustersFile->cd();  
  for (Int_t iEvent = firstEvent; iEvent < firstEvent+nEvents; iEvent++) {
    printf("Processing event %d\n", iEvent);
    tracker.SetEventNumber(iEvent);
    tracker.PropagateBack(seedsFile, tracksFile);
    if (inwards) tracker.Clusters2Tracks(tracksFile, tracksFile);
  }

  seedsFile->Close();
  delete seedsFile;
  clustersFile->Close();
  delete clustersFile;
  tracksFile->Close();
  delete tracksFile;

  if (gDEBUG>1) gBenchmark->Show(name);
  return kTRUE;
}

