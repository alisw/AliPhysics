////////////////////////////////////////////////////////////////////////
//
// AliITSTrackingV2.C 
//
// date: 18.03.2003
// author: Thomas Kuhr based on AliTPCTracking.C, AliITSFindClustersV2.C
//         and AliITSFindTracksV2.C
// version: 1.0
// description: 
//      reconstructs of tracks in ITS in the following steps:
//         ITS cluster finding
//         ITS track finding
// input parameters: 
//        Int_t nEvents      ... nr of events to process
//        Int_t firstEventNr ... first event number (starts from 0)
//        char* fileNameHits ... name of file with hits
//        char* fileNameITSDigits .. name of file with ITS digits
//        char* fileNameITSTracks .. name of file with TPC tracks
//        char* fileNameITSClusters .. name of file with ITS clusters (output)
//        char* fileNameITSTracks .. name of file with ITS tracks (output)
//
////////////////////////////////////////////////////////////////////////


#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "TFile.h"
#include "TBenchmark.h"
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSclustererV2.h"
#include "AliITStrackerV2.h"
#include "AliRun.h"
#endif

Int_t gDEBUG = 2;

Int_t AliITSTrackingV2(Int_t nEvents=1, Int_t firstEvent=0,
		       const char* fileNameHits="galice.root",
		       const char* fileNameITSDigits="its.digits.root",
		       const char* fileNameTPCTracks="tpc.tracks.root",
		       const char* fileNameITSClusters="its.clusters.root",
		       const char* fileNameITSTracks="its.tracks.root");
	     
Int_t ITSFindClustersV2(const char* fileNameITSDigits, 
			const char* fileNameITSClusters, 
			Int_t n, Int_t first);
Int_t ITSFindTracksV2(const char* fileNameTPCTracks, 
		      const char* fileNameITSClusters, 
		      const char* fileNameITSTracks, 
		      Int_t n, Int_t first);

////////////////////////////////////////////////////////////////////////
Int_t AliITSTrackingV2(Int_t nEvents, Int_t firstEvent,
		       const char* fileNameHits,
		       const char* fileNameITSDigits,
		       const char* fileNameTPCTracks,
		       const char* fileNameITSClusters,
		       const char* fileNameITSTracks) {

  AliTracker::SetFieldFactor(fileNameHits,kFALSE);

// find clusters

  if (fileNameITSDigits && fileNameITSClusters) {
    if(ITSFindClustersV2(fileNameITSDigits,fileNameITSClusters,nEvents,firstEvent)) {
      cerr << "ITS clustering failed \n";
      return 1;
    }
  }

// find tracks

  if (fileNameTPCTracks && fileNameITSClusters && fileNameITSTracks) {
    if(ITSFindTracksV2(fileNameTPCTracks,fileNameITSClusters,fileNameITSTracks,nEvents,firstEvent)) {
      cerr << "ITS tracking failed \n";
      return 2;
    }
  }

  return 0;
}
////////////////////////////////////////////////////////////////////////
Int_t ITSFindClustersV2(const char* fileNameITSDigits, const char* fileNameITSClusters, Int_t nEvents, Int_t firstEvent) {
//
// create ITS clusters, store them in the file fileNameITSClusters
// gAlice object must be in memory

  const char *name="ITSFindClustersV2";
  if (gDEBUG>1) cout<<name<<" starts...\n";
  if (gDEBUG>1) gBenchmark->Start(name);

  TFile *in =TFile::Open(fileNameITSDigits);
  if (!in->IsOpen()) {
    cerr<<"Can't open file "<<fileNameITSDigits<<endl; 
    return 1;
  }
  gAlice->SetTreeDFileName(fileNameITSDigits);
  TFile *out=TFile::Open(fileNameITSClusters,"recreate");
  if (!out->IsOpen()) {
    cerr<<"Can't open file "<<fileNameITSClusters<<endl; 
    return 1;
  }

  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
  if (!ITS) { cerr<<"Can't find the ITS !\n"; return 3; }

  AliITSgeom *geom=ITS->GetITSgeom();
  out->cd();
  geom->Write();
  gROOT->cd();

  AliITSclustererV2 clusterer(geom);
  for (Int_t iEvent = firstEvent; iEvent<firstEvent+nEvents; iEvent++){
    cout<<"ITSFindClusters: processing event "<<iEvent<<endl;
    gAlice->GetEvent(iEvent);
    in->cd();
    clusterer.SetEvent(iEvent);
    clusterer.Digits2Clusters(in,out); 
  }

  out->Close();
  in->Close();

  if (gDEBUG>1) gBenchmark->Show(name);
  return 0;
}
////////////////////////////////////////////////////////////////////////
Int_t ITSFindTracksV2(const char* fileNameTPCTracks, 
		      const char* fileNameITSClusters, 
		      const char* fileNameITSTracks, 
		      Int_t nEvents, Int_t first) {
  Int_t rc=0;
  const char *name="ITSFindTracksV2";
  if (gDEBUG>1) cout<<name<<" starts...\n";
  if (gDEBUG>1) gBenchmark->Start(name);

  TFile *out=TFile::Open(fileNameITSTracks,"recreate");
  if (!out->IsOpen()) {
    cerr<<"Can't open file "<<fileNameITSTracks<<endl; 
    return 1;
  }
  TFile *in =TFile::Open(fileNameTPCTracks);
  if (!in->IsOpen()) {
    cerr<<"Can't open file "<<fileNameTPCTracks<<endl; 
    return 1;
  }
  TFile *in2 =TFile::Open(fileNameITSClusters);
  if (!in2->IsOpen()) {
    cerr<<"Can't open file "<<fileNameITSClusters<<endl; 
    return 1;
  }

  AliITSgeom *geom=(AliITSgeom*)in2->Get("AliITSgeom");
  if (!geom) { cerr<<"can't get ITS geometry !\n"; return 1;}

  AliITStrackerV2 tracker(geom);
  for (Int_t iEvent=first;iEvent<first+nEvents;iEvent++){
    cout<<"ITSFindTracks -- event "<<iEvent<<endl;
    tracker.SetEventNumber(iEvent);
    rc=tracker.Clusters2Tracks(in,out);
  }

  in->Close();
  in2->Close();
  out->Close();

  if (gDEBUG>1) gBenchmark->Show(name);
  return rc;
}

