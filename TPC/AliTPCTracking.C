////////////////////////////////////////////////////////////////////////
//
// AliTPCTracking.C 
//
// date: 22.08.2002
// author: Jiri Chudoba based on code of Jourij Belikov
// version: 1.0
// description: 
//      reconstructs of tracks in TPC inthe following steps:
//         TPC cluster finding
//         TPC track finding
// input parameters: 
//        Int_t nEvents      ... nr of events to process
//        Int_t firstEventNr ... first event number (starts from 0)
//        const char* fileNameHits ... name of file with hits
//        const char* fileNameDigits .. name of file with TPC digits
//        const char* fileNameClusters .. name of file with TPC clusters (output)
//        const char* fileNameTracks .. name of file with TPC tracks (output)
//
//        default file names correspond to pp production (2002-04)
//
// History:
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
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TBenchmark.h"
#include "AliTPCtracker.h"
#include "AliTPCtrackerMI.h"
#include "AliTPCclusterer.h"
#include "AliTPCclustererMI.h"
#include "AliTPC.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliTracker.h"

#endif

Int_t gDEBUG = 2;

Int_t AliTPCTracking(Int_t nEvents=1, Int_t firstEvent=0,
		     const char* fileNameHits="galice.root",
		     const char* fileNameDigits="tpc.digits.root",
		     const char* fileNameClusters="tpc.clusters.root",
		     const char* fileNameTracks="tpc.tracks.root",
		     Bool_t versionMI = kFALSE);

Int_t TPCFindClusters(Int_t nEvents=1, Int_t firstEvent=0,
		      const char* fileNameDigits="tpc.digits.root", 
		      const char* fileNameClusters="tpc.clusters.root",
		      Bool_t versionMI = kFALSE);
Int_t TPCFindClusters(Int_t nEvents, Int_t firstEvent,
		      TFile* fileDigits, TFile* fileClusters, 
		      AliTPCParam* paramTPC=0);
Int_t TPCFindClustersMI(Int_t nEvents, Int_t firstEvent,
			TFile* fileDigits, TFile* fileClusters, 
			AliTPCParam* paramTPC=0);
Int_t TPCFindTracks(Int_t nEvents=1, Int_t firstEvent=0,
		    const char* fileNameClusters="tpc.clusters.root",
		    const char* fileNameTracks="tpc.tracks.root",
		    Bool_t versionMI = kFALSE);
Int_t TPCFindTracks(Int_t nEvents, Int_t firstEvent,
		    TFile* fileClusters, TFile* fileTracks,
		    AliTPCParam* paramTPC=0);
Int_t TPCFindTracksMI(Int_t nEvents, Int_t firstEvent,
		      TFile* fileClusters, TFile* fileTracks,
		      AliTPCParam* paramTPC=0);

void FindVertex(Int_t iEvent, Double_t *vertex);
void PrintVertex(TArrayF &primaryVertex);

////////////////////////////////////////////////////////////////////////
Int_t AliTPCTracking( Int_t nEvents, Int_t firstEvent,
		      const char* fileNameHits,
		      const char* fileNameDigits,
		      const char* fileNameClusters,
		      const char* fileNameTracks,
		      Bool_t versionMI) {

  AliTracker::SetFieldFactor(fileNameHits,kFALSE);

// ********** Find TPC clusters *********** //
  if (fileNameDigits && fileNameClusters) {
    if (TPCFindClusters(nEvents,firstEvent,fileNameDigits,fileNameClusters,versionMI)) {
      cerr<<"Failed to get TPC clusters: !\n";
      return 1;
    }
  }      

// ********** Find TPC tracks *********** //
  if (fileNameClusters && fileNameTracks) {
    if (TPCFindTracks(nEvents,firstEvent,fileNameClusters,fileNameTracks,versionMI)) {
      cerr<<"Failed to get TPC tracks !\n";
      return 2;
    }
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////
Int_t TPCFindClusters(Int_t nEvents, Int_t firstEvent,
		      const char* fileNameDigits, const char* fileNameClusters,
		      Bool_t versionMI) {
  
  Int_t rc;
  const Char_t *name="TPCFindClusters";
  if (gDEBUG>1) cout<<name<<" starts...\n";
  if (gDEBUG>1) gBenchmark->Start(name);
  TFile *fileClusters = TFile::Open(fileNameClusters,"recreate");
  TFile *fileDigits = TFile::Open(fileNameDigits);
  if (!fileDigits->IsOpen()) {
    cerr<<"Cannnot open "<<fileNameDigits<<" !\n"; 
    return 1;
  }
  if (!fileClusters->IsOpen()) {
    cerr<<"Cannnot open "<<fileNameClusters<<" !\n"; 
    return 1;
  }

  if (versionMI) {
    rc = TPCFindClustersMI(nEvents,firstEvent,fileDigits,fileClusters);
  } else {
    rc = TPCFindClusters(nEvents,firstEvent,fileDigits,fileClusters);
  }

  fileDigits->Close();
  fileClusters->Close();
  delete fileDigits;
  delete fileClusters;
  if (gDEBUG>1) gBenchmark->Stop(name);
  if (gDEBUG>1) gBenchmark->Show(name);

  return rc;
}
////////////////////////////////////////////////////////////////////////
Int_t TPCFindClusters(Int_t nEvents, Int_t firstEvent,
		      TFile* fileDigits, TFile* fileClusters, 
		      AliTPCParam* paramTPC) {

  fileDigits->cd();
  if (!paramTPC) paramTPC = AliTPC::LoadTPCParam(fileDigits);
  if (!paramTPC) return 1;

  for (Int_t iEvent = firstEvent; iEvent < firstEvent+nEvents; iEvent++){
    if (gDEBUG > 2) cout<<"TPCFindClusters: event "<<iEvent<<endl;
    AliTPCclusterer::Digits2Clusters(paramTPC, fileClusters, iEvent);
  }
  return 0;
}
////////////////////////////////////////////////////////////////////////
Int_t TPCFindClustersMI(Int_t nEvents, Int_t firstEvent,
			TFile* fileDigits, TFile* fileClusters, 
			AliTPCParam* paramTPC) {

  fileDigits->cd();
  if (!paramTPC) paramTPC = AliTPC::LoadTPCParam(fileDigits);
  if (!paramTPC) return 1;

  AliTPCclustererMI clusterer;
  for (Int_t iEvent = firstEvent; iEvent < firstEvent+nEvents; iEvent++){
    if (gDEBUG > 2) cout<<"TPCFindClustersMI: event "<<iEvent<<endl;
    char treeName[100];
    sprintf(treeName, "TreeD_75x40_100x60_150x60_%d", iEvent);
    TTree* input = (TTree*) fileDigits->Get(treeName);
    fileClusters->cd();
    sprintf(treeName, "TreeC_TPC_%d", iEvent);
    TTree* output = new TTree(treeName, treeName); 
    clusterer.SetInput(input);
    clusterer.SetOutput(output);
    clusterer.Digits2Clusters(paramTPC, iEvent);
  }

  return 0;
}
////////////////////////////////////////////////////////////////////////
Int_t TPCFindTracks(Int_t nEvents, Int_t firstEvent,
		    const char* fileNameClusters, const char* fileNameTracks,
		    Bool_t versionMI) {

  Int_t rc = 0;
  const Char_t *name="TPCFindTracks";
  if (gDEBUG>1) cout<<name<<" starts"<<endl;
  if (gDEBUG>1) gBenchmark->Start(name);
  TFile *fileTracks = TFile::Open(fileNameTracks,"recreate");
  TFile *fileClusters =TFile::Open(fileNameClusters);

  if (versionMI) {
    rc = TPCFindTracksMI(nEvents, firstEvent, fileClusters, fileTracks);
  } else {
    rc = TPCFindTracks(nEvents, firstEvent, fileClusters, fileTracks);
  }

  fileClusters->Close();
  fileTracks->Close();
  delete fileClusters;
  delete fileTracks;
  if (gDEBUG>1) gBenchmark->Stop(name);
  if (gDEBUG>1) gBenchmark->Show(name);
  return rc;

}
////////////////////////////////////////////////////////////////////////
Int_t TPCFindTracks(Int_t nEvents, Int_t firstEvent,
		    TFile *fileClusters, TFile * fileTracks,
		    AliTPCParam* paramTPC) {

  Int_t rc = 0;
  if (!paramTPC) paramTPC = AliTPC::LoadTPCParam(fileClusters);
  if (!paramTPC) return 1;

  for (Int_t iEvent = firstEvent; iEvent < firstEvent+nEvents; iEvent++){
    if (gDEBUG > 2) cout<<"TPCFindTracks: event "<<iEvent<<endl;
    AliTPCtracker *tracker = new AliTPCtracker(paramTPC);
    tracker->SetEventNumber(iEvent);
    Double_t vertex[3];
    FindVertex(iEvent,vertex);
    tracker->SetVertex(vertex);
    fileClusters->cd();
    rc = tracker->Clusters2Tracks(0,fileTracks);
    delete tracker;
  }
  return rc;
}
////////////////////////////////////////////////////////////////////////
Int_t TPCFindTracksMI(Int_t nEvents, Int_t firstEvent,
		      TFile *fileClusters, TFile * fileTracks,
		      AliTPCParam* paramTPC) {

  Int_t rc = 0;
  if (!paramTPC) paramTPC = AliTPC::LoadTPCParam(fileClusters);
  if (!paramTPC) return 1;

  for (Int_t iEvent = firstEvent; iEvent < firstEvent+nEvents; iEvent++){
    if (gDEBUG > 2) cout<<"TPCFindTracksMI: event "<<iEvent<<endl;
    AliTPCtrackerMI tracker(paramTPC, iEvent);
    rc = tracker.Clusters2Tracks(0, fileTracks);
  }
  return rc;
}

////////////////////////////////////////////////////////////////////////
void FindVertex(Int_t eventNr, Double_t *vertex) {

  vertex[0] = vertex[1] = vertex[2] = 0.;
  if (!gAlice) {
    cerr<<"gAlice was not found! Using vertex position (0,0,0).\n";
    return;
  }
  
  gAlice->GetEvent(eventNr);
  AliHeader *header = gAlice->GetHeader();
  if (!header) {
    cerr<<"header was not found!\n";
    return;
  } 
  AliGenEventHeader* genEventHeader = header->GenEventHeader();
  if (!genEventHeader) {
    cerr<<"AliGenEventHeader was not found!\n";
    return;
  } 

  TArrayF primaryVertex(3);
  genEventHeader->PrimaryVertex(primaryVertex);
  PrintVertex(primaryVertex);
  vertex[0] = static_cast<Double_t>(primaryVertex[0]);
  vertex[1] = static_cast<Double_t>(primaryVertex[1]);
  vertex[2] = static_cast<Double_t>(primaryVertex[2]);
//  delete header;
  delete genEventHeader;
  return;
     
}
////////////////////////////////////////////////////////////////////////
void PrintVertex(TArrayF &primaryVertex) 
{
  cout <<"Vertex: "
       <<primaryVertex[0]<<" "
       <<primaryVertex[1]<<" "
       <<primaryVertex[2]<<" "<<endl;
} 
////////////////////////////////////////////////////////////////////////
