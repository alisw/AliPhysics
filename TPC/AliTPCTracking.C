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
//        Char_t* fileNameHits ... name of file with hits
//        Char_t* fileNameDigits .. name of file with TPC digits
//        Char_t* fileNameClusters .. name of file with TPC clusters (output)
//        Char_t* fileNameTracks .. name of file with TPC tracks (output)
//
//        default file names correspond to pp production (2002-04)
//
// History:
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
#include "TTree.h"
#include "TSystem.h"
#include "TArrayF.h"
#include "TPC/alles.h"
#include "TPC/AliTPCtracker.h"
#include "TPC/AliTPCclusterer.h"
#include "TPC/AliTPC.h"
#include "STEER/AliRun.h"
#include "STEER/AliHeader.h"
#include "STEER/AliGenEventHeader.h"
#include "STEER/AliMagF.h"
#include "STEER/AliTracker.h"

#endif

Int_t gDEBUG = 2;

Int_t TPCFindClusters(Int_t nEvents=1, Int_t firstEvent=0,
		      Char_t* fileNameDigits="rfio:galiceSDR.root", 
		      Char_t* fileNameClusters="tpc.clusters.root");
Int_t TPCFindClusters(Int_t nEvents, Int_t firstEvent,
		      TFile* fileDigits, TFile* fileClusters, 
		      AliTPCParam* paramTPC=0);
Int_t TPCFindTracks(Int_t nEvents=1, Int_t firstEvent=0,
		    Char_t* fileNameClusters="tpc.clusters.root",
		    Char_t* fileNameTracks="tpc.tracks.root");
Int_t TPCFindTracks(Int_t nEvents, Int_t firstEvent,
		    TFile* fileClusters, TFile* fileTracks,
		    AliTPCParam* paramTPC=0);

void FindVertex(Int_t iEvent, Double_t *vertex);
void PrintVertex(TArrayF &primaryVertex);

Int_t AliTPCTracking(Int_t nEvents=1, Int_t firstEvent=0,
		     Char_t* fileNameHits="rfio:galice.root",
		     Char_t* fileNameDigits="rfio:galiceSDR.root",
		     Char_t* fileNameClusters="tpc.clusters.root",
		     Char_t* fileNameTracks="tpc.tracks.root");

////////////////////////////////////////////////////////////////////////
Int_t AliTPCTracking( Int_t nEvents, Int_t firstEvent,
		      Char_t* fileNameHits,
		      Char_t* fileNameDigits,
		      Char_t* fileNameClusters,
		      Char_t* fileNameTracks) {

  AliTracker::SetFieldFactor(fileNameHits,kFALSE);

// ********** Find TPC clusters *********** //
  if (TPCFindClusters(nEvents,firstEvent,fileNameDigits,fileNameClusters)) {
    cerr<<"Failed to get TPC clusters: !\n";
    return 1;
  }      

// ********** Find TPC tracks *********** //
  if (TPCFindTracks(nEvents,firstEvent,fileNameClusters,fileNameTracks)) {
    cerr<<"Failed to get TPC tracks !\n";
    return 2;
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////
Int_t TPCFindClusters(Int_t nEvents, Int_t firstEvent,
		      Char_t* fileNameDigits, Char_t* fileNameClusters) {
  
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

  rc = TPCFindClusters(nEvents,firstEvent,fileDigits,fileClusters);

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
Int_t TPCFindTracks(Int_t nEvents, Int_t firstEvent,
		    Char_t* fileNameClusters, Char_t* fileNameTracks) {

  Int_t rc = 0;
  const Char_t *name="TPCFindTracks";
  if (gDEBUG>1) cout<<name<<" starts"<<endl;
  if (gDEBUG>1) gBenchmark->Start(name);
  TFile *fileTracks = TFile::Open(fileNameTracks,"recreate");
  TFile *fileClusters =TFile::Open(fileNameClusters);

  rc = TPCFindTracks(nEvents, firstEvent, fileClusters, fileTracks);

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
  exit;
} 
////////////////////////////////////////////////////////////////////////
