////////////////////////////////////////////////////////////////////////
//
// AliITSandTPCrecoV2.C 
//
//
////////////////////////////////////////////////////////////////////////


/****************************************************************************
 * This macro is supposed to do reconstruction in the barrel ALICE trackers *
 * (Make sure you have TPC digits and ITS hits before using this macro !!!) *
 * It does the following steps                                              *
 *                   1) TPC cluster finding                                 *
 *                   2) TPC track finding                                   *
 *                   3) ITS cluster finding - fast recpoints -              *
 *                   4) ITS tracking                                        *
 * TPC part taken from AliTPCTracking.C - Jiri Chudoba                      *
 * The TPC part will be removed                                             *
 * M.Masera 30/09/2002 09:40                                                *
 * Last revision: 4/3/2003                                                  *
 *                                                                          *
 * Use:                                                                     *
 * The main function is AliITSandTPCrecoV2 All the input arguments have a   *
 * default value.                                                           *
 * nEvents: number of events to be processed.                               *
 *          If nEvents<0, all the available events are used starting from   *
 *          the first one                                                   *
 * firstEvent: is the first event to be analyzed.                           *
 *             It is set to 0 if nEvents<0 (default)                        *
 * fileNameHits: input root file with Kine and hits                         *
 * fileNameDigits: input root file with ITS and TPC digits and with ITS     *
 *                 recpoints. Only fast points are used at the moment.      *
 *                 An option to use slow points should be added             *
 * fileNameClusters: output root file with TPC clusters                     *
 * fileNameTracks: output root file with TPC tracks                         *
 * fileNameITSClusters: output root file with ITS clusters                  *
 * fileNameITSTracks: output root file with ITS tracks                      *
 *                                                                          *
 *   AliITSandTPCrecoV2(Int_t nEvents=-1, Int_t firstEvent=0,               * 
 *		     TString fileNameHits="galice.root",                    *
 *		     TString fileNameDigits="galiceSDR.root",               *
 *		     TString fileNameClusters="tpc.clusters.root",          *
 *		     TString fileNameTracks="tpc.tracks.root",              *
 *                   TString fileNameITSClusters="its.clusters.root",       *
 *                   TString fileNameITSTracks="its.tracks.root");          *
 *                                                                          *
 *  Global variables:
 *  
 * Int_t gDEBUG = 0;   // NO verbose printouts if ==0                       *
 * TArrayD *gzvtx = 0; // Z coordinate of primary vertices                  *
 * Int_t gCurrEve = 0; // current event number                              *
 * Bool_t gPPMode = kTRUE;  // it must be set to kTRUE for pp interactions  *
 * Bool_t gFastReco = kTRUE;  // it must set to kTRUE for fast reconstr.    *
 * const Double_t kgVertexNotFound = -100;  // the Z coord. of the          *
 *                                 primary vertex is set to this value      *
 *                                 when the vertexing fails                 *
 * Bool_t gUseNewClustersV2 = kFALSE;  // if kTRUE the AliITSclustereV2 class
 *                                 is used to produce V2 clusters.          *
 *                                 ITS recpoints are used otherwise         *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TBranch.h>
#include <TParticle.h>
#include <TRandom.h>
#include "alles.h"
#include "AliMagF.h"
#include "AliMagFCM.h"
#include "AliTPCtracker.h"
#include "AliTPCclusterer.h"
#include "Riostream.h"
#include "TString.h"
#include "TArrayD.h"
#include "TLine.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TParticle.h"
#include "TText.h"
#include "TSystem.h"
#include "TVector.h"
#include "TObjArray.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "AliITS.h"
#include "AliITSclustererV2.h"
#include "AliITSgeom.h"
#include "AliITSRecPoint.h"
#include "AliITSclusterV2.h"
#include "AliITStrackerV2.h"
#include "AliITSVertexerPPZ.h"
#include "AliITSVertexerIons.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "TArrayF.h"
#include "AliGenPythiaEventHeader.h"

#endif


Int_t gDEBUG = 0;
TArrayD *gzvtx = 0;
Int_t gCurrEve = 0;
Bool_t gPPMode = kTRUE;
Bool_t gFastReco = kTRUE;
Bool_t gUseNewClustersV2 = kFALSE; 
const Double_t kgVertexNotFound = -100.;


Int_t TPCFindClusters(Int_t nEvents=1, Int_t firstEvent=0,
		      TString fileNameDigits="galiceSDR.root", 
		      TString fileNameClusters="tpc.clusters.root");
Int_t TPCFindClusters(Int_t nEvents, Int_t firstEvent,
		      TFile* fileDigits, TFile* fileClusters, 
		      AliTPCParam* paramTPC=0);
Int_t TPCFindTracks(Int_t nEvents=1, Int_t firstEvent=0,
		    TString fileNameClusters="tpc.clusters.root",
		    TString fileNameTracks="tpc.tracks.root");
Int_t TPCFindTracks(Int_t nEvents, Int_t firstEvent,
		    TFile* fileClusters, TFile* fileTracks,
		    AliTPCParam* paramTPC=0);

Int_t AliTPCTracking(Int_t nEvents=1, Int_t firstEvent=0,
		     TString fileNameHits="galice.root",
		     TString fileNameDigits="galiceSDR.root",
		     TString fileNameClusters="tpc.clusters.root",
		     TString fileNameTracks="tpc.tracks.root");
Int_t AliITSandTPCrecoV2(Int_t nEvents=-1, Int_t firstEvent=0,
			TString fileNameHits="galice.root",
			TString fileNameDigits="galiceSDR.root",
			TString fileNameClusters="tpc.clusters.root",
			TString fileNameTracks="tpc.tracks.root",
			TString fileNameITSClusters="its.clusters.root",
			TString fileNameITSTracks="its.tracks.root");
Int_t ITSFindClusters(TString inname, TString recpointFileName, TString fileNameITSClusters, Int_t n, Int_t first);
Int_t ITSFindTracks(TString fileNameTracks, TString fileNameITSClusters, TString fileNameITSTracks, Int_t n, Int_t first, Int_t vc, Int_t vc2);
void FindZV(Int_t nEvents, Int_t first, TString FileNameHits, TString FileWithRecPoints);
Int_t DirectClusterFinder(Int_t eventn,Int_t first,TFile *in,TFile *out);
////////////////////////////////////////////////////////////////////////
Int_t AliITSandTPCrecoV2(Int_t nEvents, Int_t firstEvent,
			TString fileNameHits,
			TString fileNameDigits,
			TString fileNameClusters,
			TString fileNameTracks,
			TString fileNameITSClusters,
			TString fileNameITSTracks) {
  const Char_t *name="AliITSandTPCrecoV2";

  gBenchmark->Start(name);

  // Set the conversion constant for the Kalman filter
  // and set the gAlice pointer 
  Char_t *finame = (Char_t *)fileNameHits.Data();
  AliTracker::SetFieldFactor(finame,kFALSE);

  if(nEvents < 0 ) {
    nEvents = (Int_t)gAlice->TreeE()->GetEntries();
    firstEvent = 0;
    cout << "Processing events from " << firstEvent << " up to " << nEvents-1 <<endl;
  }

  if(fileNameDigits!=fileNameHits)gAlice->SetTreeDFileName(fileNameDigits.Data());
  gzvtx = new TArrayD(nEvents);

  // measure Z vertex coordinate for seeding
 
  FindZV(nEvents,firstEvent,fileNameHits,fileNameDigits); 
  
  // TPC tracking - abort macro execution if tracking fails
  if(AliTPCTracking(nEvents,firstEvent,fileNameHits,fileNameDigits,
		    fileNameClusters,fileNameTracks)) {
    cerr << "TPC tracking failed \n";
    return 1;
  }
 
  // ITS clustering
  if(ITSFindClusters(fileNameHits,fileNameDigits,fileNameITSClusters,nEvents,firstEvent)) {
    cerr << "ITS clustering failed \n";
    return 2;
  }

  Int_t firstpass = 1;  // 1=vert. constraint; 0=no vert. constr.; -1=skip pass
  Int_t secondpass = 0; // as above
  if(ITSFindTracks(fileNameTracks,fileNameITSClusters,fileNameITSTracks,nEvents,firstEvent,firstpass,secondpass)) {
    cerr << "ITS tracking failed \n";
    return 3;
  }


  gBenchmark->Stop(name);
  gBenchmark->Show(name);
  return 0;
}
////////////////////////////////////////////////////////////////////////
Int_t AliTPCTracking( Int_t nEvents, Int_t firstEvent,
		      TString fileNameHits,
		      TString fileNameDigits,
		      TString fileNameClusters,
		      TString fileNameTracks) {


  // ********** Find TPC clusters *********** //
  if (TPCFindClusters(nEvents,firstEvent,fileNameDigits,fileNameClusters)) {
    cerr<<"Failed to get TPC clusters: !\n";
    return 1;
  }      

  // ********** Find TPC tracks *********** //
  //   if (TPCFindTracks(TPCclsName,TPCtrkName,n)) {
  if (TPCFindTracks(nEvents,firstEvent,fileNameClusters,fileNameTracks)) {
    cerr<<"Failed to get TPC tracks !\n";
    return 2;
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////
Int_t TPCFindClusters(Int_t nEvents, Int_t firstEvent,
		      TString fileNameDigits, TString fileNameClusters) {
  
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
		    TString fileNameClusters, TString fileNameTracks) {

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
  Double_t vertex[3];
  for(Int_t j=0; j<3; j++){vertex[j]=0.;}
  AliTPCtracker *tracker = new AliTPCtracker(paramTPC);
  for (Int_t iEvent = firstEvent; iEvent < firstEvent+nEvents; iEvent++){
    if (gDEBUG > 2) cout<<"TPCFindTracks: event "<<iEvent<<endl;
    tracker->SetEventNumber(iEvent);
    if((*gzvtx)[iEvent] > -100){
      vertex[2] = (*gzvtx)[iEvent];
    }
    else {
      vertex[2]=0;
    }
    tracker->SetVertex(vertex);
    fileClusters->cd();
    rc = tracker->Clusters2Tracks(0,fileTracks);
  }
  delete tracker;
  return rc;
}


////////////////////////////////////////////////////////////////////////
Int_t ITSFindClusters(TString inname, TString recpointFileName, TString fileNameITSClusters, Int_t n, Int_t first) {
  Int_t rc=0;
  const Char_t *name="ITSFindClusters";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);
  TFile *out=TFile::Open(fileNameITSClusters,"recreate");
  TFile *in = (TFile*)gROOT->GetListOfFiles()->FindObject(inname.Data());
  if(gUseNewClustersV2){
    cout<<"Direct method\n";
    return DirectClusterFinder(n,first,in,out);
  }
  else {
    cout<<"Clusters through rec points \n";
  }
 
  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
  if (!ITS) { cerr<<"Can't get the ITS !\n"; return 1;}
  AliITSgeom *geom=ITS->GetITSgeom();
  AliITSclustererV2 *clusterer = new AliITSclustererV2(geom);
  out->cd();   
  geom->Write();
     
  Int_t ev;
  if(gFastReco && gDEBUG>0)cout <<"Using fast points\n";
  else cout<<"Using slow points\n";
  for (ev = first; ev<n; ev++){
    cout<<"ITSFindClusters: processing event "<<ev<<endl;
    in->cd();  
    gAlice->GetEvent(ev);
     
    TTree *pTree=gAlice->TreeR();
    if (!pTree){
      cerr << "ITSFindClusters: can't get TreeR\n";
      return 1;
    }
    TBranch *branch = 0;
    if (gFastReco) {
      branch = pTree->GetBranch("ITSRecPointsF");
    }
    else {
      branch = pTree->GetBranch("ITSRecPoints");
    }
    if (!branch) {
      cerr << "ITSFindClusters: can't get ITSRecPoints branch \n";
      return 2;
    }

    out->cd();
    TClonesArray *clusters=new TClonesArray("AliITSclusterV2",10000);
    char   cname[100];
    sprintf(cname,"TreeC_ITS_%d",ev);
  
    TTree *cTree=new TTree(cname,"ITS clusters");
    cTree->Branch("Clusters",&clusters);
 
    TClonesArray *points=new TClonesArray("AliITSRecPoint",10000);
    branch->SetAddress(&points);
     
    TClonesArray &cl=*clusters;
    Int_t nclusters=0;
    Int_t nentr=(Int_t)pTree->GetEntries();
    AliITSgeom *geom=ITS->GetITSgeom();
    Float_t lp[5];
    Double_t rot[9];
    Int_t lab[6];
    for (Int_t i=0; i<nentr; i++) {
      branch->GetEvent(i);
      clusterer->RecPoints2Clusters(points,i,clusters);
      cTree->Fill(); clusters->Delete();
      points->Clear();
    }
    cTree->Write();
    delete cTree; delete clusters; points->Delete(); delete points;
  }
  delete clusterer;

  out->Close();

  return rc;
}
////////////////////////////////////////////////////////////////////////
Int_t ITSFindTracks(TString inname, TString inname2, TString outname, Int_t n, Int_t first, Int_t vc, Int_t vc2) {
  Int_t rc=0;
  if(gPPMode){
    AliITStrackV2::SetSigmaYV(0.015);
    AliITStrackV2::SetSigmaZV(0.017);
  }
  TFile *out=TFile::Open(outname.Data(),"recreate");
  if (!out->IsOpen()) {
    cerr<<"Can't open file "<<outname.Data()<<endl; 
    return 1;
  }
  TFile *in =TFile::Open(inname.Data());
  if (!in->IsOpen()) {
    cerr<<"Can't open file "<<inname.Data()<<endl; 
    return 1;
  }
  TFile *in2 =TFile::Open(inname2.Data());
  if (!in2->IsOpen()) {
    cerr<<"Can't open file "<<inname2.Data()<<endl; 
    return 1;
  }

  AliITSgeom *geom=(AliITSgeom*)gFile->Get("AliITSgeom");
  if (!geom) { cerr<<"can't get ITS geometry !\n"; return 1;}
  Double_t vrtc[3];
  AliITStrackerV2 tracker(geom);
  for (Int_t i=first;i<n;i++){
    tracker.SetEventNumber(i);
    Double_t vtx=(*gzvtx)[i];
    cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    cout<<"ITSFindTracks -- event "<<i<<endl;
    cout<<"Computed vertex position: "<<vtx<<endl;
    in2->cd();
    vrtc[0]=0.;
    vrtc[1]=0.;
    vrtc[2]=vtx;
    if(vtx != kgVertexNotFound){
      tracker.SetVertex(vrtc);   // set vertex only if it was computed
    }
    else {
      vc = 0;  // use vertex contraint only if vertex info is available
      vc2 = -1;
    }

    tracker.SetupFirstPass(&vc);

    tracker.SetupSecondPass(&vc2);
   
    rc=tracker.Clusters2Tracks(in,out);
  }

  in->Close();
  in2->Close();
  out->Close();

  return rc;
}

//////////////////////////////////////////////////////////////////
void FindZV(Int_t nEvents, Int_t first, TString FileNameHits, TString FileWithRecPoints){
  TFile *in = (TFile*)gROOT->GetListOfFiles()->FindObject(FileNameHits.Data());
  TFile *fo = new TFile("AliITSVertices.root","recreate"); 
  if(FileNameHits!=FileWithRecPoints)gAlice->SetTreeRFileName(FileWithRecPoints);
  AliITSVertexer *findVert;
  if(gPPMode){
    findVert = new AliITSVertexerPPZ(in,fo);
  }
  else {
    findVert = new AliITSVertexerIons(in,fo);
  }
  Int_t last = first + nEvents - 1;
  AliITSVertex *vert = 0;
  for(Int_t i=first; i<=last; i++){
    gAlice->GetEvent(i);
    // The true Z coord. is fetched for comparison
    AliHeader *header = gAlice->GetHeader();
    AliGenEventHeader* genEventHeader = header->GenEventHeader();
    TArrayF primaryVertex(3);
    genEventHeader->PrimaryVertex(primaryVertex);
    vert = findVert->FindVertexForCurrentEvent(i);
    if(vert){
      Double_t pos[3];
      for(Int_t kk=0;kk<3;kk++)pos[kk]=(Double_t)primaryVertex[kk];
      vert->SetTruePos(pos);
      Double_t found = vert->GetZv();
      gzvtx->AddAt(found,i);
      findVert->WriteCurrentVertex();
    }
    else {
      gzvtx->AddAt(kgVertexNotFound,i);
    }
  }
  delete findVert;
  fo->Close();
  delete fo;
}



Int_t DirectClusterFinder(Int_t eventn,Int_t first,TFile *in,TFile* out){


   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   if (!ITS) { cerr<<"Can't find the ITS !\n"; return 3; }

   AliITSgeom *geom=ITS->GetITSgeom();

   geom->Write();

   TStopwatch timer;
   AliITSclustererV2 clusterer(geom);
   for (Int_t i=first; i<eventn; i++) {
       cerr<<"Processing event number: "<<i<<endl;
       gAlice->GetEvent(i);
       //ITS->MakeTreeC(); //To make the V1 cluster finders happy
       clusterer.SetEvent(i);
       if (gFastReco) clusterer.Digits2Clusters(in,out);
       else                 clusterer.Hits2Clusters(in,out);
   }
   timer.Stop(); timer.Print();

   out->Close();

   return 0;

}

