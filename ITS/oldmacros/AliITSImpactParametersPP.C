/****************************************************************************
 *                                                                          *
 * This macro determines the impact parameter of each track in pp event     *
 * with respect to the primaty vertex position reconstructed using all      *
 * the other tracks in the event. This procedure is necessary to avoid      *
 * biases in the impact parameter estimate.                                 *
 *                                                                          *
 * Output:  TNtuple (ntd0) on file ImpactParameters.root                    *
 *          There is one entry per track. The elements of the ntuple are:   *
 *          1) Event: Event number                                          *
 *          2) nTraks: Number of tracks for this event                      *
 *          3) pt: Transverse momentum for the current track                *
 *          4) d0rphi: transverse impact parameter of the current track     *
 *          5) vx:     x coordinate of the vertex                           *
 *          6) vy:     y coordinate of the vertex                           *
 * origin: A.Dainese    andrea.dainese@pd.infn.it                           * 
 ****************************************************************************/ 
#if !defined(__CINT__) || defined(__MAKECINT__)
//-- --- standard headers------------- 
#include <Riostream.h>
//--------Root headers ---------------
#include <TSystem.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TObject.h>
#include <TVector3.h>
#include <TTree.h>
#include <TParticle.h>
#include <TArray.h>
//----- AliRoot headers ---------------
#include "alles.h"
#include "AliRun.h"
#include "AliKalmanTrack.h"
#include "AliITStrackV2.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliV0vertex.h"
#include "AliV0vertexer.h"
#include "AliITSVertex.h"
#include "AliITSVertexer.h"
#include "AliITSVertexerTracks.h"
#endif
//-------------------------------------

// field (T)
const Double_t kBz = 0.4;

// this function ckecks file existence
Bool_t GetInputFile();
void   MoveTracks(TTree& itsTree,TObjArray& array);


void   AliITSImpactParametersPP(Int_t evFirst=0,Int_t evLast=0) {

  const Char_t *name="AliITSImpactParametersPP";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);  

  AliKalmanTrack::SetConvConst(100/0.299792458/kBz);


  // check existence of input file
  if(!GetInputFile()) { cerr<<"No tracks file found"<<endl; return; }

  TString outName("ImpactParameters.root");

  // primary vertex
  Double_t v1[3] = {0.,0.,0,};
  
  Int_t    ev,itsEntries,i;
  Int_t    nTotEv=0;
  Double_t pt,d0rphi;
  Char_t   trksName[100];
  AliITStrackV2 *itstrack = 0;

  // output ntuple
  TNtuple *ntd0 = new TNtuple("ntd0","ntd0","Event:nTrks:pt:d0rphi:vx:vy");

  // create the AliITSVertexerTracks object
  AliITSVertexerTracks *vertexer1 = new AliITSVertexerTracks;
  vertexer1->SetMinTracks(3);
  vertexer1->SetDebug(0);
  Int_t skipped[2];

  // Open file with ITS tracks
  TFile* itstrks = TFile::Open("AliITStracksV2.root");


  // loop on events in file
  for(ev=evFirst; ev<=evLast; ev++) {
    printf(" --- Processing event  %d ---\n",ev);
    sprintf(trksName,"TreeT_ITS_%d",ev);

    // tracks from ITS
    TTree *itsTree  =(TTree*)itstrks->Get(trksName);
    if(!itsTree) continue;
    itsEntries = (Int_t)itsTree->GetEntries();
    printf("+++\n+++ Number of tracks in ITS: %d\n+++\n\n",itsEntries);

    TObjArray tArray(itsEntries);
    MoveTracks(*itsTree,tArray);
    
    // count the total number of events
    nTotEv++;

    // loop on tracks
    for(i=0; i<itsEntries; i++) {
      itstrack = (AliITStrackV2*)tArray.At(i);

      // pt
      pt = itstrack->Pt();

      // primary vertex from other tracks in event
      Bool_t goodVtx = kFALSE;
      vertexer1->SetVtxStart(0.,0.);
      skipped[0] = i;
      skipped[1] = -1;
      vertexer1->SetSkipTracks(1,skipped);
      AliITSVertex *vertex1 = 
        (AliITSVertex*)vertexer1->VertexOnTheFly(*itsTree); 
      if(vertex1->GetNContributors()>0) goodVtx = kTRUE; 
      vertex1->GetXYZ(v1);
      //vertex1->PrintStatus(); 
      // impact parameter w.r.t. v1      (in microns)
      d0rphi =  10000.*itstrack->GetD(v1[0],v1[1]);
      delete vertex1;

      // fill ntuple
      if (goodVtx) ntd0->Fill(ev,itsEntries,pt,d0rphi,v1[0],v1[1]);

      itstrack = 0;
    }   // loop on tracks

    delete itsTree;

  }    // loop on events in file


  printf("\n+++\n+++ Total number of events: %d\n+++\n",nTotEv);

  delete vertexer1;

  itstrks->Close();

  // store ntuple in a file
  TFile* outroot = new TFile(outName.Data(),"recreate");
  ntd0->Write();
  outroot->Close();
  delete outroot;


  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return;
}
//___________________________________________________________________________
Bool_t   GetInputFile() {
  TString itsName("AliITStracksV2.root");

  if(gSystem->AccessPathName(itsName.Data(),kFileExists)) return kFALSE;

  return kTRUE;
}
//___________________________________________________________________________
void MoveTracks(TTree& itsTree,TObjArray& array) {
//
// this function creates two TObjArrays with tracks
//
 
  Int_t entr = (Int_t)itsTree.GetEntries();

  // trasfer tracks from tree to array
  for(Int_t i=0; i<entr; i++) {

    AliITStrackV2 *t = new AliITStrackV2; 
    itsTree.SetBranchAddress("tracks",&t);

    itsTree.GetEvent(i);

    array.AddLast(t);

  }

  return;
}



