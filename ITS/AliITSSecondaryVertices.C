#if !defined(__CINT__) || defined(__MAKECINT__)
//-- --- standard headers------------- 
#include <Riostream.h>
//--------Root headers ---------------
#include <TSystem.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TObject.h>
#include <TTree.h>
#include <TChain.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCut.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
//----- AliRoot headers ---------------
//#include "alles.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliKalmanTrack.h"
#include "AliESDVertex.h"
#include "AliITSVertexer.h"
#include "AliITSVertexerTracks.h"
#include <AliHeader.h>
#include <AliGenEventHeader.h>
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliITStrackV2.h"
#include "AliITSSimpleVertex.h"
#include "AliITSStrLine.h"
#include "AliITSLoader.h"
#include <AliESD.h>
#include <AliStack.h>
//-------------------------------------
#endif  


void AliITSSecondaryVertices(Int_t ntracks, Int_t *pos) {
 
  // example of usige AliITSVertexerTracks to get a (secodnary) vertex from a set of tracks.

  Double_t field = 0.4;
 
  AliRunLoader *rl = AliRunLoader::Open("galice.root");
  if (!rl) {
    cerr<<"Can't start session !\n";
  }

  rl->LoadgAlice();
  if (rl->GetAliRun()){
    AliMagF * magf = gAlice->Field();
    AliKalmanTrack:: SetFieldMap(magf);
    AliKalmanTrack::SetUniformFieldTracking();
  }else {
    cerr<<"Can't get AliRun !\n";
    
  }

  field=rl->GetAliRun()->Field()->SolenoidField()/10.;
  rl->LoadHeader();
 
  
  TFile *ef=TFile::Open("AliESDs.root");
  if (!ef || !ef->IsOpen()) {cerr<<"Can't AliESDs.root !n";}
  AliESD* event = new AliESD;
  TTree* tree = (TTree*) ef->Get("esdTree");
  if (!tree) {cerr<<"no ESD tree foundn";};
  tree->SetBranchAddress("ESD", &event);


  Int_t e=0;
  tree->GetEvent(e);
  Int_t ntrk=event->GetNumberOfTracks();
  cout<<"Number of ESD tracks : "<<ntrk<<endl; 

  // Open input and output files
  TFile *outFile = TFile::Open("AliITSVertices.root","recreate");
    

  TTree *trkTree = new TTree("TreeT_ITS","its tracks");
  AliITStrackV2 *itstrack = 0;
  trkTree->Branch("tracks","AliITStrackV2",&itstrack,ntrk,0);
  Int_t pipe=0;
  for(Int_t i=0; i<ntrk; i++) {
    AliESDtrack *esdTrack = (AliESDtrack*)event->GetTrack(i);
    UInt_t status=esdTrack->GetStatus();
    UInt_t flags=AliESDtrack::kTPCin|AliESDtrack::kITSin;
    if ((status&AliESDtrack::kTPCin)==0)continue;
    if ((status&AliESDtrack::kITSin)==0)continue;    
    if ((status&AliESDtrack::kITSrefit)==0) continue;
    if ((status&flags)==status) pipe=1;

    itstrack = new AliITStrackV2(*esdTrack);
    if (pipe!=0) {
      itstrack->PropagateTo(3.,0.0028,65.19);
      itstrack->PropagateToVertex();  // This is not "exactly" the vertex 
    }
    Int_t nclus=itstrack->GetNumberOfClusters();

    if(nclus<6)continue;
    trkTree->Fill();
    itstrack = 0;
    // delete esdTrack; 
  }
  // Primary vertex
  Double_t vtx[3];
  event->GetVertex()->GetXYZ(vtx);
  cout<<"Primary Vertex:"<<endl;
  event->GetVertex()->PrintStatus();
  cout<<"\n\n\n"<<endl;
  

  // Create vertexer
  AliITSVertexerTracks *vertexer = new AliITSVertexerTracks();
  vertexer->SetDCAcut(0);
  vertexer->SetDebug(0);


  Double_t vertF1[3],vertF2[3],vertF3[3],vertF4[3],vertF5[3];
  Int_t ncombi1,ncombi2,ncombi3,ncombi4,ncombi5;
  Double_t sig1,sig2,sig3,sig4,sig5; 
  AliITSSimpleVertex *vert=0;

  
  // Calculate vertex from tracks in pos
  // BEWARE: the method VertexForSelectedTracks returns the address of the data member fSimpVert 
  // which is always the same for all calls to VertexForSelectedTracks.

  cout<<"\nStraight Line Min Dist finder with weights"<<endl;
  vert=(AliITSSimpleVertex*) vertexer->VertexForSelectedTracks(event,ntracks,pos,1);
  vert->Print();
  vert->GetXYZ(vertF1);
  ncombi1=vert->GetNContributors();
  sig1=vert->GetDispersion();


  cout<<"Straight Line Min Dist finder"<<endl;
  vert=(AliITSSimpleVertex*) vertexer->VertexForSelectedTracks(event,ntracks,pos,2);   
  vert->Print();
  vert->GetXYZ(vertF2);
  ncombi2=vert->GetNContributors();
  sig2=vert->GetDispersion();
  
  cout<<"Helix finder"<<endl;
  vert=(AliITSSimpleVertex*) vertexer->VertexForSelectedTracks(event,ntracks,pos,3);   
  vert->Print();
  vert->GetXYZ(vertF3);
  ncombi3=vert->GetNContributors();
  sig3=vert->GetDispersion();
 
  cout<<"Straight Line finder with weights"<<endl;
  vert=(AliITSSimpleVertex*) vertexer->VertexForSelectedTracks(event,ntracks,pos,4);   
  vert->Print();
  vert->GetXYZ(vertF4);
  ncombi4=vert->GetNContributors();
  sig4=vert->GetDispersion();
  
  cout<<"Straight Line finder with weights"<<endl;
  vert=(AliITSSimpleVertex*) vertexer->VertexForSelectedTracks(event,ntracks,pos,5);
  vert->Print();
  vert->GetXYZ(vertF5);
  ncombi5=vert->GetNContributors();
  sig5=vert->GetDispersion();
  

  delete vertexer;

  outFile->Close();
  return;
}
