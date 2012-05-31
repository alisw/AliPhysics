#if !defined(__CINT__) || defined(__MAKECINT__)
#include<Riostream.h>
#include<AliRun.h>
#include<TTree.h>
#include<AliRunLoader.h>
#include<AliESD.h>
#include<TMath.h>
#include<AliHeader.h>
#include<AliGenEventHeader.h>
#include<AliVertexerTracks.h>
#include<AliITSVertexerTracks.h>
#include<AliVertex.h>
#include<AliESDVertex.h>
#include<TFile.h>
#endif

void TestVertexerTracks(){
  if (gAlice) {
    delete AliRunLoader::Instance();
    delete gAlice;
    gAlice = 0x0;
  }

  AliRunLoader *rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0) {
    cerr<<"Can not open session"<<endl;
    return;
  }
  
  rl->LoadgAlice();
  AliMagF * magf = gAlice->Field();
  AliTracker::SetFieldMap(magf,kTRUE);
  printf("MagneticField=%f\n",AliTracker::GetBz());

  rl->LoadHeader();
  rl->LoadKinematics();
  Int_t totev=rl->GetNumberOfEvents();
  cout<<"Number of events="<<totev<<endl;

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

  rl->GetEvent(e);
  AliESDVertex *vertESD = event->GetVertex();
  Double_t recVertex[3];
  AliGenEventHeader *header=rl->GetHeader()->GenEventHeader();
  TArrayF mcVertex(3);
  header->PrimaryVertex(mcVertex);
  cout<<"Primary vertex (MC) ";
  for(Int_t kk=0;kk<3;kk++)cout<<mcVertex[kk]<<" ";
  cout<<endl;
  vertESD->GetXYZ(recVertex);
  cout<<"Primary vertex (from SPD clusters) ";
  for(Int_t kk=0;kk<3;kk++)cout<<recVertex[kk]<<" ";
  cout<<endl;
  TObjArray *trkarr=new TObjArray();
  trkarr->Expand(3);

  TTree *trkTree = new TTree("TreeT","tracks");
  AliESDtrack *esdTrack = 0;
  trkTree->Branch("tracks","AliESDtrack",&esdTrack);

  esdTrack=(AliESDtrack*)event->GetTrack(2);
  trkarr->AddLast(esdTrack);
  trkTree->Fill();

  esdTrack=(AliESDtrack*)event->GetTrack(3);
  trkarr->AddLast(esdTrack);
  trkTree->Fill();

  esdTrack=(AliESDtrack*)event->GetTrack(4);
  trkarr->AddLast(esdTrack);
  trkTree->Fill();

  // BEWARE: the method VertexForSelectedTracks returns the address of the data member fVert 
  // which is always the same for all calls to VertexForSelectedTracks.

  printf("\n*****************************\n");
  printf("Call VertexForSelectedTracks with different arguments\n");
  printf("*****************************\n");
  AliVertexerTracks *vtx=new AliVertexerTracks();
  AliVertex *vert=vtx->VertexForSelectedTracks(trkarr);
  vert->Print();
  AliVertex *vert2=vtx->VertexForSelectedTracks(trkTree);
  vert->Print();
  vert2->Print();

  printf("\n*****************************\n");
  printf("Use different finder algorithms\n");
  printf("*****************************\n");
  vtx->SetFinderAlgorithm(1);
  AliVertex *v1=vtx->VertexForSelectedTracks(trkarr);
  v1->Print();
  vtx->SetFinderAlgorithm(2);
  AliVertex *v2=vtx->VertexForSelectedTracks(trkarr);
  v2->Print();
  vtx->SetFinderAlgorithm(3);
  AliVertex *v3=vtx->VertexForSelectedTracks(trkarr);
  v3->Print();
  vtx->SetFinderAlgorithm(4);
  AliVertex *v4=vtx->VertexForSelectedTracks(trkarr);
  v4->Print();
  vtx->SetFinderAlgorithm(5);
  AliVertex *v5=vtx->VertexForSelectedTracks(trkarr);
  v5->Print();
  delete vtx;

  printf("\n*****************************\n");
  printf("Call AliITSVertexerTracks::VertexForSelected Tracks\n");
  printf("*****************************\n");
  Int_t trkPos[3]={0,1,2};
  AliITSVertexerTracks *itsvtx=new AliITSVertexerTracks();
  itsvtx->SetDebug(0);
  AliVertex *vert3=itsvtx->VertexForSelectedTracks(event,3,trkPos,1);
  vert3->Print();
  delete itsvtx;
  
  printf("\n*****************************\n");
  printf("Find Primary Vertex\n");
  printf("*****************************\n");
  AliITSVertexerTracks *itsvtx2=new AliITSVertexerTracks();
  itsvtx2->SetDebug(0);
  AliESDVertex* evert=itsvtx2->FindPrimaryVertexForCurrentEvent(event);
  evert->Print();
  delete itsvtx2;
}



