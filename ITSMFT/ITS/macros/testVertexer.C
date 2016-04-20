#if !defined(__CINT__) || defined(__MAKECINT__)
#include <time.h>
#include <sys/time.h>
#include "AliITSUVertexerFast.h"
#include <TROOT.h>
#include <Riostream.h>
#include <TSystem.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TBranch.h>
#include <TArrayF.h>
#include <TH1F.h>
#include <TParticle.h>
#include "AliESDVertex.h"
#include "AliRun.h"
#include "AliGenEventHeader.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliLoader.h"
#include "AliStack.h"
#include "AliGeomManager.h"
#include "AliITSUVertexer.h"
#include "AliITSUClusterPix.h"
#include "AliITSURecoDet.h"
#include "AliITSURecoLayer.h"
#include "AliITSUGeomTGeo.h"
#include "AliLog.h"
#include "AliESDEvent.h"
#endif

void testVertexer(int nev=-1,int evStart=0) {
  double elapsedTime[3];
  gROOT->SetStyle("Plain");
  gSystem->CompileMacro("AliITSUVertexerFast.cxx","kg");

  //AliLog::SetGlobalDebugLevel(5);

  gAlice=NULL;
  // phi (rad), z(cm), pair(cm), cluster(cm)
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  runLoader->LoadgAlice();

  gAlice = runLoader->GetAliRun();

  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadRecPoints();

  AliLoader *dl = runLoader->GetDetectorLoader("ITS");

  AliGeomManager::LoadGeometry("geometry.root");
  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE);
  Int_t nLayers = gm->GetNLayers();
  AliITSUClusterPix::SetGeom(gm);
  AliITSURecoDet *its = new AliITSURecoDet(gm, "ITSinterface");
  its->CreateClusterArrays();
 
  TTree * cluTree = 0x0;
  
  int nevTot = (Int_t)runLoader->GetNumberOfEvents();
  printf("N Events : %i \n",nevTot);
  evStart = evStart<nevTot ? evStart : nevTot-1;
  if (evStart<0) evStart = 0;
  //
  int lastEv = nev<0 ? nevTot : evStart+nev;
  if (lastEv > nevTot) lastEv = nevTot;
  //

  TFile* esdFile = TFile::Open("AliESDs.root");
  if (!esdFile || !esdFile->IsOpen()) {
    printf("Error in opening ESD file");
    return;
  }
  
  AliESDEvent * esd = new AliESDEvent;
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    printf("Error: no ESD tree found");
    return;
  }
  esd->ReadFromTree(tree);

  TArrayF mcVertex(3); 

  float r[3],dz;
  for (int i = 0; i < 3; ++i) {
  	AliITSURecoLayer* layer = its->GetLayerActive(i);
  	r[i] = layer->GetR();
  	dz = layer->GetZMax() - layer->GetZMin();
  }

  AliITSUVertexerFast *vert = new AliITSUVertexerFast(64,256,dz);
  vert->SetRadii(r[0],r[1],r[2]);

  for (Int_t iEvent = evStart; iEvent < lastEv; iEvent++) {
    runLoader->GetEvent(iEvent);
    printf("\n Event %i \n",iEvent);
    
    runLoader->GetHeader()->GenEventHeader()->PrimaryVertex(mcVertex);
  
    AliStack *stack = runLoader->Stack();
    cluTree=dl->TreeR();
    int nlr=0;
    TClonesArray *clusters[3];
    while(1) {
      TBranch* br = cluTree->GetBranch(Form("ITSRecPoints%d",nlr));
      if (!br) break;
      br->SetAddress(its->GetLayerActive(nlr)->GetClustersAddress());
      if(nlr < 3) {
        clusters[nlr] = its->GetLayerActive(nlr)->GetClusters();
        nlr++;
      }
      else break;
    }
    cluTree->GetEntry(0);
    struct timeval start, end;
    gettimeofday(&start, NULL);
    vert->LoadClusters(clusters);
    float xyz[3];
    vert->FindVertex(xyz);   
    gettimeofday(&end, NULL);
    elapsedTime[iEvent] = (end.tv_sec - start.tv_sec) * 1000.0;      // sec to ms
    elapsedTime[iEvent] += (end.tv_usec - start.tv_usec) / 1000.0;   // us to ms
  }
  
  cout << "AAA\t " <<elapsedTime[0] <<"\t"<<elapsedTime[1]<<"\t"<<elapsedTime[2]<<endl;
}
