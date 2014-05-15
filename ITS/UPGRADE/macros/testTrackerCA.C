#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>
#include <TStopwatch.h>

#include "AliTracker.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPManager.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliMagF.h"
#include "AliESDEvent.h"
#include "AliITSURecoParam.h"
#include "AliITSUReconstructor.h"
#include "AliITSUTrackerSA.h"
#endif

extern TSystem *gSystem;

const AliESDVertex *SetMCvertex(const AliRunLoader *rl, AliTracker *tr);

void testTrackerCA() {
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
  gSystem->Load("libITSUpgradeRec");

  // TGeoGlobalMagField::Instance()->
  //   SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));
  // // TGeoManager::Import("geometry.root");
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s",gSystem->pwd()));
  man->SetSpecificStorage("ITS/Align/Data",
			  Form("local://%s",gSystem->pwd()));
  man->SetSpecificStorage("ITS/Calib/RecoParam",
			  Form("local://%s",gSystem->pwd()));
  man->SetRun(0);
  if ( !TGeoGlobalMagField::Instance()->GetField() ) {
    printf("Loading field map...\n");
    AliGRPManager grpMan;
    if( !grpMan.ReadGRPEntry() ) { 
      printf("Cannot get GRP entry\n"); 
    }
    if( !grpMan.SetMagField() ) { 
      printf("Problem with magnetic field setup\n"); 
    }
  }
  AliGeomManager::LoadGeometry("geometry.root");
  //
  //
  AliCDBEntry* ent = man->Get("ITS/Calib/RecoParam");
  AliITSURecoParam* par = (AliITSURecoParam*)((TObjArray*)ent->GetObject())->At(1);
  //
  AliITSUReconstructor *rec = new AliITSUReconstructor();
  rec->SetRecoParam(par);
  //
  rec->Init();
  AliITSUTrackerSA *tracker = new AliITSUTrackerSA();
  tracker->Init(rec);


  TFile *esdFile=TFile::Open("AliESDs.root","recreate");
  TTree *esdTree = new TTree("esdTree", "Tree with ESD objects");
  AliESDEvent *esd=new AliESDEvent();
  esd->CreateStdContent();
  esd->WriteToTree(esdTree);
    
  TFile *clsFile=TFile::Open("ITS.RecPoints.root");
    
  AliRunLoader *rl = AliRunLoader::Open("galice.root","something");
  rl->LoadHeader();
  

  TStopwatch timer;    
  Int_t nEvents=1;//rl->GetNumberOfEvents();
  for (Int_t i=0; i<nEvents; i++) {
    cout<<"\nEvent number "<<i<<endl;
    rl->GetEvent(i);

    const AliESDVertex *vtx=SetMCvertex(rl,tracker);
    esd->SetPrimaryVertexSPD(vtx);              

    TTree *cTree=(TTree *)clsFile->Get(Form("Event%d/TreeR",i));
    tracker->LoadClusters(cTree);
    tracker->Clusters2Tracks(esd);
    //tracker->PropagateBack(esd);
    tracker->RefitInward(esd);
    tracker->UnloadClusters();

    Int_t n=esd->GetNumberOfTracks();
    cout << "Number of reconstructed tracks " << n << endl;
    for (Int_t t=0; t<n; t++) {
      AliESDtrack *track=esd->GetTrack(t);
      if (!track->RelateToVertex(vtx, tracker->GetBz(), 33)) continue;
      //Double_t r[3]; track->GetXYZ(r);
      //cout<<r[0]<<' '<<r[1]<<' '<<r[2]-vtx->GetZ()<<endl;
    }

    esdTree->Fill();
    esd->Reset();
    delete vtx;
  }
  timer.Stop(); timer.Print();
 
  delete tracker;
  //delete clsFile;
  esdFile->cd();
  esdTree->Write();
  delete esd;
  delete esdFile;
}

const AliESDVertex *SetMCvertex(const AliRunLoader *rl, AliTracker *tracker) {
  AliGenEventHeader *h=rl->GetHeader()->GenEventHeader();
  TArrayF vtx(3);
  h->PrimaryVertex(vtx);
  cout<<"Vertex "<<vtx[0]<<' '<<vtx[1]<<' '<<vtx[2]<<endl;
  Double_t xyz[]={vtx[0],vtx[1],vtx[2]};
  Double_t ers[]={2.,2.,2.};
  tracker->SetVertex(xyz,ers);
  AliESDVertex *vertex=new AliESDVertex(xyz,ers);
  return vertex;
}
