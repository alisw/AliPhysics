#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>
#include <TStopwatch.h>
#include <TCanvas.h>
#include <TVirtualPad.h>

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
#include "AliITSUCATracker.h"
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
  AliITSUCATracker *tracker = new AliITSUCATracker(rec);
  tracker->SetPhiCut(0.3f);
  //tracker->SetThetaBin(0.03f);
  tracker->SetChi2Cut(3000);
  tracker->Init(rec);


  TFile *esdFile=TFile::Open("AliESDs.root","recreate");
  TTree *esdTree = new TTree("esdTree", "Tree with ESD objects");
  AliESDEvent *esd=new AliESDEvent();
  esd->CreateStdContent();
  esd->WriteToTree(esdTree);
    
  TFile *clsFile=TFile::Open("ITS.RecPoints.root");
    
  AliRunLoader *rl = AliRunLoader::Open("galice.root","something");
  rl->LoadHeader();

  TH1F fGHisto("ghisto","Good tracks #chi^{2};#chi^{2};Entries",100,0,400);
  TH1F fFHisto("fhisto","Fake tracks #chi^{2};#chi^{2};Entries",100,0,400);

  TStopwatch timer;    
  Int_t nEvents=rl->GetNumberOfEvents();
  for (Int_t i=0; i<nEvents; i++) {
    //for ( Int_t i=1; i<2; i++) {  
    cout<<"\nEvent number "<<i<<endl;
    rl->GetEvent(i);

    const AliESDVertex *vtx=SetMCvertex(rl,tracker);
    esd->SetPrimaryVertexSPD(vtx);              

    TTree *cTree=(TTree *)clsFile->Get(Form("Event%d/TreeR",i));
    tracker->LoadClusters(cTree);
    tracker->Clusters2Tracks(esd);
    tracker->PropagateBack(esd);
    tracker->RefitInward(esd);
    // fGHisto.Add(tracker->GetGHisto());
    // fFHisto.Add(tracker->GetFHisto());
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
  #ifdef _TUNING_
  TCanvas *cv1 = new TCanvas("chi2","chi2",1400,1000);
  cv1->Divide(2,2);
  for(int i=0;i<4;i++) {
    cv1->cd(i+1);
    cv1->cd(i+1)->SetLogy();
    tracker->fGoodCombChi2[i]->Draw();
    tracker->fFakeCombChi2[i]->SetLineColor(kRed);
    tracker->fFakeCombChi2[i]->Draw("same");
  }
  TCanvas *ccv1 = new TCanvas("deltan","DeltaN",1400,1000);
  ccv1->Divide(2,2);
  for(int i=0;i<4;i++) {
    ccv1->cd(i+1);
    ccv1->cd(i+1)->SetLogy();
    tracker->fGoodCombN[i]->Draw();
    tracker->fFakeCombN[i]->SetLineColor(kRed);
    tracker->fFakeCombN[i]->Draw("same");
  }
  TCanvas *cv2 = new TCanvas("chi2tracks","chi2tracks");
  cv2->cd();
  cv2->cd()->SetLogy();
  tracker->fGoodCombChi2[4]->Draw();
  tracker->fFakeCombChi2[4]->SetLineColor(kRed);
  tracker->fFakeCombChi2[4]->Draw("same");
  TCanvas *cv3 = new TCanvas();
  cv3->cd();
  tracker->fTanF->SetLineColor(kRed);
  tracker->fTanF->Draw("");
  tracker->fTan->Draw("same");
  TCanvas *cv33 = new TCanvas();
  cv33->cd();
  tracker->fPhiF->SetLineColor(kRed);
  tracker->fPhiF->Draw("");
  tracker->fPhi->Draw("same");
  // TCanvas *cv4 = new TCanvas();
  // cv4->cd();
  // tracker->fNEntries->Draw();
  for (int i = 0; i < 6; ++i)
  {
    TCanvas *cv = new TCanvas(Form("Doublets%i",i),Form("Doublets%i",i),1400,500);
    cv->Divide(2);
    cv->cd(1);
    tracker->fGDZ[i]->Draw();
    tracker->fFDZ[i]->SetLineColor(kRed);
    tracker->fFDZ[i]->Draw("same");
    cv->cd(2);
    tracker->fGDXY[i]->Draw();
    tracker->fFDXY[i]->SetLineColor(kRed);
    tracker->fFDXY[i]->Draw("same");
  }
  for (int i = 0; i < 5; ++i)
  {
    TCanvas *cv = new TCanvas(Form("Cells%i",i),Form("Cells%i",i),1400,500);
    cv->Divide(2);
    cv->cd(1);
    tracker->fGDCAZ[i]->Draw();
    tracker->fFDCAZ[i]->SetLineColor(kRed);
    tracker->fFDCAZ[i]->Draw("same");
    cv->cd(2);
    tracker->fGDCAXY[i]->Draw();
    tracker->fFDCAXY[i]->SetLineColor(kRed);
    tracker->fFDCAXY[i]->Draw("same");
  }
  // TCanvas *ccv4 = new TCanvas("iparams","Impact Parameters",1000,700);
  // ccv4->Divide(2);
  // ccv4->cd(1);
  // tracker->fIPGoodXY->Draw();
  // tracker->fIPFakeXY->SetLineColor(kRed);
  // tracker->fIPFakeXY->Draw("same");
  // ccv4->cd(2);
  // tracker->fIPGoodZ->Draw();
  // tracker->fIPFakeZ->SetLineColor(kRed);
  // tracker->fIPFakeZ->Draw("same");
  #endif
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
