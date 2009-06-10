#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TClassTable.h>
#include <TGraph.h>
#include <TTree.h>
#include <TH1.h>
#include <TFile.h>
#include <TH2.h>
#include <TInterpreter.h>
#include <TStyle.h>
#include "AliHeader.h"
#include "AliITSdEdxAnalyzer.h"
#include "AliRunLoader.h"
#include "AliESDEvent.h"
#endif

Bool_t MakeITSdEdxAnalysis(TString path=".",Int_t nmombins=16){

  if (gClassTable->GetID("AliRun") < 0) {
    gInterpreter->ExecuteMacro("loadlibs.C");
  }
 
  TString galicefile = path+"/galice.root";
  TString esdfile = path+"/AliESDs.root";

  AliRunLoader* runLoader = AliRunLoader::Open(galicefile.Data());
  if (!runLoader) {
    printf("Error in getting run loader");
    return kFALSE;
  }
  runLoader->LoadHeader();
  runLoader->LoadKinematics();

  TFile* esdFile = TFile::Open(esdfile.Data());
  if (!esdFile || !esdFile->IsOpen()) {
    printf("Error in opening ESD file");
    return kFALSE;
  }

  AliESDEvent * esd = new AliESDEvent;
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    printf("Error: no ESD tree found");
    return kFALSE;
  }
  esd->ReadFromTree(tree);
  AliITSdEdxAnalyzer* enan=new AliITSdEdxAnalyzer(nmombins,0.1,3.1);
  enan->SetMIPdEdx(79.);

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    runLoader->GetEvent(iEvent);
    tree->GetEvent(iEvent);
    if (!esd) {
      printf("Error: no ESD object found for event %d", iEvent);
      return kFALSE;
    }
    AliStack* stack = runLoader->Stack();
    enan->ReadEvent(esd,stack);
  }
  enan->WriteHistos();



  TCanvas* c1=new TCanvas("c1","dEdx - pions - layer 4");
  c1->Divide(4,4);
  for(Int_t i=0; i<nmombins;i++){
    TH1F* h=enan->GetSingleLayerdEdxHisto(4,211,i);
    c1->cd(i+1);
    h->Draw();
  }


  TCanvas* c1b=new TCanvas("c1b","dEdx - pions - Truncated mean");
  c1b->Divide(4,4);
  for(Int_t i=0; i<nmombins;i++){
    TH1F* h=enan->GetTruncatedMeandEdxHisto(211,i);
    c1b->cd(i+1);
    h->Draw();
  }

  TCanvas* c2=new TCanvas("c2","Delta dEdx - pions - layer 4");
  c2->Divide(4,4);
  for(Int_t i=0; i<nmombins;i++){
    TH1F* h=enan->GetSingleLayerDeltadEdxHisto(4,211,i);
    c2->cd(i+1);
    h->Draw();
  }

  TCanvas* c2b=new TCanvas("c2b","Delta dEdx - pions - Truncated Mean");
  c2b->Divide(4,4);
  for(Int_t i=0; i<nmombins;i++){
    TH1F* h=enan->GetTruncatedMeanDeltadEdxHisto(211,i);
    c2b->cd(i+1);
    h->Draw();
  }

  gStyle->SetPalette(1);

  TCanvas* c3=0;
  c3=new TCanvas("c3","dEdx vs. p - pions - layer 4");
  TH2F* h2=enan->GetSingleLayerdEdxVsPHisto(4,211);
  h2->Draw("colz");
  enan->SetUseBBFromAliExternalTrackParam();
  TGraph* gpion1=enan->GetBetheBlochGraph(211);
  gpion1->Draw("LSAME");
  enan->SetUseBBFromAliITSpidESD();
  TGraph* gpion2=enan->GetBetheBlochGraph(211);
  gpion2->SetLineColor(2);
  gpion2->Draw("LSAME");

  TCanvas* c3b=0;
  c3b=new TCanvas("c3b","dEdx vs. p - pions - Truncated Mean");
  TH2F* h2b=enan->GetTruncatedMeandEdxVsPHisto(211);
  h2b->Draw("colz");
  gpion1->Draw("LSAME");
  gpion2->Draw("LSAME");

  return kTRUE;
}
