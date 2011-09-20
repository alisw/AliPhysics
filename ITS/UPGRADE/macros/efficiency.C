
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TFile.h>
#include <Riostream.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TArrayF.h>
#include <TStyle.h>
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliTrackReference.h"
#include "AliITSsegmentationUpgrade.h"
#endif

Int_t minNpoints=3;
Bool_t IsTrackable(TClonesArray *trackRef);
TArrayF radii;

void efficiency(char *title="",Bool_t isPrimary=kTRUE){

  gSystem->Load("libITSUpgradeBase.so");
  gSystem->Load("libITSUpgradeRec.so");
  gSystem->Load("libITSUpgradeSim.so");

  // setting up the configuration to be considered
  AliITSsegmentationUpgrade *seg = new AliITSsegmentationUpgrade();
  if(!seg){
    printf("no segmentation info available... Exiting");
    return;
  }
  radii.Set(seg->GetNLayers());
  for(Int_t l=0; l<seg->GetNLayers(); l++){
    radii.AddAt(seg->GetRadius(l),l);
  }
  delete seg;

  TString titlet(title);
  titlet.ReplaceAll(" ","");

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  //booking plots
  // pt
  Int_t nbins= 64;
  Double_t ptmax = 1.6;

  TH1F *hPtRec = new TH1F("hPtRec"," p_{T} distribution of reconstructed tracks",nbins,0,ptmax);
  hPtRec->SetXTitle("p_{T}");
  hPtRec->SetYTitle(Form("entries / %3.0f MeV/c",1000*(ptmax/nbins)));

  TH1F *hPtRecMc = new TH1F("hPtRecMc"," MC p_{T} distribution of reconstructed tracks",nbins,0,ptmax);
  hPtRecMc->SetXTitle("p_{T} GeV/c");
  hPtRecMc->SetYTitle(Form("entries / %3.0f MeV/c",1000*(ptmax/nbins)));
 
  TH1F *hPtRecGood = new TH1F("hPtRecGood"," P_{t} distribution of reconstructed tracks [positive labels]",nbins,0,ptmax);
  hPtRecGood->SetXTitle("p_{T} GeV/c");
  hPtRecGood->SetYTitle(Form("entries / %3.0f MeV/c",1000*(ptmax/nbins)));

  TH1F *hPtMc  = new TH1F("hPtMc", " p_{T} distribution of MC tracks",nbins,0,ptmax);
  hPtMc->SetXTitle("p_{T} GeV/c");
  hPtMc->SetYTitle(Form("entries / %3.0f MeV/c",1000*(ptmax/nbins)));

  // Eta 
  Int_t nbinsEta = 24;
  Double_t eta[2]={-1.2,1.2}; 

  TH1F *hEtaRec = new TH1F("hEtaRec"," #eta distribution of reconstructed tracks",nbinsEta,eta[0],eta[1]);
  hEtaRec->SetXTitle("#eta");
  hEtaRec->SetYTitle(Form("entries / %3.3f ",((eta[1]-eta[0])/nbinsEta)));

  TH1F *hEtaRecMc = new TH1F("hEtaRecMc"," #eta Mc distribution of reconstructed tracks",nbinsEta,eta[0],eta[1]);
  hEtaRecMc->SetXTitle("#eta");
  hEtaRecMc->SetYTitle(Form("entries / %3.3f ",((eta[1]-eta[0])/nbinsEta)));
 
  TH1F *hEtaRecGood = new TH1F("hEtaRecGood"," #eta distribution of Good reconstructed tracks",nbinsEta,eta[0],eta[1]);
  hEtaRecGood->SetXTitle("#eta");
  hEtaRecGood->SetYTitle(Form("entries / %3.3f ",((eta[1]-eta[0])/nbinsEta)));

  TH1F *hEtaMc  = new TH1F("hEtaMc", " #eta distribution of MC tracks",nbinsEta,eta[0],eta[1]);
  hEtaMc->SetXTitle("#eta");
  hEtaMc->SetYTitle(Form("entries / %3.3f ",((eta[1]-eta[0])/nbinsEta)));


  Double_t ptlimit=0.0;

  Double_t cont=0; 
  Double_t ntrack=0;

  // Init simulation 
  gAlice=NULL;
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  if (!runLoader) {
    Error("OpengAlice", "no files, please check the path");
    return ;
  }
  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadTrackRefs();

  //Trackref
  TTree *trackRefTree = 0x0; 
  TClonesArray *trackRef = new TClonesArray("AliTrackReference",1000);
  
  // ESD
  TFile* esdFile = TFile::Open("AliESDs.root");
  if (!esdFile || !esdFile->IsOpen()) {
    Error("CheckESD", "opening ESD file %s failed", "AliESDs.root");
    return ;
  }
  AliESDEvent * esd = new AliESDEvent;
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("CheckESD", "no ESD tree found");
    return ;
  }
  esd->ReadFromTree(tree);

  printf("setting the low pt value to %f \n",ptlimit);
  // Loop over events

  for(Int_t i=0; i<runLoader->GetNumberOfEvents(); i++){
    //for(Int_t i=0; i< 1; i++){
    cout << "Event "<< i << endl;
    runLoader->GetEvent(i); 
 
    AliStack *stack = runLoader->Stack(); 
    trackRefTree=runLoader->TreeTR();
    TBranch *br = trackRefTree->GetBranch("TrackReferences");
    if(!br) {
      printf("no TR branch available , exiting \n");
      return;
    }
    br->SetAddress(&trackRef);

    Int_t nParticles = stack->GetNtrack();
    trackRefTree=runLoader->TreeTR();
    trackRefTree->SetBranchAddress("TrackReferences",&trackRef);
    for(Int_t p=0; p<nParticles; p++){ 
      if(stack->Particle(p)->Pt()<ptlimit) continue;
      if(isPrimary && !stack->IsPhysicalPrimary(p)) continue;
      trackRefTree->GetEntry(stack->TreeKEntry(p));

      if(IsTrackable(trackRef)) {
	cont++;
        hPtMc->Fill((stack->Particle(p))->Pt());
        hEtaMc->Fill((stack->Particle(p))->Eta());
      }
    }//loop on kine particle
    cout << " # trackable "<< cont; 

    // ESD
    tree->GetEvent(i);
    for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {
      AliESDtrack* track = esd->GetTrack(iTrack);
      if(track->Pt()<ptlimit) continue;      
      if(isPrimary && !stack->IsPhysicalPrimary(TMath::Abs(track->GetLabel()))) continue;
      if(track->GetNcls(0)<minNpoints) continue; // selection compatible with IsSelected
      ntrack++;
      hPtRec->Fill(track->Pt());
      hEtaRec->Fill(track->Eta());
      if(track->GetLabel()>=0) {
	hPtRecGood->Fill(stack->Particle(track->GetLabel())->Pt());
	hEtaRecGood->Fill(stack->Particle(track->GetLabel())->Eta());
      }
      hPtRecMc->Fill(stack->Particle(TMath::Abs(track->GetLabel()))->Pt());
      hEtaRecMc->Fill(stack->Particle(TMath::Abs(track->GetLabel()))->Eta());
    }// loop on tracks 
    cout << " tracks " << ntrack << endl;
   
  }//loop events
  
  cout << "Global Efficiency  "<< ntrack/cont << endl;

  TCanvas *c = new TCanvas();
  c->cd();
  c->cd()->SetLogy();
  c->cd()->SetGridy();
  c->cd()->SetGridx();
  c->cd()->SetObjectStat(0);
  hPtMc->SetMarkerStyle(20);
  hPtMc->Draw("err");
  hPtRec->SetLineColor(kRed);
  hPtRec->SetMarkerStyle(20);
  hPtRec->SetMarkerColor(kRed);
  hPtRec->Sumw2();
  hPtRec->Draw("sameerr");
  hPtRecMc->SetLineColor(kBlue);
  hPtRecMc->SetMarkerStyle(20);
  hPtRecMc->SetMarkerColor(kBlue);
  hPtRecMc->Sumw2();
  hPtRecMc->Draw("sameerr");
  hPtRecGood->SetLineColor(kGreen);
  hPtRecGood->SetMarkerStyle(20);
  hPtRecGood->SetMarkerColor(kGreen);
  hPtRecGood->Sumw2();
  hPtRecGood->Draw("sameerr");

  TLegend *leg = new TLegend(0.5,0.7,0.99,0.95);
  leg->SetHeader(Form("%s",title));
  leg->AddEntry(hPtRecGood,"Good tracks (positive labels)","p");
  leg->AddEntry(hPtRec,"Reconstructed tracks (Rec p_{T})","p");
  leg->AddEntry(hPtRecMc,"Reconstructed tracks ","p");
  leg->AddEntry(hPtMc,"MC trackable particles","p");
  leg->Draw();
  c->SaveAs(Form("ptdist_%s.png",titlet.Data()));

  TCanvas *cc = new TCanvas();
  cc->cd();
  cc->cd()->SetGridy();
  cc->cd()->SetGridx();
 
  Double_t max = 1.2;
 
  TH1F *effrec=new TH1F("hEffRec","tracking efficiency of reconstructed tracks",nbins,0,ptmax);
  effrec->Divide(hPtRecMc,hPtMc,1,1,"b");
  effrec->SetXTitle("p_{T} GeV/c");
  effrec->SetLineColor(kBlue);
  effrec->SetMarkerStyle(20);
  effrec->SetMarkerColor(kBlue);
  effrec->SetMaximum(max);
  effrec->SetMinimum(0);
  effrec->DrawCopy("err");


  TH1F *effrecgood=new TH1F("hEffRecGood","tracking efficiency of reconstructed tracks [positive labels]",nbins,0,ptmax);
  effrecgood->SetXTitle("P_{T}");
  effrecgood->Divide(hPtRecGood,hPtMc,1,1,"b");
  effrecgood->SetLineColor(kGreen);
  effrecgood->SetMarkerStyle(20);
  effrecgood->SetMarkerColor(kGreen);
  effrecgood->SetMaximum(max);
  effrecgood->SetMinimum(0);
  effrecgood->DrawCopy("sameerr");

  TH1F *purity=new TH1F("purity","track purity [positive labels / all tracks]",nbins,0,ptmax);
  purity->Divide(hPtRecGood,hPtRecMc,1,1,"b");
  purity->SetLineColor(kGray);
  purity->SetMarkerStyle(28);
  purity->SetMarkerColor(kOrange);
  purity->SetMaximum(max);
  purity->SetMinimum(0);
  purity->DrawCopy("sameerr");

  TH1F *fakes=new TH1F("fakes","fake ratio [negative labels / reconstructable]",nbins,0,ptmax);
  
  fakes->Add(hPtRecMc,hPtRecGood,1,-1);
  fakes->Divide(hPtMc);
  fakes->SetLineColor(kViolet);
  fakes->SetLineStyle(3);
  fakes->SetLineWidth(3);
  fakes->SetMarkerColor(kViolet);
  fakes->SetMaximum(max);
  fakes->SetMinimum(0);
  fakes->DrawCopy("samehist");

  TLegend *legEff = new TLegend(0.1,0.83,0.65,0.99);
  legEff->SetHeader(Form("%s",title));
  legEff->AddEntry(effrec,"Overall Efficiency ","p");
  legEff->AddEntry(effrecgood,"Efficiency good tracks (positive labels)","p");
  legEff->AddEntry(purity,"Purity (positive labels / All tracks)","p");
  legEff->AddEntry(fakes,"Fakes","l");
  legEff->Draw();

  cc->SaveAs(Form("efficiency_%s.png",titlet.Data()));


  TCanvas *cEta = new TCanvas();
  cEta->cd();
  cEta->cd()->SetGridy();
  cEta->cd()->SetGridx();
  cEta->cd()->SetObjectStat(0);
  hEtaMc->SetMarkerStyle(20);
  hEtaMc->Draw("err");
  hEtaMc->SetMinimum(0);
  hEtaRec->SetLineColor(kRed);
  hEtaRec->SetMarkerStyle(20);
  hEtaRec->SetMarkerColor(kRed);
  hEtaRec->Sumw2();
  hEtaRec->Draw("sameerr");
  hEtaRecMc->SetLineColor(kBlue);
  hEtaRecMc->SetMarkerStyle(20);
  hEtaRecMc->SetMarkerColor(kBlue);
  hEtaRecMc->Sumw2();
  hEtaRecMc->Draw("sameerr");
  hEtaRecGood->SetLineColor(kGreen);
  hEtaRecGood->SetMarkerStyle(20);
  hEtaRecGood->SetMarkerColor(kGreen);
  hEtaRecGood->Sumw2();
  hEtaRecGood->Draw("sameerr");

  TLegend *legEta = new TLegend(0.2,0.3,0.59,0.55);
  legEta->SetHeader(Form("%s",title));
  legEta->AddEntry(hEtaRecGood,"Good tracks (positive labels)","p");
  legEta->AddEntry(hEtaRec,"Reconstructed tracks (Rec p_{T})","p");
  legEta->AddEntry(hEtaRecMc,"Reconstructed tracks ","p");
  legEta->AddEntry(hEtaMc,"MC trackable particles","p");
  legEta->Draw();
  cEta->SaveAs(Form("etadist_%s.png",title));

  TCanvas *ccEta = new TCanvas();
  ccEta->cd();
  ccEta->cd()->SetGridy();
  ccEta->cd()->SetGridx();

  TH1F *effrecEta=new TH1F("hEffRecEta","#eta tracking efficiency of reconstructed tracks",nbinsEta,eta[0],eta[1]);
  effrecEta->Divide(hEtaRecMc,hEtaMc,1,1,"b");
  effrecEta->SetXTitle("#eta");
  effrecEta->SetLineColor(kBlue);
  effrecEta->SetMarkerStyle(20);
  effrecEta->SetMarkerColor(kBlue);
  effrecEta->SetMaximum(max);
  effrecEta->SetMinimum(0);
  effrecEta->DrawCopy("err");


  TH1F *effrecgoodEta=new TH1F("hEffRecGoodEta","#eta tracking efficiency of reconstructed tracks [positive labels]",nbinsEta,eta[0],eta[1]);
  effrecgoodEta->SetXTitle("eta");
  effrecgoodEta->Divide(hEtaRecGood,hEtaMc,1,1,"b");
  effrecgoodEta->SetLineColor(kGreen);
  effrecgoodEta->SetMarkerStyle(20);
  effrecgoodEta->SetMarkerColor(kGreen);
  effrecgoodEta->SetMaximum(max);
  effrecgoodEta->SetMinimum(0);
  effrecgoodEta->DrawCopy("sameerr");

  TH1F *purityEta=new TH1F("purityEta","#eta track purity [positive labels / all tracks]",nbinsEta,eta[0],eta[1]);
  purityEta->Divide(hEtaRecGood,hEtaRecMc,1,1,"b");
  purityEta->SetLineColor(kGray);
  purityEta->SetMarkerStyle(28);
  purityEta->SetMarkerColor(kOrange);
  purityEta->SetMaximum(max);
  purityEta->SetMinimum(0);
  purityEta->DrawCopy("sameerr");

  TH1F *fakesEta=new TH1F("fakesEta","#eta fake ratio [negative labels / reconstructable]",nbinsEta,eta[0],eta[1]);
  fakesEta->Add(hEtaRecMc,hEtaRecGood,1,-1);
  fakesEta->Divide(hEtaMc);
  fakesEta->SetLineColor(kViolet);
  fakesEta->SetLineStyle(3);
  fakesEta->SetLineWidth(3);
  fakesEta->SetMarkerColor(kViolet);
  fakesEta->SetMaximum(max);
  fakesEta->SetMinimum(0);
  fakesEta->DrawCopy("samehist");

  TLegend *legEffEta = new TLegend(0.3,0.3,0.75,0.55);
  legEffEta->SetHeader(Form("%s",title));
  legEffEta->AddEntry(effrecEta,"Overall Efficiency ","p");
  legEffEta->AddEntry(effrecgoodEta,"Efficiency good tracks (positive labels)","p");
  legEffEta->AddEntry(purityEta,"Purity (positive labels / All tracks)","p");
  legEffEta->AddEntry(fakesEta,"Fakes","l");
  legEffEta->Draw();

  ccEta->SaveAs(Form("efficiencyEta_%s.png",titlet.Data()));


}// main 

//___________________________________________________
Bool_t IsTrackable(TClonesArray *trackRef){

  Bool_t isOk=kFALSE;

  Int_t nTrackRef =0;
  TArrayF isInLayer(radii.GetSize());
  for(Int_t l=0; l<radii.GetSize(); l++ ) isInLayer.AddAt(0,l);
  for(Int_t nT =0; nT<trackRef->GetEntries(); nT++){
    AliTrackReference *trR = (AliTrackReference*)trackRef->At(nT);
    if(!trR) continue;
    if(trR->DetectorId()!=0)continue;
    Double_t rPart = TMath::Sqrt(trR->X()*trR->X()+trR->Y()*trR->Y());
    for(Int_t iLay=0;iLay<radii.GetSize();iLay++){
      if(TMath::Abs(rPart-radii.At(iLay))<0.01) isInLayer.AddAt(1,iLay);
    }
  }

  for(Int_t iLayer=0;iLayer<radii.GetSize();iLayer++){
    if(isInLayer.At(iLayer)) nTrackRef++;
  }

  if(nTrackRef>=minNpoints) isOk=kTRUE;
  return isOk; 
}

