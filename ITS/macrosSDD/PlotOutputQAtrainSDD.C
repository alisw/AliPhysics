#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include "AliITSgeomTGeo.h"
#endif

Double_t LangausFun(Double_t *x, Double_t *par);


void PlotOutputQAtrainSDD(TString option="local",
			  Int_t nRun=0,
			  TString period="LHC11a",
			  TString qaTrain="",
			  TString fileName="QAresults.root"){

  gStyle->SetOptStat(0);
  //  gStyle->SetOptTitle(0);
  TFile *f;
  TString path;
  Int_t year=2011;
  if(period.Contains("LHC10")) year=2010;
  else if(period.Contains("LHC09")) year=2009;

  if(option.Contains("local")){
    f=new TFile(fileName.Data());
    printf("Opened file %s\n",f->GetName());
  }else{
    TGrid::Connect("alien:");
    path=Form("/alice/data/%d/%s/%09d/ESDs/",year,period.Data(),nRun);
    printf("search in path %s\n",path.Data());
    if(!gGrid||!gGrid->IsConnected()) {
      printf("gGrid not found! exit macro\n");
      return;
    }
    TGridResult *gr = gGrid->Query(path,fileName);
    Int_t nFiles = gr->GetEntries();
    if (nFiles < 1){
      printf("QA file for run %d not found\n",nRun);
      return;
    }
    printf("================>%d files found\n", nFiles);
    if(qaTrain.Length()>0){
      Int_t found=kFALSE;
      for (Int_t iFil = 0; iFil <nFiles ; iFil++) { 
 	fileName=Form("%s",gr->GetKey(iFil,"turl"));
	TString isMerged=fileName;
	isMerged.Remove(isMerged.Sizeof()-16); 
	isMerged.Remove(0,isMerged.Sizeof()-5);
	if(!isMerged.Contains("QA")) continue;
	if(fileName.Contains(qaTrain.Data())){
	  found=kTRUE;
	  break;
	}
      }
      if(!found){
	printf(" File from %s train not found\n",qaTrain.Data());
	return;
      }
    }else{
      Int_t theFile=0;
      Int_t maxVer=0;
      if (nFiles > 1){
	for (Int_t iFil = 0; iFil <nFiles ; iFil++) { 
	  fileName=Form("%s",gr->GetKey(iFil,"turl"));
	  TString isMerged=fileName;
	  isMerged.Remove(isMerged.Sizeof()-16); 
	  isMerged.Remove(0,isMerged.Sizeof()-5);
	  if(!isMerged.Contains("QA")) continue;
	  TString cutFilename=fileName.ReplaceAll("/QAresults.root","");
	  Int_t last=cutFilename.Sizeof()-3;
	  cutFilename.Remove(0,last);
	  Int_t qaver=atoi(cutFilename.Data());
	  if(qaver>maxVer){
	    maxVer=qaver;
	    theFile=iFil;
	  }
	}
      }
      fileName=Form("%s",gr->GetKey(theFile,"turl"));
    }
    f=TFile::Open(fileName.Data());
  }

  TDirectoryFile* df=(TDirectoryFile*)f->Get("SDD_Performance");
  if(!df){
    printf("SDD_Performance MISSING -> Exit\n");
    return;
  }
  TList* l=(TList*)df->Get("coutputRP");
  if(!df){
    printf("coutputRP TList MISSING -> Exit\n");
    return;
  }

  TH1F* hcllay=(TH1F*)l->FindObject("hCluInLay");
  TH1F* hapmod=(TH1F*)l->FindObject("hAllPmod");
  TH1F* hgpmod=(TH1F*)l->FindObject("hGoodPmod");
  TH1F* hmpmod=(TH1F*)l->FindObject("hMissPmod");
  TH1F* hbrmod=(TH1F*)l->FindObject("hBadRegmod");
  TH1F* hskmod=(TH1F*)l->FindObject("hSkippedmod");
  TH1F* hoamod=(TH1F*)l->FindObject("hOutAccmod");
  TH1F* hnrmod=(TH1F*)l->FindObject("hNoRefitmod");

  //  TH1F* hapxl=(TH1F*)l->FindObject("hAllPxloc");
  TH1F* hgpxl=(TH1F*)l->FindObject("hGoodPxloc");
  TH1F* hmpxl=(TH1F*)l->FindObject("hMissPxloc");
  TH1F* hbrxl=(TH1F*)l->FindObject("hBadRegxloc");
  //  TH1F* hapzl=(TH1F*)l->FindObject("hAllPzloc");
  TH1F* hgpzl=(TH1F*)l->FindObject("hGoodPzloc");
  TH1F* hmpzl=(TH1F*)l->FindObject("hMissPzloc");
  TH1F* hbrzl=(TH1F*)l->FindObject("hBadRegzloc");

  TH2F* hClSizAn=(TH2F*)l->FindObject("hCluSizAn");
  TH2F* hClSizTb=(TH2F*)l->FindObject("hCluSizTb");

  TH2F* hdedx3=(TH2F*)l->FindObject("hdEdxL3VsP");
  TH2F* hdedx4=(TH2F*)l->FindObject("hdEdxL4VsP");
  TH2F* hdedxmod=(TH2F*)l->FindObject("hdEdxVsMod");
  

  TH1F* hmodR=(TH1F*)l->FindObject("hRPMod");
  TH1F* hmodT=(TH1F*)l->FindObject("hTPMod");
  TH1F* hgamod=(TH1F*)l->FindObject("hGAMod");

  TH2F* h2dmodR3=new TH2F("h2dmodR3","Rec Points, Layer 3",6,0.5,6.5,14,0.5,14.5);
  TH2F* h2dmodR4=new TH2F("h2dmodR4","Rec Points, Layer 4",8,0.5,8.5,22,0.5,22.5);
  TH2F* h2dmodT3=new TH2F("h2dmodT3","Track Points, Layer 3",6,0.5,6.5,14,0.5,14.5);
  TH2F* h2dmodT4=new TH2F("h2dmodT4","Track Points, Layer 4",8,0.5,8.5,22,0.5,22.5);
  TH2F* h2dmodR3N=new TH2F("h2dmodR3N","Rec Points/GoodAnode/Event, Layer 3",6,0.5,6.5,14,0.5,14.5);
  TH2F* h2dmodR4N=new TH2F("h2dmodR4N","Rec Points/GoodAnode/Event, Layer 4",8,0.5,8.5,22,0.5,22.5);
  TH2F* h2dmodT3N=new TH2F("h2dmodT3N","Track Points/GoodAnode/Event, Layer 3",6,0.5,6.5,14,0.5,14.5);
  TH2F* h2dmodT4N=new TH2F("h2dmodT4N","Track Points/GoodAnode/Event, Layer 4",8,0.5,8.5,22,0.5,22.5);
  TH1F* hmodRN=new TH1F("hmodRN","Normalized Rec Points per Module",260,239.5,499.5);
  TH1F* hmodTN=new TH1F("hmodTN","Normalized Track Points per Module",260,239.5,499.5);

  TH1F* hev=(TH1F*)l->FindObject("hNEvents");
  Int_t nTotEvents=hev->GetBinContent(2);
  Int_t nTrigEvents=hev->GetBinContent(3);
  Int_t nEvents=nTotEvents;
  printf("---- Statistics ----\n");
  printf("Number of Events = %d\n",nTotEvents);
  if(nTrigEvents>0){ 
    printf("Number of Triggered Events = %d\n",nTrigEvents);
    nEvents=nTrigEvents;
  }else{
    printf("No request on the trigger done when running the task\n");
  }
  Int_t bestMod=0;
  for(Int_t iMod=0; iMod<260;iMod++){
    Int_t gda=(Int_t)hgamod->GetBinContent(iMod+1);
    if(gda>bestMod) bestMod=gda;
  }
  Int_t nChunks=1;
  if(bestMod>512){
    nChunks=(Int_t)(bestMod/512.+0.5);
  }
  printf("Chunks merged = %d\n",nChunks);
  hgamod->Scale(1./nChunks);
  TCanvas* cgan=new TCanvas("cgan","Good Anodes");
  cgan->SetTickx();
  cgan->SetTicky();
  hgamod->SetMarkerStyle(20);
  hgamod->SetMarkerSize(0.6);
  hgamod->SetMinimum(-10.);
  hgamod->Draw("P");
  hgamod->GetXaxis()->SetTitle("SDD Module Id");
  hgamod->GetYaxis()->SetTitle("Number of good anodes");
  cgan->Update();

  printf("---- Modules with > 2%% of bad anodes ----\n");
  for(Int_t iMod=0; iMod<260; iMod++){
    Int_t idMod=iMod+240;
    Float_t rps=hmodR->GetBinContent(iMod+1);
    Float_t tps=hmodT->GetBinContent(iMod+1);
    Float_t ga=hgamod->GetBinContent(iMod+1);
    if(ga<500){
      printf("Module %d - Good Anodes = %d\n",idMod,(Int_t)ga);
    }
    Float_t rpsN=0.;
    Float_t tpsN=0.;
    Float_t erpsN=0.;
    Float_t etpsN=0.;
    if(ga>0 && nEvents>0){
      rpsN=rps/ga/(Float_t)nEvents;
      tpsN=tps/ga/(Float_t)nEvents;
      erpsN=TMath::Sqrt(rps)/ga/(Float_t)nEvents;
      etpsN=TMath::Sqrt(tps)/ga/(Float_t)nEvents;
    }
    hmodRN->SetBinContent(iMod+1,rpsN);
    hmodTN->SetBinContent(iMod+1,tpsN);
    hmodRN->SetBinError(iMod+1,erpsN);
    hmodTN->SetBinError(iMod+1,etpsN);
    Int_t iLay,iLad,iDet;
    AliITSgeomTGeo::GetModuleId(idMod,iLay,iLad,iDet);
    if(iLay==3){
      h2dmodR3->SetBinContent(iDet,iLad,rps);
      h2dmodT3->SetBinContent(iDet,iLad,tps);
      h2dmodR3N->SetBinContent(iDet,iLad,rpsN);
      h2dmodT3N->SetBinContent(iDet,iLad,tpsN);
    }
    else if(iLay==4){
      h2dmodR4->SetBinContent(iDet,iLad,rps);
      h2dmodT4->SetBinContent(iDet,iLad,tps);
      h2dmodR4N->SetBinContent(iDet,iLad,rpsN);
      h2dmodT4N->SetBinContent(iDet,iLad,tpsN);
    }
  }
  if(nEvents<1) return;

  gStyle->SetPalette(1);

  if(hmodR->GetEntries()>0){
    TCanvas* cmodR=new TCanvas("cmodR","RecPoint Occup",1200,1200);
    cmodR->Divide(2,3);
    cmodR->cd(1);
    gPad->SetLeftMargin(0.14);
    hmodR->Draw();
    hmodR->GetXaxis()->SetTitle("SDD Module Id");
    hmodR->GetYaxis()->SetTitle("RecPoints");
    hmodR->GetYaxis()->SetTitleOffset(1.55);
    cmodR->cd(2);
    gPad->SetLeftMargin(0.14);
    hmodRN->Draw("E");
    hmodRN->GetXaxis()->SetTitle("SDD Module Id");
    hmodRN->GetYaxis()->SetTitle("RecPoints/GoodAnode/Event");
    hmodRN->GetYaxis()->SetTitleOffset(1.55);
    cmodR->cd(3);
    gPad->SetLeftMargin(0.14);
    h2dmodR3->Draw("colz");
    h2dmodR3->GetXaxis()->SetTitle("Detector");
    h2dmodR3->GetYaxis()->SetTitle("Ladder");
    cmodR->cd(4);
    gPad->SetLeftMargin(0.14);
    h2dmodR3N->Draw("colz");
    h2dmodR3N->GetXaxis()->SetTitle("Detector");
    h2dmodR3N->GetYaxis()->SetTitle("Ladder");
    cmodR->cd(5);
    gPad->SetLeftMargin(0.14);
    h2dmodR4->Draw("colz");
    h2dmodR4->GetXaxis()->SetTitle("Detector");
    h2dmodR4->GetYaxis()->SetTitle("Ladder");
    cmodR->cd(6);
    gPad->SetLeftMargin(0.14);
    gPad->SetLeftMargin(0.14);
    h2dmodR4N->Draw("colz");
    h2dmodR4N->GetXaxis()->SetTitle("Detector");
    h2dmodR4N->GetYaxis()->SetTitle("Ladder");
    cmodR->Update();
  }


  if(hmodT->GetEntries()>0){
    TCanvas* cmodT=new TCanvas("cmodT","TrackPoint Occup",1200,1200);
    cmodT->Divide(2,3);
    cmodT->cd(1);
    hmodT->Draw();
    hmodT->GetXaxis()->SetTitle("SDD Module Id");
    hmodT->GetYaxis()->SetTitle("TrackPoints");
    hmodT->GetYaxis()->SetTitleOffset(1.4);
    cmodT->cd(2);
    gPad->SetLeftMargin(0.14);
    hmodTN->Draw("E");
    hmodTN->GetXaxis()->SetTitle("SDD Module Id");
    hmodTN->GetYaxis()->SetTitle("TrackPoints");
    hmodTN->GetYaxis()->SetTitleOffset(1.4);
    cmodT->cd(3);
    gPad->SetLeftMargin(0.14);
    h2dmodT3->Draw("colz");
    h2dmodT3->GetXaxis()->SetTitle("Detector");
    h2dmodT3->GetYaxis()->SetTitle("Ladder");
    cmodT->cd(4);
    gPad->SetLeftMargin(0.14);
    h2dmodT3N->Draw("colz");
    h2dmodT3N->GetXaxis()->SetTitle("Detector");
    h2dmodT3N->GetYaxis()->SetTitle("Ladder");  
    cmodT->cd(5);
    gPad->SetLeftMargin(0.14);
    h2dmodT4->Draw("colz");
    h2dmodT4->GetXaxis()->SetTitle("Detector");
    h2dmodT4->GetYaxis()->SetTitle("Ladder");
    cmodT->cd(6);
    gPad->SetLeftMargin(0.14);
    h2dmodT4N->Draw("colz");
    h2dmodT4N->GetXaxis()->SetTitle("Detector");
    h2dmodT4N->GetYaxis()->SetTitle("Ladder");
    cmodT->Update();
  }

  TH1F* htplad3=(TH1F*)l->FindObject("hTPLad3");
  TH1F* htplad4=(TH1F*)l->FindObject("hTPLad4");
  TH1F* hgalad3=(TH1F*)l->FindObject("hGALad3");
  TH1F* hgalad4=(TH1F*)l->FindObject("hGALad4");
  TH1F* hnormOcc3=new TH1F("hnormOcc3","",14,-0.5,13.5);
  TH1F* hnormOcc4=new TH1F("hnormOcc4","",22,-0.5,21.5);
  Bool_t tpok=kFALSE;
  for(Int_t ilad=0;ilad<14;ilad++){ 
    Float_t occ=0.;
    Float_t eocc=0.;
    Int_t gd3=hgalad3->GetBinContent(ilad+1);
    if(gd3>0){
      occ=(Float_t)htplad3->GetBinContent(ilad+1)/(Float_t)gd3/(Float_t)nEvents;
      eocc=TMath::Sqrt((Float_t)htplad3->GetBinContent(ilad+1))/(Float_t)gd3/(Float_t)nEvents;
    }
    hnormOcc3->SetBinContent(ilad+1,occ);
    hnormOcc3->SetBinError(ilad+1,eocc);
  }
  for(Int_t ilad=0;ilad<22;ilad++){ 
    Float_t occ=0.;
    Float_t eocc=0.;
    Int_t gd4=hgalad4->GetBinContent(ilad+1);
    if(gd4>0){
      occ=(Float_t)htplad4->GetBinContent(ilad+1)/(Float_t)gd4/(Float_t)nEvents;
      eocc=TMath::Sqrt((Float_t)htplad4->GetBinContent(ilad+1))/(Float_t)gd4/(Float_t)nEvents;
    }
    hnormOcc4->SetBinContent(ilad+1,occ);
    hnormOcc4->SetBinError(ilad+1,eocc);
  }


  if(tpok){
    TCanvas* cn0=new TCanvas("cn0","Normalized Ladder Occupancy",1400,600);
    cn0->Divide(2,1);
    cn0->cd(1);
    gPad->SetLeftMargin(0.14);
    hnormOcc3->Draw();
    hnormOcc3->GetXaxis()->SetTitle("Ladder number (layer 3)");
    hnormOcc3->GetYaxis()->SetTitle("TrackPoints/GoodAnodes/Events");
    hnormOcc3->GetYaxis()->SetTitleOffset(1.35);
    cn0->cd(2);
    gPad->SetLeftMargin(0.14);
    hnormOcc4->Draw();
    hnormOcc4->GetXaxis()->SetTitle("Ladder number (layer 4)");
    hnormOcc4->GetYaxis()->SetTitle("TrackPoints/GoodAnode/Events");
    hnormOcc4->GetYaxis()->SetTitleOffset(1.35);
    cn0->Update();
  }

  if(hcllay){
    Double_t norm=hcllay->GetBinContent(1);
    if(norm>0.){
      hcllay->Scale(1./norm);
      hcllay->SetTitle("");
      hcllay->GetXaxis()->SetRange(2,7);
      hcllay->SetMinimum(0.);
      hcllay->SetMaximum(1.1);
      hcllay->SetMarkerStyle(23);
      TCanvas* ceffL=new TCanvas("ceffL","PointPerLayer",1000,600);
      ceffL->SetGridy();
      hcllay->Draw();
      hcllay->GetXaxis()->SetTitle("Layer");
      hcllay->GetYaxis()->SetTitle("Fraction of tracks with point in layer");
      ceffL->Update();
    }
  }

  hgpmod->SetTitle("");
  TCanvas* ceff0=new TCanvas("ceff0","ModuleIndexInfo",1000,600);
  hgpmod->Draw();
  hgpmod->GetXaxis()->SetTitle("SDD Module Id");
  hgpmod->GetYaxis()->SetTitle("Number of tracks");
  hmpmod->SetLineColor(2);
  hmpmod->SetMarkerColor(2);
  hmpmod->SetMarkerStyle(22);
  hmpmod->SetMarkerSize(0.5);
  hmpmod->Draw("psame");
  hbrmod->SetLineColor(kGreen+1);
  hbrmod->SetMarkerColor(kGreen+1);
  hbrmod->SetMarkerStyle(20);
  hbrmod->SetMarkerSize(0.5);
  hbrmod->Draw("same");
  hskmod->SetLineColor(kYellow);
  hskmod->Draw("same");
  hoamod->SetLineColor(4);
  hoamod->Draw("same");
  hnrmod->SetLineColor(6);
  hnrmod->Draw("same");
  TLatex* t1=new TLatex(0.7,0.85,"Good Point");
  t1->SetNDC();
  t1->SetTextColor(1);
  t1->Draw();
  TLatex* t2=new TLatex(0.7,0.8,"Missing Point");
  t2->SetNDC();
  t2->SetTextColor(2);
  t2->Draw();
  TLatex* t3=new TLatex(0.7,0.75,"Bad Region");
  t3->SetNDC();
  t3->SetTextColor(kGreen+1);
  t3->Draw();
  ceff0->Update();

  TH1F* heff=new TH1F("heff","",260,239.5,499.5);
  TH1F* hfracskip=new TH1F("hfracskip","",260,239.5,499.5);

  for(Int_t imod=0; imod<260;imod++){
    Float_t numer=hgpmod->GetBinContent(imod+1)+hbrmod->GetBinContent(imod+1)+hoamod->GetBinContent(imod+1)+hnrmod->GetBinContent(imod+1);
    Float_t denom=hapmod->GetBinContent(imod+1)-hskmod->GetBinContent(imod+1);
    Float_t eff=0.;
    Float_t erreff=0.;
    if(denom>0){
      eff=numer/denom;
      erreff=TMath::Sqrt(eff*(1-eff)/denom);
    }
    heff->SetBinContent(imod+1,eff);
    heff->SetBinError(imod+1,erreff);
    Float_t numer2=hskmod->GetBinContent(imod+1);
    Float_t denom2=hapmod->GetBinContent(imod+1);
    Float_t frac=0.;
    Float_t efrac=0.;
    if(denom2>0.){
      frac=numer2/denom2;
      efrac=TMath::Sqrt(frac*(1-frac)/denom2);
    }
    hfracskip->SetBinContent(imod+1,frac);
    hfracskip->SetBinError(imod+1,efrac);

  }

  printf("---- Modules with efficiency < 90%% ----\n");
  TCanvas* ceff1=new TCanvas("ceff1","Efficiency",1000,1000);
  ceff1->Divide(1,2);
  ceff1->cd(1);
  heff->Draw();
  heff->GetXaxis()->SetTitle("SDD Module Id");
  heff->GetYaxis()->SetTitle("Fraction of tracks with point in good region");
  for(Int_t ibin=1; ibin<=heff->GetNbinsX(); ibin++){
    Float_t e=heff->GetBinContent(ibin);
    if(e<0.9){
      Int_t iMod=(Int_t)heff->GetBinCenter(ibin);
      Int_t lay,lad,det;
      AliITSgeomTGeo::GetModuleId(iMod,lay,lad,det);
      printf("Module %d - Layer %d Ladder %2d Det %d  -   Eff. %.3f\n",iMod,lay,lad,det,heff->GetBinContent(ibin));
    }
  }
  ceff1->cd(2);
  hfracskip->Draw();
  hfracskip->GetXaxis()->SetTitle("SDD Module Id");
  hfracskip->GetYaxis()->SetTitle("Fraction of tracks with skipped SDD");
  

  if(hgpxl){
    hgpxl->SetTitle("");
    hgpzl->SetTitle("");
    TCanvas* ceff2=new TCanvas("ceff2","LocalCoord",1000,600);
    ceff2->Divide(2,1);
    ceff2->cd(1);
    hgpxl->Draw();
    hgpxl->GetXaxis()->SetTitle("Xlocal (cm)");
    hgpxl->GetYaxis()->SetTitle("Number of tracks");
    hmpxl->SetLineColor(2);
    hmpxl->SetMarkerColor(2);
    hmpxl->SetMarkerStyle(22);
    hmpxl->SetMarkerSize(0.5);
    hmpxl->Draw("psame");
    hbrxl->SetLineColor(kGreen+1);
    hbrxl->SetMarkerColor(kGreen+1);
    hbrxl->SetMarkerStyle(20);
    hbrxl->SetMarkerSize(0.5);
    hbrxl->Draw("same");
    t1->Draw();
    t2->Draw();
    t3->Draw();
    ceff2->cd(2);
    hgpzl->Draw();
    hgpzl->GetXaxis()->SetTitle("Zlocal (cm)");
    hgpzl->GetYaxis()->SetTitle("Number of tracks");
    hmpzl->SetLineColor(2);
    hmpzl->SetMarkerColor(2);
    hmpzl->SetMarkerStyle(22);
    hmpzl->SetMarkerSize(0.5);
    hmpzl->Draw("psame");
    hbrzl->SetLineColor(kGreen+1);
    hbrzl->SetMarkerColor(kGreen+1);
    hbrzl->SetMarkerStyle(20);
    hbrzl->SetMarkerSize(0.5);
    hbrzl->Draw("same");
    t1->Draw();
    t2->Draw();
    t3->Draw();
    ceff2->Update();
  }


  if(hClSizAn && hClSizTb){
    TCanvas* ccs=new TCanvas("ccs","Cluster Size",1200,600);
    ccs->Divide(2,1);
    ccs->cd(1);
    gPad->SetLogz();
    gPad->SetRightMargin(0.12);
    hClSizAn->Draw("colz");
    hClSizAn->GetXaxis()->SetTitle("Drift Time (ns)");
    hClSizAn->GetYaxis()->SetTitle("Cluster Size - Anodes");
    ccs->cd(2);
    gPad->SetLogz();
    gPad->SetRightMargin(0.12);
    hClSizTb->Draw("colz");
    hClSizTb->GetXaxis()->SetTitle("Drift Time (ns)");
    hClSizTb->GetYaxis()->SetTitle("Cluster Size - Time Bins");
    ccs->Update();
  }

  TH1F* htimR=(TH1F*)l->FindObject("hDrTimRP");
  TH1F* htimT=(TH1F*)l->FindObject("hDrTimTPAll");
  TH1F* htimTe=(TH1F*)l->FindObject("hDrTimTPExtra");
  TH1F* htimTne=(TH1F*)l->FindObject("hDrTimTPNoExtra");
  htimR->Rebin(4);
  htimT->Rebin(4);
  htimTe->Rebin(4);
  htimTne->Rebin(4);
  htimR->SetLineWidth(2);
  htimT->SetLineWidth(2);
  htimTe->SetLineWidth(2);
  htimTne->SetLineWidth(2);

  TCanvas* ctim=new TCanvas("ctim","DriftTime",1400,600);
  ctim->Divide(2,1);
  ctim->cd(1);
  htimR->Draw(); 
  htimR->GetYaxis()->SetTitleOffset(1.2);
  htimR->GetXaxis()->SetTitle("Drift Time (ns)");
  htimR->GetYaxis()->SetTitle("RecPoints");
  ctim->cd(2);
  htimT->Draw();
  htimTe->SetLineColor(2);
  htimTe->Draw("same");
  htimTne->SetLineColor(4);
  htimTne->Draw("same");
  htimT->GetXaxis()->SetTitle("Drift Time (ns)");
  htimT->GetYaxis()->SetTitle("TrackPoints");
  htimT->GetYaxis()->SetTitleOffset(1.2);
  TLatex* ta=new TLatex(0.5,0.85,"All Clusters");
  ta->SetNDC();
  ta->SetTextColor(1);
  ta->Draw();
  TLatex* te=new TLatex(0.5,0.8,"Extra Clusters");
  te->SetNDC();
  te->SetTextColor(2);
  te->Draw();
  TLatex* tn=new TLatex(0.5,0.75,"Non-Extra Clusters");
  tn->SetNDC();
  tn->SetTextColor(4);
  tn->Draw();
  ctim->Update();

  TCanvas* cdedx=new TCanvas("cdedx","dedx",1400,600);
  cdedx->Divide(3,1);
  cdedx->cd(1);
  gPad->SetLogz();
  hdedx3->Draw("col");
  hdedx3->GetXaxis()->SetTitle("P (GeV/c)");
  hdedx3->GetYaxis()->SetTitle("dE/dx (keV/300 #mum) Layer 3");
  hdedx3->GetYaxis()->SetTitleOffset(1.25);
  cdedx->cd(2);
  gPad->SetLogz();
  hdedx4->Draw("col");
  hdedx4->GetXaxis()->SetTitle("P (GeV/c)");
  hdedx4->GetYaxis()->SetTitle("dE/dx (keV/300 #mum) Layer 4");
  hdedx4->GetYaxis()->SetTitleOffset(1.25);
  cdedx->cd(3);
  gPad->SetLogz();
  hdedxmod->Draw("col"); 
  hdedxmod->GetXaxis()->SetTitle("SDD Module Id");
  hdedxmod->GetYaxis()->SetTitle("dE/dx (keV/300 #mum)");
  hdedxmod->GetYaxis()->SetTitleOffset(1.25);
  cdedx->Update();

  printf("---- dE/dx vs.DriftTime ----\n");
  TCanvas* csig=new TCanvas("csig","dedx vs. DriftTime",1000,700);
  csig->Divide(4,2);
  TH1F* hSigTim[8];
  TGraphErrors* gmpv=new TGraphErrors(0);
  TGraphErrors* gsigg=new TGraphErrors(0);
  TGraphErrors* gsigl=new TGraphErrors(0);
  gmpv->SetTitle("");
  gsigg->SetTitle("");
  gsigl->SetTitle("");
  Int_t iPoint=0;
  TF1 *lfun = new TF1("LangausFun",LangausFun,50.,300.,4);
  for(Int_t it=0; it<8; it++){
    hSigTim[it]=(TH1F*)l->FindObject(Form("hSigTimeInt%d",it));
    csig->cd(it+1);
    hSigTim[it]->Draw();
    if(hSigTim[it]->GetEntries()>200){
      lfun->SetLineWidth(2);
      lfun->SetParameter(0,5.);
      lfun->SetParameter(1,80.);
      lfun->SetParameter(2,hSigTim[it]->GetEntries()/10.);
      lfun->SetParameter(3,10.);
      lfun->SetParLimits(3,0.,20);

      hSigTim[it]->Fit("LangausFun","QLR");
      hSigTim[it]->GetXaxis()->SetTitle(Form("dE/dx, time interval %d",it+1));
      hSigTim[it]->GetYaxis()->SetTitle("Events");
      Float_t mpv=lfun->GetParameter(1);
      Float_t empv=lfun->GetParError(1);
      Float_t sig=lfun->GetParameter(3);
      Float_t esig=lfun->GetParError(3);
      Float_t sigl=lfun->GetParameter(0);
      Float_t esigl=lfun->GetParError(0);
      gmpv->SetPoint(iPoint,(Float_t)it,mpv);
      gmpv->SetPointError(iPoint,0.,empv);
      gsigg->SetPoint(iPoint,(Float_t)it,sig);
      gsigg->SetPointError(iPoint,0.,esig);
      gsigl->SetPoint(iPoint,(Float_t)it,sigl);
      gsigl->SetPointError(iPoint,0.,esigl);
      ++iPoint;
      gPad->Update();
      printf("Bin %d - MPV=%.3f  \t SigmaLandau=%.3f  \t SigmaGaus=%.3f\n",it,mpv,sigl,sig);
    }
  }

  TCanvas* cpars=new TCanvas("cpars","Params",800,900);
  cpars->Divide(1,3,0.01,0.);
  cpars->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetFrameLineWidth(2);
  gPad->SetTickx();
  gPad->SetTicky();
  gmpv->SetMarkerStyle(20);
  //  gmpv->SetMinimum(0);
  //  gmpv->SetMaximum(120);
  gmpv->GetXaxis()->SetLimits(-0.2,6.8);
  gmpv->Draw("AP");
  //  gmpv->GetXaxis()->SetTitle("Drift Time interval number");
  gmpv->GetYaxis()->SetTitle("Landau MPV (keV)");
  gmpv->GetXaxis()->SetTitleSize(0.05);
  gmpv->GetYaxis()->SetTitleSize(0.05);
  gmpv->GetYaxis()->SetTitleOffset(1.2);
  cpars->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetFrameLineWidth(2);
  gPad->SetTickx();
  gPad->SetTicky();
  gsigl->SetMarkerStyle(20);
  gsigl->GetXaxis()->SetLimits(-0.2,6.8);
  gsigl->Draw("AP");
  //  gsigl->GetXaxis()->SetTitle("Drift Time interval number");
  gsigl->GetYaxis()->SetTitle("#sigma_{Landau} (keV)");
  gsigl->GetXaxis()->SetTitleSize(0.05);
  gsigl->GetYaxis()->SetTitleSize(0.05);
  gsigl->GetYaxis()->SetTitleOffset(1.2);
  cpars->cd(3);
  gPad->SetLeftMargin(0.14);
  gPad->SetFrameLineWidth(2);
  gPad->SetTickx();
  gPad->SetTicky();
  gsigg->SetMarkerStyle(20);
  gsigg->GetXaxis()->SetLimits(-0.2,6.8);
  gsigg->Draw("AP");
  gsigg->GetXaxis()->SetTitle("Drift Time interval number");
  gsigg->GetYaxis()->SetTitle("#sigma_{Gauss} (keV)");
  gsigg->GetXaxis()->SetTitleSize(0.05);
  gsigg->GetYaxis()->SetTitleSize(0.05);
  gsigg->GetYaxis()->SetTitleOffset(1.2);
  cpars->Update();

}


Double_t LangausFun(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

