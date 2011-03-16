#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include "AliITSgeomTGeo.h"
#include "AliESDEvent.h"
#endif

enum {kAll, kTPCITS, kITSsa, kITSpureSA};

void CheckSDDInESD(TString filename="AliESDs.root", Int_t optTracks=kAll){


  TFile* esdFile = TFile::Open(filename.Data());
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
  TH1F* hpt=new TH1F("hpt","",100,0.,10.);
  TH1F* hphi=new TH1F("hphi","",100,-1,1);
  TH1F* hlam=new TH1F("hlam","",100,-2.,2.);
  TH1F* halpha=new TH1F("halpha","",100,-7,7);
  TH1F* hitscl=new TH1F("hitscl","",7,-0.5,6.5);
  TH1F* htpccl=new TH1F("htpccl","",200,-0.5,199.5);
  TH1F* hitsmap=new TH1F("hitsmap","",64,-0.5,63.5);
  TH1F* hclulay=new TH1F("hclulay","",7,-1.5,5.5);

  TH1F* hvx=new TH1F("hvx","",100,-1.,1.);
  TH1F* hvy=new TH1F("hvy","",100,-1.,1.);
  TH1F* hvz=new TH1F("hvz","",100,-20.,20.);
  TH1F* hdedx3=new TH1F("hdedx3","",100,0.,300.);
  TH1F* hdedx4=new TH1F("hdedx4","",100,0.,300.);
  TH1F* hdedx5=new TH1F("hdedx5","",100,0.,300.);
  TH1F* hdedx6=new TH1F("hdedx6","",100,0.,300.);
  TH1F* hStatus=new TH1F("hStatus","",11,-1.5,9.5);


  // -- Local coordinates


  // -- Module histos

  TH1F* hAllPMod  = new TH1F("hAllPmod","Crossing Tracks vs. Module",260,239.5,499.5);
  TH1F* hGoodPMod  = new TH1F("hGoodPmod","PointsAssocToTrack per Module",260,239.5,499.5);
  TH1F* hBadRegMod  = new TH1F("hBadRegmod","Tracks in BadRegion per Module",260,239.5,499.5);
  TH1F* hMissPMod  = new TH1F("hMissPmod","Missing Points per Module",260,239.5,499.5);
  TH1F* hSkippedMod  = new TH1F("hSkippedmod","Tracks in Skipped Module",260,239.5,499.5);
  TH1F* hOutAccMod  = new TH1F("hOutAccmod","Tracks outside zAcc per Module",260,239.5,499.5);
  TH1F* hNoRefitMod  = new TH1F("hNoRefitmod","Points rejected in refit per Module",260,239.5,499.5);

  TH1F* hAllPXloc  = new TH1F("hAllPxloc","Crossing Tracks vs. Xloc",75, -3.75, 3.75);
  TH1F* hGoodPXloc  = new TH1F("hGoodPxloc","PointsAssocToTrack vs. Xloc",75, -3.75, 3.75);
  TH1F* hBadRegXloc  = new TH1F("hBadRegxloc","Tracks in BadRegion vs. Xloc",75, -3.75, 3.75);
  TH1F* hMissPXloc  = new TH1F("hMissPxloc","Missing Points vs. Xloc",75, -3.75, 3.75);
  TH1F* hAllPZloc  = new TH1F("hAllPzloc","Crossing Tracks vs. Zloc",77, -3.85, 3.85);
  TH1F* hGoodPZloc  = new TH1F("hGoodPzloc","PointsAssocToTrack vs. Zloc",77, -3.85, 3.85);
  TH1F* hBadRegZloc  = new TH1F("hBadRegzloc","Tracks in BadRegion vs. Zloc",77, -3.85, 3.85);
  TH1F* hMissPZloc  = new TH1F("hMissPzloc","Missing Points vs. Zloc",77, -3.85, 3.85);
  TH2F* hdEdxVsMod=new TH2F("hdEdxVsMod","dE/dx vs. mod",260,239.5,499.5,100,0.,500.);

  gStyle->SetPalette(1);
  

  for (Int_t iEvent = 0; iEvent < tree->GetEntries(); iEvent++) {
    tree->GetEvent(iEvent);
    if (!esd) {
      printf("Error: no ESD object found for event %d", iEvent);
      return;
    }
    cout<<"-------- Event "<<iEvent<<endl;
    printf(" Tracks # = %d\n",esd->GetNumberOfTracks());
    const AliESDVertex *spdv=esd->GetVertex();
    printf(" SPD Primary Vertex in %f %f %f with %d contributors\n",spdv->GetXv(),spdv->GetYv(),spdv->GetZv(),spdv->GetNContributors());
    const AliESDVertex *trkv=esd->GetPrimaryVertex();
    printf(" Track Primary Vertex with %d contributors\n",trkv->GetNContributors());
    if(spdv->IsFromVertexer3D()){
      hvx->Fill(spdv->GetXv());
      hvy->Fill(spdv->GetYv());
      hvz->Fill(spdv->GetZv());
    }
    Double_t itss[4];
    for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {
      AliESDtrack* track = esd->GetTrack(iTrack);
      Int_t nITSclus=track->GetNcls(0);
      UChar_t clumap=track->GetITSClusterMap();
      Int_t nPointsForPid=0;
      for(Int_t i=2; i<6; i++){
	if(clumap&(1<<i)) ++nPointsForPid;
      }
      //      track->PropagateTo(4.,5.);
      htpccl->Fill(track->GetNcls(1)<70);
      Int_t status=track->GetStatus();
      Bool_t tpcin=0;
      hStatus->Fill(-1.);
      if(status & AliESDtrack::kTPCin){
	tpcin=1;
	hStatus->Fill(0.);
      }
      if(status & AliESDtrack::kTPCout){
	hStatus->Fill(1.);
      }
      if(status & AliESDtrack::kTPCrefit){
	hStatus->Fill(2.);
      }
      Bool_t itsin=0;
      if(status & AliESDtrack::kITSin){
	itsin=1;
	hStatus->Fill(3.);
      }
      if(status & AliESDtrack::kITSout){
	hStatus->Fill(4.);
      }
      if(status & AliESDtrack::kITSrefit){
	hStatus->Fill(5.);
      }
      if(!tpcin && itsin){
	hStatus->Fill(6.);
      }
      if(status & AliESDtrack::kITSpureSA){
	hStatus->Fill(7.);
      }

      if(status & AliESDtrack::kITSrefit){
	if((optTracks==kTPCITS) && !(status & AliESDtrack::kTPCin)) continue;
	if((optTracks==kITSsa) && (status & AliESDtrack::kTPCin)) continue;
	if((optTracks==kITSsa) && (status & AliESDtrack::kITSpureSA)) continue;
	if((optTracks==kITSpureSA) && (status & AliESDtrack::kITSpureSA)) continue;

	 track->GetITSdEdxSamples(itss);
	//	printf("Track %d (label %d) in ITS with %d clusters clumap %d pointspid= %d\n",iTrack,track->GetLabel(),nITSclus,clumap,nPointsForPid);
	//printf("   dedx=%f %f %f %f\n",itss[0],itss[1],itss[2],itss[3]);
	hitscl->Fill(nITSclus);
	hdedx3->Fill(itss[0]);
	hdedx4->Fill(itss[1]);
	hdedx5->Fill(itss[2]);
	hdedx6->Fill(itss[3]);
	hitsmap->Fill(clumap);
	hclulay->Fill(-1.);
	for(Int_t iLay=0;iLay<6;iLay++){
	  if(clumap&1<<iLay) hclulay->Fill(iLay);
	}
	hpt->Fill(track->Pt());
	hphi->Fill(TMath::ASin(track->GetSnp()));
	hlam->Fill(TMath::ATan(track->GetTgl()));
	halpha->Fill(track->GetAlpha());
	Int_t iMod,status;
	Float_t xloc,zloc;
	for(Int_t iLay=2; iLay<=3; iLay++){
	  Bool_t ok=track->GetITSModuleIndexInfo(iLay,iMod,status,xloc,zloc);
	  if(ok){
	    iMod+=240;
	    hAllPMod->Fill(iMod);
	    hAllPXloc->Fill(xloc);
	    hAllPZloc->Fill(zloc);
	    if(status==1){
	      hGoodPMod->Fill(iMod);
	      hGoodPXloc->Fill(xloc);
	      hGoodPZloc->Fill(zloc);
	      if(track->Pt()>1.) hdEdxVsMod->Fill(iMod,itss[iLay-2]);
	    }
	    else if(status==2){ 
	      hBadRegMod->Fill(iMod);
	      hBadRegXloc->Fill(xloc);
	      hBadRegZloc->Fill(zloc);
	    }
	    else if(status==3) hSkippedMod->Fill(iMod);
	    else if(status==4) hOutAccMod->Fill(iMod);
	    else if(status==5){
	      hMissPMod->Fill(iMod);
	      hMissPXloc->Fill(xloc);
	      hMissPZloc->Fill(zloc);
	    }
	    else if(status==6) hNoRefitMod->Fill(iMod);
	  }
	}
      }
    }
  }
  Float_t norm=hclulay->GetBinContent(1);
  if(norm<1.) norm=1.;
  hclulay->Scale(1./norm);
  gStyle->SetLineWidth(2);

  TCanvas* c1=new TCanvas("c1","Track quantities",900,900);
  c1->Divide(2,2);
  c1->cd(1);
  htpccl->Draw();
  htpccl->GetXaxis()->SetTitle("Clusters in TPC ");
  c1->cd(2);
  hitscl->Draw();
  hitscl->GetXaxis()->SetTitle("Clusters in ITS ");
  c1->cd(3);
  hclulay->Draw();
  hclulay->GetXaxis()->SetRange(2,7);
  hclulay->GetXaxis()->SetTitle("# ITS Layer");
  hclulay->GetYaxis()->SetTitle("Fraction of tracks with point in Layer x");
  c1->cd(4);

  TCanvas* c2=new TCanvas("c2","dedx per Layer",900,900);
  c2->Divide(2,2);
  c2->cd(1);
  hdedx3->Draw();
  hdedx3->GetXaxis()->SetTitle("dE/dx Lay3");
  c2->cd(2);
  hdedx4->Draw();
  hdedx4->GetXaxis()->SetTitle("dE/dx Lay4");
  c2->cd(3);
  hdedx5->Draw();
  hdedx5->GetXaxis()->SetTitle("dE/dx Lay5");
  c2->cd(4);
  hdedx6->Draw();
  hdedx6->GetXaxis()->SetTitle("dE/dx Lay6");

  hdEdxVsMod->SetStats(0);
  TCanvas* cdedx=new TCanvas("cdedx","dedx SDD",1400,600);
  cdedx->SetLogz();
  hdEdxVsMod->Draw("col"); 
  hdEdxVsMod->GetXaxis()->SetTitle("SDD Module Id");
  hdEdxVsMod->GetYaxis()->SetTitle("dE/dx (keV/300 #mum)");
  hdEdxVsMod->GetYaxis()->SetTitleOffset(1.25);



  TCanvas* cv=new TCanvas("cv","Vertex",600,900);
  cv->Divide(1,3);
  cv->cd(1);
  hvx->Draw();
  hvx->GetXaxis()->SetTitle("Xv (cm)");
  cv->cd(2);
  hvy->Draw();
  hvy->GetXaxis()->SetTitle("Yv (cm)");
  cv->cd(3);
  hvz->Draw();
  hvz->GetXaxis()->SetTitle("Xv (cm)");

  hGoodPMod->SetStats(0);
  hGoodPMod->SetTitle("");
  TCanvas* ceff0=new TCanvas("ceff0","ModuleIndexInfo",1000,600);
  hGoodPMod->Draw("e");
  hGoodPMod->GetXaxis()->SetTitle("SDD Module Id");
  hGoodPMod->GetYaxis()->SetTitle("Number of tracks");
  hMissPMod->SetLineColor(2);
  hMissPMod->SetMarkerColor(2);
  hMissPMod->SetMarkerStyle(22);
  hMissPMod->SetMarkerSize(0.5);
  hMissPMod->Draw("psame");
  hBadRegMod->SetLineColor(kGreen+1);
  hBadRegMod->SetMarkerColor(kGreen+1);
  hBadRegMod->SetMarkerStyle(20);
  hBadRegMod->SetMarkerSize(0.5);
  hBadRegMod->Draw("esame");
  hSkippedMod->SetLineColor(kYellow);
  hSkippedMod->Draw("esame");
  hOutAccMod->SetLineColor(4);
  hOutAccMod->Draw("esame");
  hNoRefitMod->SetLineColor(6);
  hNoRefitMod->Draw("esame");
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
  for(Int_t imod=0; imod<260;imod++){
    Float_t numer=hGoodPMod->GetBinContent(imod+1)+hBadRegMod->GetBinContent(imod+1)+hOutAccMod->GetBinContent(imod+1)+hNoRefitMod->GetBinContent(imod+1);
    Float_t denom=hAllPMod->GetBinContent(imod+1);
    Float_t eff=0.;
    Float_t erreff=0.;
    if(denom>0){
      eff=numer/denom;
      erreff=TMath::Sqrt(eff*(1-eff)/denom);
    }
    heff->SetBinContent(imod+1,eff);
    heff->SetBinError(imod+1,erreff);
  }

  printf("---- Modules with efficiency < 90%% ----\n");
  heff->SetStats(0);
  TCanvas* ceff1=new TCanvas("ceff1","Efficiency",1000,600);
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
  ceff1->Update();



  hGoodPXloc->SetTitle("");
  hGoodPZloc->SetTitle("");
  hGoodPXloc->SetStats(0);
  hGoodPZloc->SetStats(0);
  hGoodPXloc->SetMinimum(0);
  hGoodPZloc->SetMinimum(0);
  TCanvas* ceff2=new TCanvas("ceff2","LocalCoord",1000,600);
  ceff2->Divide(2,1);
  ceff2->cd(1);
  hGoodPXloc->Draw("e");
  hGoodPXloc->GetXaxis()->SetTitle("Xlocal (cm)");
  hGoodPXloc->GetYaxis()->SetTitle("Number of tracks");
  hMissPXloc->SetLineColor(2);
  hMissPXloc->SetMarkerColor(2);
  hMissPXloc->SetMarkerStyle(22);
  hMissPXloc->SetMarkerSize(0.5);
  hMissPXloc->Draw("psame");
  hBadRegXloc->SetLineColor(kGreen+1);
  hBadRegXloc->SetMarkerColor(kGreen+1);
  hBadRegXloc->SetMarkerStyle(20);
  hBadRegXloc->SetMarkerSize(0.5);
  hBadRegXloc->Draw("psame");
  t1->Draw();
  t2->Draw();
  t3->Draw();
  ceff2->cd(2);
  hGoodPZloc->Draw("e");
  hGoodPZloc->GetXaxis()->SetTitle("Zlocal (cm)");
  hGoodPZloc->GetYaxis()->SetTitle("Number of tracks");
  hMissPZloc->SetLineColor(2);
  hMissPZloc->SetMarkerColor(2);
  hMissPZloc->SetMarkerStyle(22);
  hMissPZloc->SetMarkerSize(0.5);
  hMissPZloc->Draw("psame");
  hBadRegZloc->SetLineColor(kGreen+1);
  hBadRegZloc->SetMarkerColor(kGreen+1);
  hBadRegZloc->SetMarkerStyle(20);
  hBadRegZloc->SetMarkerSize(0.5);
  hBadRegZloc->Draw("psame");
  t1->Draw();
  t2->Draw();
  t3->Draw();
  ceff2->Update();


  
}





