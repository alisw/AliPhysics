#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLine.h>
#include <TSystem.h>
#include <TGrid.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TObjArray.h>
#include "AliCDBEntry.h"
#include "AliITSDriftSpeedArraySDD.h"
#include "AliITSgeomTGeo.h"
#endif

/*  $Id: PlotDriftSpeedSDDVsTime.C 41510 2010-06-01 09:21:24Z prino $    */


// Macro to plot the drift speed vs. time from the OCDB files 
// created from INJECTOR runs (OCDB/ITS/Calib/DriftSpeedSDD)
// Origin: F. Prino (prino@to.infn.it)

void FillErrors(Float_t errSpeed[260]);

void PlotDriftSpeedSDDVsTime(Int_t year=2011, Int_t firstRun=142600, 
			     Int_t lastRun=999999999,
			     Int_t anode=128){
  TGrid::Connect("alien:",0,0,"t");
  Float_t errSpeed[260];
  FillErrors(errSpeed);
  Int_t iAn=anode;
  if(anode>256) iAn=anode-256;

  TString cmd=Form("gbbox find \"/alice/data/%d/OCDB/ITS/Calib/DriftSpeedSDD\" \"Run*.root\" > runSpeedAlien.txt",year);
  gSystem->Exec(cmd.Data());
  FILE* listruns=fopen("runSpeedAlien.txt","r");
  Char_t filnam[200],filnamalien[200];
  TGraphErrors** gvdrvstime=new TGraphErrors*[520];
  TGraphErrors** gvdrvsrun=new TGraphErrors*[520];
  TGraphErrors** gstatusinjvstime=new TGraphErrors*[520];
  TGraphErrors** gstatusinjvsrun=new TGraphErrors*[520];


  TGraph* gGoodInjVsRun=new TGraph(0);
  gGoodInjVsRun->SetName("gGoodInjVsRun");
  TGraph* gRescaledSpeedVsRun=new TGraph(0);
  gRescaledSpeedVsRun->SetName("gRescaledSpeedVsRun");
  TGraph* gAverSpeedVsRun=new TGraph(0);
  gAverSpeedVsRun->SetName("gAverSpeedVsRun");
  TGraph* gGoodInjVsTime=new TGraph(0);
  gGoodInjVsTime->SetName("gGoodInjVsTime");
  TGraph* gRescaledSpeedVsTime=new TGraph(0);
  gRescaledSpeedVsTime->SetName("gRescaledSpeedVsTime");
  TGraph* gAverSpeedVsTime=new TGraph(0);
  gAverSpeedVsTime->SetName("gAverSpeedVsIime");

  TGraph* gGoodInjVsRunL3=new TGraph(0);
  gGoodInjVsRunL3->SetName("gGoodInjVsRunL3");
  TGraph* gRescaledSpeedVsRunL3=new TGraph(0);
  gRescaledSpeedVsRunL3->SetName("gRescaledSpeedVsRunL3");
  TGraph* gAverSpeedVsRunL3=new TGraph(0);
  gAverSpeedVsRunL3->SetName("gAverSpeedVsRunL3");
  TGraph* gGoodInjVsTimeL3=new TGraph(0);
  gGoodInjVsTimeL3->SetName("gGoodInjVsTimeL3");
  TGraph* gRescaledSpeedVsTimeL3=new TGraph(0);
  gRescaledSpeedVsTimeL3->SetName("gRescaledSpeedVsTimeL3");
  TGraph* gAverSpeedVsTimeL3=new TGraph(0);
  gAverSpeedVsTimeL3->SetName("gAverSpeedVsIimeL3");

  TGraph* gGoodInjVsRunL4=new TGraph(0);
  gGoodInjVsRunL4->SetName("gGoodInjVsRunL4");
  TGraph* gRescaledSpeedVsRunL4=new TGraph(0);
  gRescaledSpeedVsRunL4->SetName("gRescaledSpeedVsRunL4");
  TGraph* gAverSpeedVsRunL4=new TGraph(0);
  gAverSpeedVsRunL4->SetName("gAverSpeedVsRunL4");
  TGraph* gGoodInjVsTimeL4=new TGraph(0);
  gGoodInjVsTimeL4->SetName("gGoodInjVsTimeL4");
  TGraph* gRescaledSpeedVsTimeL4=new TGraph(0);
  gRescaledSpeedVsTimeL4->SetName("gRescaledSpeedVsTimeL4");
  TGraph* gAverSpeedVsTimeL4=new TGraph(0);
  gAverSpeedVsTimeL4->SetName("gAverSpeedVsIimeL4");

  TGraph* gFracGoodInjVsRun=new TGraph(0);
  gFracGoodInjVsRun->SetName("gFracGoodInjVsRun");
  TGraph* gFracRescaledSpeedVsRun=new TGraph(0);
  gFracRescaledSpeedVsRun->SetName("gFracRescaledSpeedVsRun");
  TGraph* gFracAverSpeedVsRun=new TGraph(0);
  gFracAverSpeedVsRun->SetName("gFracAverSpeedVsRun");
  TGraph* gFracGoodInjVsTime=new TGraph(0);
  gFracGoodInjVsTime->SetName("gFracGoodInjVsTime");
  TGraph* gFracRescaledSpeedVsTime=new TGraph(0);
  gFracRescaledSpeedVsTime->SetName("gFracRescaledSpeedVsTime");
  TGraph* gFracAverSpeedVsTime=new TGraph(0);
  gFracAverSpeedVsTime->SetName("gAverSpeedVsIime");

  TGraph* gFracGoodInjVsRunL3=new TGraph(0);
  gFracGoodInjVsRunL3->SetName("gFracGoodInjVsRunL3");
  TGraph* gFracRescaledSpeedVsRunL3=new TGraph(0);
  gFracRescaledSpeedVsRunL3->SetName("gFracRescaledSpeedVsRunL3");
  TGraph* gFracAverSpeedVsRunL3=new TGraph(0);
  gFracAverSpeedVsRunL3->SetName("gFracAverSpeedVsRunL3");
  TGraph* gFracGoodInjVsTimeL3=new TGraph(0);
  gFracGoodInjVsTimeL3->SetName("gFracGoodInjVsTimeL3");
  TGraph* gFracRescaledSpeedVsTimeL3=new TGraph(0);
  gFracRescaledSpeedVsTimeL3->SetName("gFracRescaledSpeedVsTimeL3");
  TGraph* gFracAverSpeedVsTimeL3=new TGraph(0);
  gFracAverSpeedVsTimeL3->SetName("gFracAverSpeedVsIimeL3");

  TGraph* gFracGoodInjVsRunL4=new TGraph(0);
  gFracGoodInjVsRunL4->SetName("gFracGoodInjVsRunL4");
  TGraph* gFracRescaledSpeedVsRunL4=new TGraph(0);
  gFracRescaledSpeedVsRunL4->SetName("gFracRescaledSpeedVsRunL4");
  TGraph* gFracAverSpeedVsRunL4=new TGraph(0);
  gFracAverSpeedVsRunL4->SetName("gFracAverSpeedVsRunL4");
  TGraph* gFracGoodInjVsTimeL4=new TGraph(0);
  gFracGoodInjVsTimeL4->SetName("gFracGoodInjVsTimeL4");
  TGraph* gFracRescaledSpeedVsTimeL4=new TGraph(0);
  gFracRescaledSpeedVsTimeL4->SetName("gFracRescaledSpeedVsTimeL4");
  TGraph* gFracAverSpeedVsTimeL4=new TGraph(0);
  gFracAverSpeedVsTimeL4->SetName("gFracAverSpeedVsIimeL4");
  
  for(Int_t iMod=0; iMod<260;iMod++){
    for(Int_t iSide=0; iSide<2; iSide++){
      Int_t index=2*iMod+iSide;
      gvdrvstime[index]=new TGraphErrors(0);    
      gvdrvstime[index]->SetTitle(Form("Module %d Side %d",iMod+240,iSide));
      gvdrvstime[index]->SetName(Form("gspmod%ds%dt",iMod+240,iSide));
      gvdrvsrun[index]=new TGraphErrors(0);    
      gvdrvsrun[index]->SetTitle(Form("Module %d Side %d",iMod+240,iSide));
      gvdrvsrun[index]->SetName(Form("gspmod%ds%dr",iMod+240,iSide));
      gstatusinjvstime[index]=new TGraphErrors(0);    
      gstatusinjvstime[index]->SetTitle(Form("Module %d Side %d",iMod+240,iSide));
      gstatusinjvstime[index]->SetName(Form("gstinmod%ds%dt",iMod+240,iSide));
      gstatusinjvsrun[index]=new TGraphErrors(0);    
      gstatusinjvsrun[index]->SetTitle(Form("Module %d Side %d",iMod+240,iSide));
      gstatusinjvsrun[index]->SetName(Form("gstinmod%ds%dr",iMod+240,iSide));
    }
  }

  Float_t driftField=(1800-45)/291/0.012;  
  Int_t nrun,nrun2,nv,ns;
  UInt_t timeZero;
  if(year==2009) timeZero=1247762992;
  else if(year==2010) timeZero=1262300400;
  else timeZero=1293836400; // 1/1/2011 at 0:00 CEST

  while(!feof(listruns)){
    fscanf(listruns,"%s\n",filnam);
    Char_t directory[100];
    sprintf(directory,"/alice/data/%d",year);
    if(!strstr(filnam,directory)) continue;
       sscanf(filnam,"/alice/data/%d/OCDB/ITS/Calib/DriftSpeedSDD/Run%d_%d_v%d_s%d.root",&year,&nrun,&nrun2,&nv,&ns);

       if(year==2009 && (nrun<85639 && nrun2> 85639)) continue;// protection for files with swapped ladders 4-5 of layer 3 
       if(year==2009 && (nrun>100000 && nv< 325)) continue; // protection for files with swapped ladder 0-1 of layer 4
  
    if(nrun<firstRun) continue;
    if(nrun>lastRun) continue;
    sprintf(filnamalien,"alien://%s",filnam);
    printf("Open file: %s\n",filnam);
    TFile *f= TFile::Open(filnamalien);
    if(f==0x0)continue;
    AliCDBEntry *ent=(AliCDBEntry*)f->Get("AliCDBEntry");
    TObjArray *drspSDD = (TObjArray *)ent->GetObject();
    
    Int_t iGoodInj=0;
    Int_t iRescaledSpeed=0;
    Int_t iAverSpeed=0;
   
    Int_t iGoodInjL3=0;
    Int_t iRescaledSpeedL3=0;
    Int_t iAverSpeedL3=0;
   
    Int_t iGoodInjL4=0;
    Int_t iRescaledSpeedL4=0;
    Int_t iAverSpeedL4=0;

    Int_t totalgoodinj=520;
    Int_t totalgoodinjL3=168;
    Int_t totalgoodinjL4=352;

    Float_t fracGoodInj=0.;
    Float_t fracGoodInjL3=0.;
    Float_t fracGoodInjL4=0.;

    Float_t fracRescaledSpeed=0.;
    Float_t fracRescaledSpeedL3=0.;
    Float_t fracRescaledSpeedL4=0.;

    Float_t fracAverSpeed=0.;
    Float_t fracAverSpeedL3=0.;
    Float_t fracAverSpeedL4=0.;
    
    AliITSDriftSpeedArraySDD *vdriftarr;
    AliITSDriftSpeedArraySDD *vdriftarr0;
    AliITSDriftSpeedArraySDD *vdriftarr1;

    UInt_t timest=0;
    Float_t timeday=0;
    Bool_t goodTime=kFALSE;

    for(Int_t iHyb=0; iHyb<520;iHyb++){
      if(!goodTime){
	vdriftarr=(AliITSDriftSpeedArraySDD*)drspSDD->At(iHyb);
	Int_t statusInj=vdriftarr->GetInjectorStatus();
	if(statusInj>0){
	  timest=vdriftarr->GetTimestamp(0);
	  if(timest>0 && timest>timeZero){
	    timeday=float(timest-timeZero)/60./60./24.;
	    goodTime=kTRUE;
	  }
	}
      }
    }

    for(Int_t iMod=0; iMod<260;iMod++){

      Int_t i0=2*iMod;
      Int_t i1=1+2*iMod;
      vdriftarr0=(AliITSDriftSpeedArraySDD*)drspSDD->At(i0);
      vdriftarr1=(AliITSDriftSpeedArraySDD*)drspSDD->At(i1);


      Int_t statusInj0=vdriftarr0->GetInjectorStatus();
      Int_t statusInj1=vdriftarr1->GetInjectorStatus();
      Int_t npt=gstatusinjvsrun[i0]->GetN();
      gstatusinjvsrun[i0]->SetPoint(npt,(Float_t)nrun,statusInj0);
      gstatusinjvsrun[i1]->SetPoint(npt,(Float_t)nrun,statusInj1);
      gstatusinjvsrun[i0]->SetPointError(npt,0,0);
      gstatusinjvsrun[i1]->SetPointError(npt,0,0);
      if(goodTime){
	Int_t npt2=gstatusinjvstime[i0]->GetN();
	gstatusinjvstime[i0]->SetPoint(npt2,timeday,statusInj0);
	gstatusinjvstime[i1]->SetPoint(npt2,timeday,statusInj1);
	gstatusinjvstime[i0]->SetPointError(npt2,0,0);
	gstatusinjvstime[i1]->SetPointError(npt2,0,0);
      }

      Float_t vdrift0=vdriftarr0->GetDriftSpeed(0,iAn);
      Float_t vdrift1=vdriftarr1->GetDriftSpeed(0,iAn);
      Float_t mob=vdrift0*1.E5/driftField;  
      Float_t temper=293.15*TMath::Power((mob/1350.),-1/2.4); 
      if(iMod==497-240) printf("Run %s   Time %d Day %f Speed=%f Temp=%f\n",filnam,timest,timeday,vdrift0,temper);

      if(statusInj0>1){
	iGoodInj++;
	if(iMod<84)iGoodInjL3++;
	else iGoodInjL4++;
	npt=gvdrvsrun[i0]->GetN();
	gvdrvsrun[i0]->SetPoint(npt,(Float_t)nrun,vdrift0);
	gvdrvsrun[i0]->SetPointError(npt,0,errSpeed[iMod]);
	if(goodTime){
	  Int_t npt2=gvdrvstime[i0]->GetN();
	  gvdrvstime[i0]->SetPoint(npt2,timeday,vdrift0);
	  gvdrvstime[i0]->SetPointError(npt2,0,errSpeed[iMod]);
	}
      }else if(statusInj0==1){
	iRescaledSpeed++;
	if(iMod<84)iRescaledSpeedL3++;
	else iRescaledSpeedL4++;
      }else{ 
	iAverSpeed++;
	if(iMod<84)iAverSpeedL3++;
	else iAverSpeedL4++;
      }
      if(statusInj1>1){ 
	iGoodInj++;
	if(iMod<84)iGoodInjL3++;
	else iGoodInjL4++;
	npt=gvdrvsrun[i1]->GetN();
	gvdrvsrun[i1]->SetPoint(npt,(Float_t)nrun,vdrift1);
	gvdrvsrun[i1]->SetPointError(npt,0,errSpeed[iMod]);
	if(goodTime){
	  Int_t npt2=gvdrvstime[i1]->GetN();
	  gvdrvstime[i1]->SetPoint(npt2,timeday,vdrift1);
	  gvdrvstime[i1]->SetPointError(npt2,0,errSpeed[iMod]);
	}
      }else if(statusInj1==1){
	iRescaledSpeed++;
	if(iMod<84)iRescaledSpeedL3++;
	else iRescaledSpeedL4++;
      }else{
	iAverSpeed++;
	if(iMod<84)iAverSpeedL3++;
	else iAverSpeedL4++;
      }
    }

    Int_t npt=gGoodInjVsRun->GetN();

    fracGoodInj=(Float_t)iGoodInj/(Float_t)totalgoodinj;
    fracGoodInjL3=(Float_t)iGoodInjL3/(Float_t)totalgoodinjL3;
    fracGoodInjL4=(Float_t)iGoodInjL4/(Float_t)totalgoodinjL4;

    fracRescaledSpeed   = (Float_t)iRescaledSpeed/(Float_t)totalgoodinj;
    fracRescaledSpeedL3 = (Float_t)iRescaledSpeedL3/(Float_t)totalgoodinjL3;
    fracRescaledSpeedL4 = (Float_t)iRescaledSpeedL3/(Float_t)totalgoodinjL4;

    fracAverSpeed   = (Float_t)iAverSpeed/(Float_t)totalgoodinj;
    fracAverSpeedL3 = (Float_t)iAverSpeedL3/(Float_t)totalgoodinjL3;
    fracAverSpeedL4 = (Float_t)iAverSpeedL4/(Float_t)totalgoodinjL4;

    gGoodInjVsRun->SetPoint(npt,(Float_t)nrun,iGoodInj);
    gRescaledSpeedVsRun->SetPoint(npt,(Float_t)nrun,iRescaledSpeed);
    gAverSpeedVsRun->SetPoint(npt,(Float_t)nrun,iAverSpeed);
	
    gGoodInjVsRunL3->SetPoint(npt,(Float_t)nrun,iGoodInjL3);
    gRescaledSpeedVsRunL3->SetPoint(npt,(Float_t)nrun,iRescaledSpeedL3);
    gAverSpeedVsRunL3->SetPoint(npt,(Float_t)nrun,iAverSpeedL3);

    gGoodInjVsRunL4->SetPoint(npt,(Float_t)nrun,iGoodInjL4);
    gRescaledSpeedVsRunL4->SetPoint(npt,(Float_t)nrun,iRescaledSpeedL4);
    gAverSpeedVsRunL4->SetPoint(npt,(Float_t)nrun,iAverSpeedL4);


    gFracGoodInjVsRun->SetPoint(npt,(Float_t)nrun,(Double_t)fracGoodInj);
    gFracRescaledSpeedVsRun->SetPoint(npt,(Float_t)nrun,(Double_t)fracRescaledSpeed);
    gFracAverSpeedVsRun->SetPoint(npt,(Float_t)nrun,(Double_t)fracAverSpeed);
	
    gFracGoodInjVsRunL3->SetPoint(npt,(Float_t)nrun,(Double_t)fracGoodInjL3);
    gFracRescaledSpeedVsRunL3->SetPoint(npt,(Float_t)nrun,(Double_t)fracRescaledSpeedL3);
    gFracAverSpeedVsRunL3->SetPoint(npt,(Float_t)nrun,(Double_t)fracAverSpeedL3);

    gFracGoodInjVsRunL4->SetPoint(npt,(Float_t)nrun,(Double_t)fracGoodInjL4);
    gFracRescaledSpeedVsRunL4->SetPoint(npt,(Float_t)nrun,(Double_t)fracRescaledSpeedL4);
    gFracAverSpeedVsRunL4->SetPoint(npt,(Float_t)nrun,(Double_t)fracAverSpeedL4);

    npt=gGoodInjVsTime->GetN();

    gGoodInjVsTime->SetPoint(npt,timeday,iGoodInj);
    gRescaledSpeedVsTime->SetPoint(npt,timeday,iRescaledSpeed);
    gAverSpeedVsTime->SetPoint(npt,timeday,iAverSpeed);

    gGoodInjVsTimeL3->SetPoint(npt,timeday,iGoodInjL3);
    gRescaledSpeedVsTimeL3->SetPoint(npt,timeday,iRescaledSpeedL3);
    gAverSpeedVsTimeL3->SetPoint(npt,timeday,iAverSpeedL3);

    gGoodInjVsTimeL4->SetPoint(npt,timeday,iGoodInjL4);
    gRescaledSpeedVsTimeL4->SetPoint(npt,timeday,iRescaledSpeedL4);
    gAverSpeedVsTimeL4->SetPoint(npt,timeday,iAverSpeedL4);


    gFracGoodInjVsTime->SetPoint(npt,timeday,(Double_t)fracGoodInj);
    gFracRescaledSpeedVsTime->SetPoint(npt,timeday,(Double_t)fracRescaledSpeed);
    gFracAverSpeedVsTime->SetPoint(npt,timeday,(Double_t)fracAverSpeed);

    gFracGoodInjVsTimeL3->SetPoint(npt,timeday,(Double_t)fracGoodInjL3);
    gFracRescaledSpeedVsTimeL3->SetPoint(npt,timeday,(Double_t)fracRescaledSpeedL3);
    gFracAverSpeedVsTimeL3->SetPoint(npt,timeday,(Double_t)fracAverSpeedL3);

    gFracGoodInjVsTimeL4->SetPoint(npt,timeday,(Double_t)fracGoodInjL4);
    gFracRescaledSpeedVsTimeL4->SetPoint(npt,timeday,(Double_t)fracRescaledSpeedL4);
    gFracAverSpeedVsTimeL4->SetPoint(npt,timeday,(Double_t)fracAverSpeedL4);

    printf("Number of half-modules with drift speed from injectors = %d\n",iGoodInj);
    printf("Number of half-modules with average drift speed        = %d\n",iAverSpeed);
    printf("Number of half-modules with drift speed from injectors L3     = %d\n",iGoodInjL3);
    printf("Number of half-modules with drift speed from golden module L3 = %d\n",iRescaledSpeedL3);
    printf("Number of half-modules with drift speed from injectors L4     = %d\n",iGoodInjL4);
    printf("Number of half-modules with drift speed from golden module L4 = %d\n",iRescaledSpeedL4);

    f->Close();
  }

  Char_t filout[100];
  sprintf(filout,"DriftSpVsTime_%d.root",year);
  TFile *ofil=new TFile(filout,"recreate");
  for(Int_t iHyb=0; iHyb<520;iHyb++){
    gvdrvstime[iHyb]->Write();
    gvdrvsrun[iHyb]->Write();
    gstatusinjvstime[iHyb]->Write();
    gstatusinjvsrun[iHyb]->Write();
  }
  gGoodInjVsRun->Write();
  gGoodInjVsRunL3->Write();
  gGoodInjVsRunL4->Write();
  gGoodInjVsTime->Write();
  gGoodInjVsTimeL3->Write();
  gGoodInjVsTimeL4->Write();
  gAverSpeedVsRun->Write();
  gAverSpeedVsRunL3->Write();
  gAverSpeedVsRunL4->Write();
  gAverSpeedVsTime->Write();
  gAverSpeedVsTimeL3->Write();
  gAverSpeedVsTimeL4->Write();
  gRescaledSpeedVsRun->Write();
  gRescaledSpeedVsRunL3->Write();
  gRescaledSpeedVsRunL4->Write();
  gRescaledSpeedVsTime->Write();
  gRescaledSpeedVsTimeL3->Write();
  gRescaledSpeedVsTimeL4->Write();
  ofil->Close();

  //  Int_t mod1=244-240;
  Int_t mod1 = 243-240;
  Int_t mod2=277-240;
  //  Int_t mod2=259-240;
//   Int_t mod2=274-240;
  Int_t mod3=327-240;
  //  Int_t mod4=453-240;
  Int_t mod4=422-240;
   //  Int_t mod4=497-240;
  Int_t lay1,lad1,det1;
  Int_t lay2,lad2,det2;
  Int_t lay3,lad3,det3;
  Int_t lay4,lad4,det4;
  AliITSgeomTGeo::GetModuleId(mod1+240,lay1,lad1,det1);
  AliITSgeomTGeo::GetModuleId(mod2+240,lay2,lad2,det2);
  AliITSgeomTGeo::GetModuleId(mod3+240,lay3,lad3,det3);
  AliITSgeomTGeo::GetModuleId(mod4+240,lay4,lad4,det4);

  gStyle->SetOptTitle(0);
  TCanvas* c0=new TCanvas("c0","Vdrift vs. time");
  c0->SetGridx();
  c0->SetGridy();
  gvdrvstime[2*mod1]->SetMarkerStyle(20);
  gvdrvstime[2*mod2]->SetMarkerStyle(22);
  gvdrvstime[2*mod2]->SetMarkerColor(2);
  gvdrvstime[2*mod2]->SetLineColor(2);
  gvdrvstime[2*mod3]->SetMarkerStyle(29);
  gvdrvstime[2*mod3]->SetMarkerColor(3);
  gvdrvstime[2*mod3]->SetLineColor(3);
  gvdrvstime[2*mod4]->SetMarkerStyle(27);
  gvdrvstime[2*mod4]->SetMarkerColor(4);
  gvdrvstime[2*mod4]->SetLineColor(4);
  gvdrvstime[2*mod1]->Draw("AP");
  gvdrvstime[2*mod1]->SetMinimum(6.3);
  gvdrvstime[2*mod1]->SetMaximum(6.75);
  Char_t title[100];
  if(year==2009){
    sprintf(title,"Time (days since July 16th 2009)");
  }else if (year==2010){
    sprintf(title,"Time (days since January 1st 2010)");
  }else{
    sprintf(title,"Time (days since January 1st 2011)");
  }
  gvdrvstime[2*mod1]->GetXaxis()->SetTitle(title);
  gvdrvstime[2*mod1]->GetYaxis()->SetTitle("Drift speed (#mum/ns)");
  gvdrvstime[2*mod2]->Draw("PSAME");
  gvdrvstime[2*mod3]->Draw("PSAME");
  gvdrvstime[2*mod4]->Draw("PSAME");
  TLegend* leg=new TLegend(0.6,0.7,0.89,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  TLegendEntry* lent=leg->AddEntry(gvdrvstime[2*mod1],Form("Lay %d Lad %d Det %d",lay1,lad1,det1),"P");
  lent=leg->AddEntry(gvdrvstime[2*mod2],Form("Lay %d Lad %d Det %d",lay2,lad2,det2),"P");
  lent->SetTextColor(2);
  lent=leg->AddEntry(gvdrvstime[2*mod3],Form("Lay %d Lad %d Det %d",lay3,lad3,det3),"P");
  lent->SetTextColor(3);
  lent=leg->AddEntry(gvdrvstime[2*mod4],Form("Lay %d Lad %d Det %d",lay4,lad4,det4),"P");
  lent->SetTextColor(4);
  leg->Draw();

  TCanvas* c1=new TCanvas("c1","Vdrift vs. run");
  c1->SetGridx();
  c1->SetGridy();
  gvdrvsrun[2*mod1]->SetMarkerStyle(20);
  gvdrvsrun[2*mod2]->SetMarkerStyle(22);
  gvdrvsrun[2*mod2]->SetMarkerColor(2);
  gvdrvsrun[2*mod2]->SetLineColor(2);
  gvdrvsrun[2*mod3]->SetMarkerStyle(29);
  gvdrvsrun[2*mod3]->SetMarkerColor(3);
  gvdrvsrun[2*mod3]->SetLineColor(3);
  gvdrvsrun[2*mod4]->SetMarkerStyle(27);
  gvdrvsrun[2*mod4]->SetMarkerColor(4);
  gvdrvsrun[2*mod4]->SetLineColor(4);
  gvdrvsrun[2*mod1]->Draw("AP");
  gvdrvsrun[2*mod1]->SetMinimum(6.3);
  gvdrvsrun[2*mod1]->SetMaximum(6.75);

  gvdrvsrun[2*mod1]->GetXaxis()->SetTitle("Run number");
  gvdrvsrun[2*mod1]->GetYaxis()->SetTitle("Drift speed (#mum/ns)");
  gvdrvsrun[2*mod2]->Draw("PSAME");
  gvdrvsrun[2*mod3]->Draw("PSAME");
  gvdrvsrun[2*mod4]->Draw("PSAME");
  leg->Draw();


  TH2F* hlay3=new TH2F("hlay3","Variation of the drift speed (%) Layer 3",6,-0.5,5.5,14,-0.5,13.5);
  hlay3->GetXaxis()->SetTitle("Detector");
  hlay3->GetYaxis()->SetTitle("Ladder");
  hlay3->GetXaxis()->SetTickLength(0);
  hlay3->GetYaxis()->SetTickLength(0);
  hlay3->SetStats(0);

  TH2F* hlay4=new TH2F("hlay4","Variation of the drift speed (%) Layer 4",8,-0.5,7.5,22,-0.5,21.5);
  hlay4->GetXaxis()->SetTitle("Detector");
  hlay4->GetYaxis()->SetTitle("Ladder");
  hlay4->GetXaxis()->SetTickLength(0);
  hlay4->GetYaxis()->SetTickLength(0);
  hlay4->GetYaxis()->SetTitle("Ladder");
  hlay4->SetStats(0);

  Double_t run1,run2,vdr1,vdr2;
  Int_t lay,lad,det;
  for(Int_t iMod=0; iMod<260; iMod++){
    Int_t lastPoint=gvdrvsrun[2*iMod]->GetN()-1;
    gvdrvsrun[2*iMod]->GetPoint(lastPoint,run2,vdr2);
    gvdrvsrun[2*iMod]->GetPoint(lastPoint-1,run1,vdr1);
    Float_t diff=0.;
    if(vdr1>0.) diff=100*(vdr2-vdr1)/vdr1;
    AliITSgeomTGeo::GetModuleId(iMod+240,lay,lad,det);
    if(lay==3) hlay3->SetBinContent(det,lad,diff);
    if(lay==4) hlay4->SetBinContent(det,lad,diff);
  }
  TLine* lin=new TLine(0,0,0,23);  
  gStyle->SetPalette(1);

  TCanvas *c0b=new TCanvas("c0b","Percent difference Last Run - Previous Run",900,600);
  c0b->Divide(2,1);
  c0b->cd(1);
  hlay3->DrawCopy("colz");
  for(Int_t i=0;i<6;i++){
    lin->SetY1(-0.5);
    lin->SetY2(13.5);
    lin->SetX1(i+0.5);
    lin->SetX2(i+0.5);
    lin->DrawClone();
  }
  for(Int_t i=0;i<14;i++){
    lin->SetX1(-0.5);
    lin->SetX2(5.5);
    lin->SetY1(i+0.5);
    lin->SetY2(i+0.5);
    lin->DrawClone();
  }
  c0b->cd(2);
  hlay4->DrawCopy("colz");
  for(Int_t i=0;i<8;i++){
    lin->SetY1(-0.5);
    lin->SetY2(21.5);
    lin->SetX1(i+0.5);
    lin->SetX2(i+0.5);
    lin->DrawClone();
  }
  for(Int_t i=0;i<22;i++){
    lin->SetX1(-0.5);
    lin->SetX2(7.5);
    lin->SetY1(i+0.5);
    lin->SetY2(i+0.5);
    lin->DrawClone();
  }

  TCanvas* c4=new TCanvas("c4","GoodMod vs. run");
  c4->SetGridx();
  c4->SetGridy();
  gGoodInjVsRun->SetMarkerStyle(20);
  gGoodInjVsRun->SetMinimum(50.);
  gGoodInjVsRun->SetMaximum(370.);
  gGoodInjVsRunL3->SetMarkerStyle(22);
  gGoodInjVsRunL3->SetMarkerColor(2);
  gGoodInjVsRunL3->SetLineColor(2);
  gGoodInjVsRunL4->SetMarkerStyle(23);
  gGoodInjVsRunL4->SetMarkerColor(4);
  gGoodInjVsRunL4->SetLineColor(4);
  gGoodInjVsRun->Draw("AP");
  gGoodInjVsRunL3->Draw("PSAME");
  gGoodInjVsRunL4->Draw("PSAME");
  gGoodInjVsRun->GetXaxis()->SetTitle("Run number");
  gGoodInjVsRun->GetYaxis()->SetTitle("Half-modules with drift speed from injectors");
  TLegend* leg2=new TLegend(0.6,0.3,0.89,0.5);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  TLegendEntry* lent2=leg2->AddEntry(gGoodInjVsRun,"All","P");
  lent2->SetTextColor(1);
  lent2=leg2->AddEntry(gGoodInjVsRunL3,"Layer 3 ","P");
  lent2->SetTextColor(2);
  lent2=leg2->AddEntry(gGoodInjVsRunL4,"Layer 4","P");
  lent2->SetTextColor(4);
  leg2->Draw();


  TCanvas* c4bis=new TCanvas("c4bis"," Frac GoodMod vs. run");
  c4bis->SetGridx();
  c4bis->SetGridy();
  gFracGoodInjVsRun->SetMarkerStyle(20);
  gFracGoodInjVsRun->SetMinimum(0.);
  gFracGoodInjVsRun->SetMaximum(0.9);
  gFracGoodInjVsRunL3->SetMarkerStyle(22);
  gFracGoodInjVsRunL3->SetMarkerColor(2);
  gFracGoodInjVsRunL3->SetLineColor(2);
  gFracGoodInjVsRunL4->SetMarkerStyle(23);
  gFracGoodInjVsRunL4->SetMarkerColor(4);
  gFracGoodInjVsRunL4->SetLineColor(4);
  gFracGoodInjVsRun->Draw("AP");
  gFracGoodInjVsRunL3->Draw("PSAME");
  gFracGoodInjVsRunL4->Draw("PSAME");
  gFracGoodInjVsRun->GetXaxis()->SetTitle("Run number");
  gFracGoodInjVsRun->GetYaxis()->SetTitle("Fraction of Half-modules with drift speed from injectors");
  gFracGoodInjVsRun->GetYaxis()->SetTitleSize(0.03);
  gFracGoodInjVsRun->GetYaxis()->SetTitleOffset(1.5);
  leg2->Draw();
  
  TCanvas* c4ter=new TCanvas("c4ter","RescaledMod vs. run");
  c4ter->SetGridx();
  c4ter->SetGridy();
  gRescaledSpeedVsRun->SetMarkerStyle(20);
  gRescaledSpeedVsRun->SetMinimum(0.);
  gRescaledSpeedVsRun->SetMaximum(120.);
  gRescaledSpeedVsRunL3->SetMarkerStyle(22);
  gRescaledSpeedVsRunL3->SetMarkerColor(2);
  gRescaledSpeedVsRunL3->SetLineColor(2);
  gRescaledSpeedVsRunL4->SetMarkerStyle(23);
  gRescaledSpeedVsRunL4->SetMarkerColor(4);
  gRescaledSpeedVsRunL4->SetLineColor(4);
  gRescaledSpeedVsRun->Draw("AP");
  gRescaledSpeedVsRunL3->Draw("PSAME");
  gRescaledSpeedVsRunL4->Draw("PSAME");
  gRescaledSpeedVsRun->GetXaxis()->SetTitle("Run number");
  gRescaledSpeedVsRun->GetYaxis()->SetTitle("Half-modules with drift speed from golden module");
  leg2->Draw();

  TCanvas* c5=new TCanvas("c5","GoodMod vs. time");
  c5->SetGridx();
  c5->SetGridy();
  gGoodInjVsTime->SetMarkerStyle(20);
  gGoodInjVsTime->SetMinimum(50.);
  gGoodInjVsTime->SetMaximum(370.);
  gGoodInjVsTimeL3->SetMarkerStyle(22);
  gGoodInjVsTimeL3->SetMarkerColor(2);
  gGoodInjVsTimeL3->SetLineColor(2);
  gGoodInjVsTimeL4->SetMarkerStyle(23);
  gGoodInjVsTimeL4->SetMarkerColor(4);
  gGoodInjVsTimeL4->SetLineColor(4);
  gGoodInjVsTime->Draw("AP");
  gGoodInjVsTimeL3->Draw("PSAME");
  gGoodInjVsTimeL4->Draw("PSAME");
  gGoodInjVsTime->GetXaxis()->SetTitle(title);
  gGoodInjVsTime->GetYaxis()->SetTitle("Half-modules with drift speed from injectors");
  leg2->Draw();
  
  TCanvas* c5bis=new TCanvas("c5bis","Frac GoodMod vs. time");
  c5bis->SetGridx();
  c5bis->SetGridy();
  gFracGoodInjVsTime->SetMarkerStyle(20);
  gFracGoodInjVsTime->SetMinimum(0.);
  gFracGoodInjVsTime->SetMaximum(0.9);
  gFracGoodInjVsTimeL3->SetMarkerStyle(22);
  gFracGoodInjVsTimeL3->SetMarkerColor(2);
  gFracGoodInjVsTimeL3->SetLineColor(2);
  gFracGoodInjVsTimeL4->SetMarkerStyle(23);
  gFracGoodInjVsTimeL4->SetMarkerColor(4);
  gFracGoodInjVsTimeL4->SetLineColor(4);
  gFracGoodInjVsTime->Draw("AP");
  gFracGoodInjVsTimeL3->Draw("PSAME");
  gFracGoodInjVsTimeL4->Draw("PSAME");
  gFracGoodInjVsTime->GetXaxis()->SetTitle(title);
  gFracGoodInjVsTime->GetYaxis()->SetTitleSize(0.03);
  gFracGoodInjVsTime->GetYaxis()->SetTitleOffset(1.5);
  gFracGoodInjVsTime->GetYaxis()->SetTitle("Fraction of Half-modules with drift speed from injectors");
  leg2->Draw();

  TCanvas* c5ter=new TCanvas("c5ter","RescaledMod vs. time");
  c5ter->SetGridx();
  c5ter->SetGridy();
  gRescaledSpeedVsTime->SetMarkerStyle(20);
  gRescaledSpeedVsTime->SetMinimum(0.);
  gRescaledSpeedVsTime->SetMaximum(120.);
  gRescaledSpeedVsTimeL3->SetMarkerStyle(22);
  gRescaledSpeedVsTimeL3->SetMarkerColor(2);
  gRescaledSpeedVsTimeL3->SetLineColor(2);
  gRescaledSpeedVsTimeL4->SetMarkerStyle(23);
  gRescaledSpeedVsTimeL4->SetMarkerColor(4);
  gRescaledSpeedVsTimeL4->SetLineColor(4);
  gRescaledSpeedVsTime->Draw("AP");
  gRescaledSpeedVsTimeL3->Draw("PSAME");
  gRescaledSpeedVsTimeL4->Draw("PSAME");
  gRescaledSpeedVsTime->GetXaxis()->SetTitle(title);
  gRescaledSpeedVsTime->GetYaxis()->SetTitleSize(0.03);
  gRescaledSpeedVsTime->GetYaxis()->SetTitleOffset(1.5);
  gRescaledSpeedVsTime->GetYaxis()->SetTitle("Half-modules with drift speed from golden module");
  leg2->Draw();
}

void FillErrors(Float_t errSpeed[260]){
  Float_t err[260]={
    0.002308,0.005120,0.004632,0.001000,0.001735,
    0.001000,0.001000,0.002667,0.004237,0.005297,
    0.001000,0.005460,0.005149,0.003921,0.001000,
    0.003906,0.001000,0.004871,0.001000,0.001000,
    0.001000,0.001000,0.002261,0.002986,0.002056,
    0.002848,0.001000,0.001777,0.002822,0.004651,
    0.001000,0.003551,0.006466,0.001000,0.002083,
    0.004531,0.001000,0.002213,0.001000,0.001000,
    0.001000,0.001000,0.001000,0.003223,0.002800,
    0.002147,0.001000,0.003364,0.001000,0.001000,
    0.002515,0.003229,0.002552,0.005765,0.002368,
    0.003473,0.002363,0.001000,0.003413,0.001000,
    0.004906,0.001000,0.004346,0.004887,0.007138,
    0.007242,0.004289,0.003970,0.002914,0.002199,
    0.001000,0.003483,0.002154,0.002914,0.003097,
    0.006034,0.003101,0.001000,0.002425,0.002651,
    0.002771,0.002409,0.002260,0.003109,0.001000,
    0.003384,0.003374,0.002212,0.004441,0.001000,
    0.001000,0.001000,0.003578,0.001000,0.001000,
    0.003517,0.003590,0.001787,0.003329,0.001000,
    0.002770,0.001000,0.004032,0.003059,0.001000,
    0.001000,0.001000,0.001000,0.001000,0.001000,
    0.001000,0.004556,0.001000,0.001000,0.001000,
    0.001000,0.001000,0.001000,0.004819,0.002100,
    0.002624,0.003784,0.003772,0.002483,0.002792,
    0.001000,0.004713,0.003214,0.003180,0.002145,
    0.002470,0.003078,0.001000,0.007131,0.002770,
    0.002533,0.001000,0.004362,0.002819,0.001000,
    0.003630,0.004215,0.002975,0.001000,0.003790,
    0.002345,0.001000,0.003999,0.004555,0.003989,
    0.001000,0.001000,0.001000,0.003136,0.002426,
    0.005144,0.002844,0.002310,0.002467,0.002503,
    0.003811,0.003440,0.004773,0.003114,0.001000,
    0.000583,0.001000,0.001000,0.003385,0.001000,
    0.001000,0.001000,0.001000,0.003108,0.002109,
    0.005325,0.003750,0.002810,0.003559,0.001000,
    0.001000,0.003262,0.003903,0.001000,0.003622,
    0.002533,0.002121,0.003733,0.005353,0.002221,
    0.004767,0.003267,0.004892,0.002152,0.003398,
    0.001000,0.003146,0.001000,0.002952,0.003310,
    0.002644,0.002573,0.001000,0.003989,0.001000,
    0.005294,0.003095,0.003479,0.002250,0.001000,
    0.001000,0.005221,0.001000,0.001653,0.004330,
    0.013188,0.007375,0.003226,0.003875,0.001000,
    0.003653,0.001000,0.002655,0.001000,0.001000,
    0.001000,0.001000,0.004718,0.001000,0.001000,
    0.001000,0.002780,0.003680,0.001000,0.002787,
    0.001000,0.004617,0.001000,0.001000,0.003231,
    0.001887,0.002090,0.003326,0.129970,0.004907,
    0.004334,0.001000,0.001000,0.003489,0.002573,
    0.002566,0.002982,0.001000,0.001000,0.003436,
    0.004016,0.003736,0.001784,0.004775,0.008090};
  for(Int_t i=0;i<260;i++) errSpeed[i]=err[i];
  
  
}
