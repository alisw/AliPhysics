#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
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

/*  $Id$    */


// Macro to plot the drift speed vs. time from the OCDB files 
// created from INJECTOR runs (OCDB/ITS/Calib/DriftSpeedSDD)
// Origin: F. Prino (prino@to.infn.it)

void FillErrors(Float_t errSpeed[260]);

void PlotDriftSpeedSDDVsTime(Int_t firstRun=76172, 
			     Int_t lastRun=999999999, 
			     Int_t anode=128){

  TGrid::Connect("alien:",0,0,"t");

  Float_t errSpeed[260];
  FillErrors(errSpeed);

  gSystem->Exec("gbbox find \"/alice/data/2009/OCDB/ITS/Calib/DriftSpeedSDD\" \"Run*.root\" > runSpeedAlien.txt");
  FILE* listruns=fopen("runSpeedAlien.txt","r");
  Char_t filnam[200],filnamalien[200];
  TGraphErrors** gvdrvstime=new TGraphErrors*[520];
  TGraphErrors** gvdrvsrun=new TGraphErrors*[520];
  Int_t iPt[260];
  for(Int_t iMod=0; iMod<260;iMod++){
    gvdrvstime[iMod]=new TGraphErrors(0);    
    gvdrvstime[iMod]->SetTitle(Form("Module %d",iMod+240));
    gvdrvstime[iMod]->SetName(Form("gspmod%dt",iMod+240));
    gvdrvsrun[iMod]=new TGraphErrors(0);    
    gvdrvsrun[iMod]->SetTitle(Form("Module %d",iMod+240));
    gvdrvsrun[iMod]->SetName(Form("gspmod%dr",iMod+240));
    iPt[iMod]=0;
  }
  Float_t Edrift=(1800-45)/291/0.012;  
  Int_t nrun,nrun2,dum;

  while(!feof(listruns)){
    fscanf(listruns,"%s\n",filnam);
    if(!strstr(filnam,"/alice/data/2009")) continue;
    sscanf(filnam,"/alice/data/2009/OCDB/ITS/Calib/DriftSpeedSDD/Run%d_%d_v%d_s%d.root",&nrun,&nrun2,&dum,&dum);
    if(nrun<85639 && nrun2> 85639) continue;
    if(nrun<firstRun) continue;
    if(nrun>lastRun) continue;
    sprintf(filnamalien,"alien://%s",filnam);
    printf("Open file: %s\n",filnam);
    TFile *f= TFile::Open(filnamalien);
    AliCDBEntry *ent=(AliCDBEntry*)f->Get("AliCDBEntry");
    TObjArray *drspSDD = (TObjArray *)ent->GetObject();
    
    AliITSDriftSpeedArraySDD *vdriftarr;
    for(Int_t iMod=0; iMod<260;iMod++){
      Int_t index=-1;
      if(anode<256) index=2*iMod;
      else index=2*iMod+1;
      vdriftarr=(AliITSDriftSpeedArraySDD*)drspSDD->At(index);
      Int_t iAn=anode;
      if(anode>256) iAn=anode-256;
      Float_t vdrift=vdriftarr->GetDriftSpeed(0,iAn);
      UInt_t timest=vdriftarr->GetTimestamp(0);
      if(timest==0) continue;
      Float_t timeday=float(timest-1247762992)/60./60./24.;
      Float_t mob=vdrift*1.E5/Edrift;  
      Float_t temper=293.15*TMath::Power((mob/1350.),-1/2.4); 
      if(iMod==497-240) printf("Run %s   Time %d Day %f Speed=%f Temp=%f\n",filnam,timest,timeday,vdrift,temper);
      gvdrvstime[iMod]->SetPoint(iPt[iMod],timeday,vdrift);
      gvdrvstime[iMod]->SetPointError(iPt[iMod],0,errSpeed[iMod]);
      gvdrvsrun[iMod]->SetPoint(iPt[iMod],(Float_t)nrun,vdrift);
      gvdrvsrun[iMod]->SetPointError(iPt[iMod],0,errSpeed[iMod]);
      ++iPt[iMod];
    }
    f->Close();
  }

  Int_t mod1=244-240;
  Int_t mod2=277-240;
  //Int_t mod1=268-240;
  //Int_t mod2=274-240;
  Int_t mod3=327-240;
  Int_t mod4=453-240;
  //  Int_t mod4=497-240;
  Int_t lay1,lad1,det1;
  Int_t lay2,lad2,det2;
  Int_t lay3,lad3,det3;
  Int_t lay4,lad4,det4;
  AliITSgeomTGeo::GetModuleId(mod1+240,lay1,lad1,det1);
  AliITSgeomTGeo::GetModuleId(mod2+240,lay2,lad2,det2);
  AliITSgeomTGeo::GetModuleId(mod3+240,lay3,lad3,det3);
  AliITSgeomTGeo::GetModuleId(mod4+240,lay4,lad4,det4);

  TFile *ofil=new TFile("DriftSp2009VsTime.root","recreate");
  for(Int_t iMod=0; iMod<260;iMod++){
    gvdrvstime[iMod]->Write();
    gvdrvsrun[iMod]->Write();
  }
  ofil->Close();

  gStyle->SetOptTitle(0);
  TCanvas* c0=new TCanvas();
  c0->SetGridx();
  c0->SetGridy();
  gvdrvstime[mod1]->SetMarkerStyle(20);
  gvdrvstime[mod2]->SetMarkerStyle(22);
  gvdrvstime[mod2]->SetMarkerColor(2);
  gvdrvstime[mod2]->SetLineColor(2);
  gvdrvstime[mod3]->SetMarkerStyle(29);
  gvdrvstime[mod3]->SetMarkerColor(3);
  gvdrvstime[mod3]->SetLineColor(3);
  gvdrvstime[mod4]->SetMarkerStyle(27);
  gvdrvstime[mod4]->SetMarkerColor(4);
  gvdrvstime[mod4]->SetLineColor(4);
  gvdrvstime[mod1]->Draw("AP");
  gvdrvstime[mod1]->SetMinimum(6.3);
  gvdrvstime[mod1]->SetMaximum(6.75);

  gvdrvstime[mod1]->GetXaxis()->SetTitle("Time (days since July 16th 2009)");
  gvdrvstime[mod1]->GetYaxis()->SetTitle("Drift speed (#mum/ns)");
  gvdrvstime[mod2]->Draw("PSAME");
  gvdrvstime[mod3]->Draw("PSAME");
  gvdrvstime[mod4]->Draw("PSAME");
  TLegend* leg=new TLegend(0.6,0.7,0.89,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  TLegendEntry* lent=leg->AddEntry(gvdrvstime[mod1],Form("Lay %d Lad %d Det %d",lay1,lad1,det1),"P");
  lent=leg->AddEntry(gvdrvstime[mod2],Form("Lay %d Lad %d Det %d",lay2,lad2,det2),"P");
  lent->SetTextColor(2);
  lent=leg->AddEntry(gvdrvstime[mod3],Form("Lay %d Lad %d Det %d",lay3,lad3,det3),"P");
  lent->SetTextColor(3);
  lent=leg->AddEntry(gvdrvstime[mod4],Form("Lay %d Lad %d Det %d",lay4,lad4,det4),"P");
  lent->SetTextColor(4);
  leg->Draw();

  TCanvas* c1=new TCanvas();
  c1->SetGridx();
  c1->SetGridy();
  gvdrvsrun[mod1]->SetMarkerStyle(20);
  gvdrvsrun[mod2]->SetMarkerStyle(22);
  gvdrvsrun[mod2]->SetMarkerColor(2);
  gvdrvsrun[mod2]->SetLineColor(2);
  gvdrvsrun[mod3]->SetMarkerStyle(29);
  gvdrvsrun[mod3]->SetMarkerColor(3);
  gvdrvsrun[mod3]->SetLineColor(3);
  gvdrvsrun[mod4]->SetMarkerStyle(27);
  gvdrvsrun[mod4]->SetMarkerColor(4);
  gvdrvsrun[mod4]->SetLineColor(4);
  gvdrvsrun[mod1]->Draw("AP");
  gvdrvsrun[mod1]->SetMinimum(6.3);
  gvdrvsrun[mod1]->SetMaximum(6.75);

  gvdrvsrun[mod1]->GetXaxis()->SetTitle("Run number");
  gvdrvsrun[mod1]->GetYaxis()->SetTitle("Drift speed (#mum/ns)");
  gvdrvsrun[mod2]->Draw("PSAME");
  gvdrvsrun[mod3]->Draw("PSAME");
  gvdrvsrun[mod4]->Draw("PSAME");
  leg->Draw();
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
