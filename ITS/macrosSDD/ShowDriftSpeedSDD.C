#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TF1.h>
#include <TFile.h>
#include <TLine.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <TGrid.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TObjArray.h>
#include "AliCDBEntry.h"
#include "AliITSDriftSpeedArraySDD.h"
#include "AliITSDriftSpeedSDD.h"
#include "AliITSgeomTGeo.h"
#endif

// Macro to plot the calibration parameters from the OCDB file 
// created from an INJECTOR run (OCDB/ITS/Calib/DriftSpeedSDD)
// Two methods ShowDriftSpeedSDD:
//  - the first takes the name of the file to be displayed
//  - the second builds the alien path+name from run number and file version
//
// Origin: F. Prino (prino@to.infn.it)

Bool_t kNoDraw = kFALSE; // set to kTRUE to eliminate module dependent plots


void ShowDriftSpeedSDD(Char_t filnam[150]="$ALICE_ROOT/ITS/Calib/DriftSpeedSDD/Run0_9999999_v0_s0.root", Int_t firstmod=0, Int_t lastmod=260,Int_t nrun=0){
  TFile *f= TFile::Open(filnam);
  AliCDBEntry *ent=(AliCDBEntry*)f->Get("AliCDBEntry");
  TObjArray *drspSDD = (TObjArray *)ent->GetObject();
  AliITSDriftSpeedArraySDD *vdriftarr0;
  AliITSDriftSpeedArraySDD *vdriftarr1;

  TH2F* hlay3=new TH2F("hlay3","Injector Status Layer 3",12,-0.5,5.5,14,-0.5,13.5);
  hlay3->GetXaxis()->SetTitle("Detector");
  hlay3->GetYaxis()->SetTitle("Ladder");
  hlay3->GetXaxis()->SetTickLength(0);
  hlay3->GetYaxis()->SetTickLength(0);
  hlay3->SetStats(0);
  hlay3->SetMinimum(-0.01);
  hlay3->SetMaximum(7.);
  TH2F* hlay4=new TH2F("hlay4","Injector Status Layer 4",16,-0.5,7.5,22,-0.5,21.5);
  hlay4->GetXaxis()->SetTitle("Detector");
  hlay4->GetYaxis()->SetTitle("Ladder");
  hlay4->GetXaxis()->SetTickLength(0);
  hlay4->GetYaxis()->SetTickLength(0);
  hlay4->GetYaxis()->SetTitle("Ladder");
  hlay4->SetStats(0);
  hlay4->SetMinimum(-0.01);
  hlay4->SetMaximum(7.);

  TGraph **gvdr0=new TGraph*[260];
  TGraph **gvdr1=new TGraph*[260];
  TCanvas *c0=0x0;
  if(!kNoDraw)c0=new TCanvas("c0","Module Drift Speed",700,1000);

  TGraph *tempvsmod0=new TGraph(0);
  TGraph *tempvsmod1=new TGraph(0);
  TGraph *vvsmod0=new TGraph(0);
  TGraph *vvsmod1=new TGraph(0);
  TGraph *poldegvsmod0=new TGraph(0); 
  TGraph *poldegvsmod1=new TGraph(0); 
  TGraph *anmaxvsmod0=new TGraph(0); 
  TGraph *anmaxvsmod1=new TGraph(0); 
  TGraph *dvcevsmod0=new TGraph(0);
  TGraph *dvcevsmod1=new TGraph(0);
  TGraph *dveevsmod0=new TGraph(0);
  TGraph *dveevsmod1=new TGraph(0);

  char tit0[100];
  sprintf(tit0,"Temperature vs. mod. number");
  if(nrun!=0)sprintf(tit0,"Temperature vs. mod. number - Run %d",nrun);
  tempvsmod0->SetTitle(tit0);
  tempvsmod1->SetTitle(tit0);

  sprintf(tit0,"Drift Speed vs. mod. number");
  if(nrun!=0)sprintf(tit0,"Drift Speed vs. mod. number - Run %d",nrun);
  vvsmod0->SetTitle(tit0);
  vvsmod1->SetTitle(tit0);

  sprintf(tit0,"Degree of poly fit vs. mod. number");
  if(nrun!=0)sprintf(tit0,"Degree of poly fit vs. mod. number - Run %d",nrun);
  poldegvsmod0->SetTitle(tit0);
  poldegvsmod1->SetTitle(tit0);

  sprintf(tit0,"Anode with max. vdrift vs. mod. number");
  if(nrun!=0)sprintf(tit0,"Anode with max. vdrift vs. mod. number - Run %d",nrun);
  anmaxvsmod0->SetTitle(tit0);
  anmaxvsmod1->SetTitle(tit0);

  sprintf(tit0,"Delta Vdrift 128-0 vs. mod. number");
  if(nrun!=0)sprintf(tit0,"Delta Vdrift 128-0 vs. mod. number - Run %d",nrun);
  dvcevsmod0->SetTitle(tit0);
  dvcevsmod1->SetTitle(tit0);

  sprintf(tit0,"Delta Vdrift 256-0 vs. mod. number");
  if(nrun!=0)sprintf(tit0,"Delta Vdrift 256-0 vs. mod. number - Run %d",nrun);
  dveevsmod0->SetTitle(tit0);
  dveevsmod1->SetTitle(tit0);

  TF1* fPoly=new TF1("fPoly","pol3",0.,256.);
  Char_t tit[100];
  TString psnm0 = "vdriftSDD.ps[";
  TString psnm1 = "vdriftSDD.ps";
  TString psnm2 = "vdriftSDD.ps]";
  if(!kNoDraw) c0->Print(psnm0.Data());
  Int_t cntpad = 0;
  Int_t iGoodInj=0;
  Int_t iRescaledSpeed=0;
  Int_t iAverSpeed=0;
  TLatex* tleft=new TLatex(0.2,0.82,"Side 0");
  tleft->SetNDC();
  tleft->SetTextColor(1);
  TLatex* tright=new TLatex(0.2,0.75,"Side 1");
  tright->SetNDC();
  tright->SetTextColor(2);

  for(Int_t i=firstmod; i<lastmod; i++){
    Int_t iMod=i+240;
    Int_t lay,lad,det;
    AliITSgeomTGeo::GetModuleId(iMod,lay,lad,det);
    Int_t i0=2*i;
    Int_t i1=1+2*i;
    vdriftarr0=(AliITSDriftSpeedArraySDD*)drspSDD->At(i0);
    vdriftarr1=(AliITSDriftSpeedArraySDD*)drspSDD->At(i1);
    AliITSDriftSpeedSDD* vdrift0=0x0;
    if(vdriftarr0) vdrift0=vdriftarr0->GetDriftSpeedObject(0);
    AliITSDriftSpeedSDD* vdrift1=0x0;
    if(vdriftarr1) vdrift1=vdriftarr1->GetDriftSpeedObject(0);

    gvdr0[i]=new TGraph(0);
    gvdr1[i]=new TGraph(0);
    gvdr0[i]->SetMarkerStyle(7);
    gvdr1[i]->SetMarkerStyle(7);
    gvdr0[i]->SetMarkerColor(kBlack);
    gvdr0[i]->SetLineColor(kBlack);
    gvdr0[i]->SetMinimum(5.);
    gvdr1[i]->SetMinimum(5.);
    gvdr0[i]->SetMaximum(7.5);
    gvdr1[i]->SetMaximum(7.5);
    sprintf(tit,"Mod %d\n",iMod);
    gvdr0[i]->SetTitle(tit);
    gvdr1[i]->SetTitle(tit);
    gvdr1[i]->SetMarkerColor(kRed);
    gvdr1[i]->SetLineColor(kRed);

    for(Int_t iAn=0; iAn<256; iAn++){
      Float_t vel0=0;
      if(vdrift0) vel0=vdrift0->GetDriftSpeedAtAnode(iAn);
      Float_t vel1=0;
      if(vdrift1) vel1=vdrift1->GetDriftSpeedAtAnode(iAn);
      gvdr0[i]->SetPoint(iAn,(Float_t)iAn,vel0);
      gvdr1[i]->SetPoint(iAn,(Float_t)iAn,vel1);
    }
    Int_t statusInj0=vdriftarr0->GetInjectorStatus();
    Int_t statusInj1=vdriftarr1->GetInjectorStatus();
    if(statusInj0>1) iGoodInj++;
    else if(statusInj0==1) iRescaledSpeed++;
    else iAverSpeed++;
    if(statusInj1>0) iGoodInj++;
    else if(statusInj1==1) iRescaledSpeed++;
    else iAverSpeed++;

    printf(" Mod. %d \tStatusLR=%X %X \t TimeStamp=%d \t v(an 128l)= %f",iMod,statusInj0,statusInj1,vdrift0->GetEventTimestamp(),vdriftarr0->GetDriftSpeed(0,128));
    printf("        \t v(an 128r)= %f  Degree=%d %d\n",vdriftarr1->GetDriftSpeed(0,128),vdrift0->GetDegreeofPoly(),vdrift1->GetDegreeofPoly());

    Int_t n7=(statusInj0&(0x1F<<25))>>25;
    Int_t n6=(statusInj0&(0x1F<<20))>>20;
    Int_t n5=(statusInj0&(0x1F<<15))>>15;
    Int_t n4=(statusInj0&(0x1F<<5))>>10;
    Int_t n3=(statusInj0&(0x1F<<5))>>5;
    Int_t n2=statusInj0&0x1F;
    Float_t aveStatus0=(7.*n7+6.*n6+5.*n5+4.*n4+3.*n3+2.*n2)/32.;
    n7=(statusInj1&(0x1F<<25))>>25;
    n6=(statusInj1&(0x1F<<20))>>20;
    n5=(statusInj1&(0x1F<<15))>>15;
    n4=(statusInj1&(0x1F<<5))>>10;
    n3=(statusInj1&(0x1F<<5))>>5;
    n2=statusInj1&0x1F;
    Float_t aveStatus1=(7.*n7+6.*n6+5.*n5+4.*n4+3.*n3+2.*n2)/32.;

    Int_t index=1+(det-1)*2;
    if(lay==3){ 
      hlay3->SetBinContent(index,lad,aveStatus0);
      hlay3->SetBinContent(index+1,lad,aveStatus1);
    }
    if(lay==4){ 
      hlay4->SetBinContent(index,lad,aveStatus0);
      hlay4->SetBinContent(index+1,lad,aveStatus1);
    }


    if(!kNoDraw){
      if (i%12==0 ) {
	c0->cd();
	c0->Modified();
	c0->Update();
	if (i) c0->Print(psnm1.Data());
	c0->Clear();
	c0->Divide(3,4);
	cntpad = 0;
      }
      c0->cd(++cntpad);
      gvdr0[i]->Draw("AP");
      gvdr0[i]->GetXaxis()->SetTitle("Anode");
      gvdr0[i]->GetYaxis()->SetTitle("Vdrift (#mum/ns)");      
      gvdr1[i]->Draw("P same");
      gvdr1[i]->GetXaxis()->SetTitle("Anode");
      gvdr1[i]->GetYaxis()->SetTitle("Vdrift (#mum/ns)");
      tleft->Draw();
      tright->Draw();
    }

    Float_t vel0=0;
    Float_t pd0=0;
    if(vdrift0){ 
      vel0=vdrift0->GetDriftSpeedAtAnode(128);
     pd0=vdrift0->GetDegreeofPoly();
    }
    Float_t vel1=0;
    Float_t pd1=0;
    if(vdrift1){ 
      vel1=vdrift1->GetDriftSpeedAtAnode(128);
      pd1=vdrift1->GetDegreeofPoly();
    }
    Float_t Edrift=(1800-45)/291/0.012;  
    Float_t mob0=vel0*1.E5/Edrift;  
    Float_t temper0=293.15*TMath::Power((mob0/1350.),-1/2.4); 
    Float_t mob1=vel1*1.E5/Edrift;  
    Float_t temper1=293.15*TMath::Power((mob1/1350.),-1/2.4); 
    tempvsmod0->SetPoint(tempvsmod0->GetN(),(Float_t)iMod,temper0);
    tempvsmod1->SetPoint(tempvsmod1->GetN(),(Float_t)iMod,temper1);
    vvsmod0->SetPoint(vvsmod0->GetN(),(Float_t)iMod,vel0);
    vvsmod1->SetPoint(vvsmod1->GetN(),(Float_t)iMod,vel1);
    poldegvsmod0->SetPoint(poldegvsmod0->GetN(),(Float_t)iMod,pd0);
    poldegvsmod1->SetPoint(poldegvsmod1->GetN(),(Float_t)iMod,pd1);

    for(Int_t ipar=0; ipar<=vdrift0->GetDegreeofPoly(); ipar++){
      fPoly->SetParameter(ipar,vdrift0->GetDriftSpeedParameter(ipar));
    }
    if(vdrift0->GetDegreeofPoly()<3){
      for(Int_t ipar=vdrift0->GetDegreeofPoly()+1; ipar<=3; ipar++) fPoly->SetParameter(ipar,0.);
    }

    anmaxvsmod0->SetPoint(anmaxvsmod0->GetN(),(Float_t)iMod,fPoly->GetMaximumX(0.,256.));
    dvcevsmod0->SetPoint(dvcevsmod0->GetN(),(Float_t)iMod,fPoly->Eval(128)-fPoly->Eval(0));
    dveevsmod0->SetPoint(dveevsmod0->GetN(),(Float_t)iMod,fPoly->Eval(256)-fPoly->Eval(0));
    
    for(Int_t ipar=0; ipar<=vdrift1->GetDegreeofPoly(); ipar++){
      fPoly->SetParameter(ipar,vdrift1->GetDriftSpeedParameter(ipar));
    }
    if(vdrift1->GetDegreeofPoly()<3){
      for(Int_t ipar=vdrift1->GetDegreeofPoly()+1; ipar<=3; ipar++) fPoly->SetParameter(ipar,0.);
    }
    anmaxvsmod1->SetPoint(anmaxvsmod1->GetN(),(Float_t)iMod,fPoly->GetMaximumX(0.,256.));
    dvcevsmod1->SetPoint(dvcevsmod1->GetN(),(Float_t)iMod,fPoly->Eval(128)-fPoly->Eval(0));
    dveevsmod1->SetPoint(dveevsmod1->GetN(),(Float_t)iMod,fPoly->Eval(256)-fPoly->Eval(0));
    //    getchar();
  }

  if(!kNoDraw){
    c0->cd();
    c0->Modified();
    c0->Update();
    c0->Print(psnm2.Data());
  }
  printf("Number of half-modules with drift speed from injectors               = %d\n",iGoodInj);
  printf("Number of half-modules with drift speed rewscaled from golden module = %d\n",iRescaledSpeed);
  printf("Number of half-modules with average drift speed                      = %d\n",iAverSpeed);

  gStyle->SetPalette(59);

  TCanvas* clay=new TCanvas("clay","Injector Status",900,600);
  clay->Divide(2,1);
  clay->cd(1);
  hlay3->Draw("colz");
  TLine** linv3=new TLine*[5];
  for(Int_t i=0;i<5;i++){
    linv3[i]=new TLine(i+0.5,-0.5,i+0.5,13.5);
    linv3[i]->SetLineColor(kGray+1);
    linv3[i]->Draw();
  }
  TLine** linh3=new TLine*[13];
  for(Int_t i=0;i<13;i++){
    linh3[i]=new TLine(-0.5,i+0.5,5.5,i+0.5);
    linh3[i]->SetLineColor(kGray+1);
    linh3[i]->Draw();
  }
  clay->cd(2);
  hlay4->Draw("colz");
  TLine** linv4=new TLine*[7];
  for(Int_t i=0;i<7;i++){
    linv4[i]=new TLine(i+0.5,-0.5,i+0.5,21.5);
    linv4[i]->SetLineColor(kGray+1);
    linv4[i]->Draw();
  }
  TLine** linh4=new TLine*[21];
  for(Int_t i=0;i<21;i++){
    linh4[i]=new TLine(-0.5,i+0.5,7.5,i+0.5);
    linh4[i]->SetLineColor(kGray+1);
    linh4[i]->Draw();
  }
  clay->Modified();


  TCanvas* c2;
  c2=new TCanvas("c2","Vdrift vs. mod",1000,700);
  vvsmod0->SetMarkerStyle(20);
  vvsmod0->Draw("AP");
  vvsmod0->GetXaxis()->SetTitle("Module Number");
  vvsmod0->GetYaxis()->SetTitle("Vdrift (#mum/ns)");
  vvsmod1->SetMarkerStyle(21);
  vvsmod1->SetMarkerColor(2);
  vvsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();

  TCanvas* c2t;
  c2t=new TCanvas("c2t","Temper vs. mod",1000,700);
  tempvsmod0->SetMarkerStyle(20);
  tempvsmod0->Draw("AP");
  tempvsmod0->GetXaxis()->SetTitle("Module Number");
  tempvsmod0->GetYaxis()->SetTitle("Temperature (K)");
  tempvsmod1->SetMarkerStyle(21);
  tempvsmod1->SetMarkerColor(2);
  tempvsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();

  TCanvas* c3;
  c3=new TCanvas("c3","Params vs. mod",900,900);
  c3->Divide(2,2);
  
  c3->cd(1);
  gPad->SetLeftMargin(0.14);
  poldegvsmod0->SetMarkerStyle(20);
  poldegvsmod0->Draw("AP");
  poldegvsmod0->GetXaxis()->SetTitle("Module Number");
  poldegvsmod0->GetYaxis()->SetTitle("Degree of Polynomial fit");
  poldegvsmod0->GetYaxis()->SetTitleOffset(1.4);
  poldegvsmod1->SetMarkerStyle(21);
  poldegvsmod1->SetMarkerColor(2);
  poldegvsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();
  c3->cd(2);
  gPad->SetLeftMargin(0.14);
  anmaxvsmod0->SetMarkerStyle(20);
  anmaxvsmod0->Draw("AP");
  anmaxvsmod0->GetXaxis()->SetTitle("Module Number");
  anmaxvsmod0->GetYaxis()->SetTitle("Anode with max. drift speed");
  anmaxvsmod0->GetYaxis()->SetTitleOffset(1.4);
  anmaxvsmod1->SetMarkerStyle(21);
  anmaxvsmod1->SetMarkerColor(2);
  anmaxvsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();
  c3->cd(3);
  gPad->SetLeftMargin(0.14);
  dvcevsmod0->SetMarkerStyle(20);
  dvcevsmod0->Draw("AP");
  dvcevsmod0->GetXaxis()->SetTitle("Module Number");
  dvcevsmod0->GetYaxis()->SetTitle("vdrift(anode128)-vdrift(anode0)");
  dvcevsmod0->GetYaxis()->SetTitleOffset(1.4);
  dvcevsmod1->SetMarkerStyle(21);
  dvcevsmod1->SetMarkerColor(2);
  dvcevsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();
  c3->cd(4);
  gPad->SetLeftMargin(0.14);
  dveevsmod0->SetMarkerStyle(20);
  dveevsmod0->Draw("AP");
  dveevsmod0->GetYaxis()->SetTitleOffset(1.4);
  dveevsmod0->GetXaxis()->SetTitle("Module Number");
  dveevsmod0->GetYaxis()->SetTitle("vdrift(anode256)-vdrift(anode0)");
  dveevsmod1->SetMarkerStyle(21);
  dveevsmod1->SetMarkerColor(2);
  dveevsmod1->Draw("SAMEP");
  tleft->Draw();
  tright->Draw();

  
}



void ShowDriftSpeedSDD(Int_t nrun, Int_t year=2011, Int_t nv=-1){
  TGrid::Connect("alien:",0,0,"t");
  TString cmd=Form("gbbox find \"/alice/data/%d/OCDB/ITS/Calib/DriftSpeedSDD\" \"Run%d*.root\" > run.txt",year,nrun);
  if(nv>0){
    cmd.Form("gbbox find \"/alice/data/%d/OCDB/ITS/Calib/DriftSpeedSDD\" \"Run%d_999999999_v%d_s0.root\" > run.txt",year,nrun,nv);  
  }
  gSystem->Exec(cmd.Data());
  
  Char_t filnam[200],filnamalien[200];
  FILE* runtxt=fopen("run.txt","r");
  fscanf(runtxt,"%s\n",filnam);    
  if(!strstr(filnam,"/alice/data/")){
    printf("Bad run number\n");
    gSystem->Exec("rm run.txt");
    return;
  }  
  sprintf(filnamalien,"alien://%s",filnam);
  
  printf("Open file: %s\n",filnamalien);
  ShowDriftSpeedSDD(filnamalien,0,260,nrun);
  fclose(runtxt);
  gSystem->Exec("rm run.txt");
  
}
