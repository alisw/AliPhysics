#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TH1F.h>
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
  TGraph **gvdr0=new TGraph*[260];
  TGraph **gvdr1=new TGraph*[260];
  TCanvas *c0=NULL;
  if(!kNoDraw)c0=new TCanvas("c0","Module Drift Speed",1100,500);

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

  TF1* fPoly=new TF1("fPoly","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.,256.);
  Char_t tit[100];
  Int_t iGoodInj=0;
  Int_t iAverSpeed=0;
  for(Int_t i=firstmod; i<lastmod; i++){
    Int_t iMod=i+240;
    if(!kNoDraw){
    }
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
    gvdr0[i]->SetMinimum(5.);
    gvdr1[i]->SetMinimum(5.);
    gvdr0[i]->SetMaximum(9.);
    gvdr1[i]->SetMaximum(9.);
    sprintf(tit,"Mod %d\n",iMod);
    gvdr0[i]->SetTitle(tit);
    gvdr1[i]->SetTitle(tit);

    for(Int_t iAn=0; iAn<256; iAn++){
      Float_t vel0=0;
      if(vdrift0) vel0=vdrift0->GetDriftSpeedAtAnode(iAn);
      Float_t vel1=0;
      if(vdrift1) vel1=vdrift1->GetDriftSpeedAtAnode(iAn);
      gvdr0[i]->SetPoint(iAn,(Float_t)iAn,vel0);
      gvdr1[i]->SetPoint(iAn,(Float_t)iAn,vel1);
    }
    if(vdriftarr0->GetInjectorStatus()>0) iGoodInj++;
    else iAverSpeed++;
    if(vdriftarr1->GetInjectorStatus()>0) iGoodInj++;
    else iAverSpeed++;

    printf(" Mod. %d \tStatusLR=%X %X \t v(an 128l)= %f",iMod,vdriftarr0->GetInjectorStatus(),vdriftarr1->GetInjectorStatus(),vdriftarr0->GetDriftSpeed(0,128));
    printf("        \t v(an 128r)= %f  Degree=%d %d\n",vdriftarr1->GetDriftSpeed(0,128),vdrift0->GetDegreeofPoly(),vdrift1->GetDegreeofPoly());
    c0->Clear();
    c0->Divide(2,1);
    c0->cd(1);
    gvdr0[i]->Draw("AP");
    gvdr0[i]->GetXaxis()->SetTitle("Anode");
    gvdr0[i]->GetYaxis()->SetTitle("Vdrift (#mum/ns)");
    c0->cd(2);
    gvdr1[i]->Draw("AP");
    gvdr1[i]->GetXaxis()->SetTitle("Anode");
    gvdr1[i]->GetYaxis()->SetTitle("Vdrift (#mum/ns)");
    c0->Update();
    
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

  printf("Number of half-modules with drift speed from injectors = %d\n",iGoodInj);
  printf("Number of half-modules with average drift speed        = %d\n",iAverSpeed);

  TCanvas* c2;
  c2=new TCanvas("c2","",1000,700);
  vvsmod0->SetMarkerStyle(20);
  vvsmod0->Draw("AP");
  vvsmod0->GetXaxis()->SetTitle("Module Number");
  vvsmod0->GetYaxis()->SetTitle("Vdrift (#mum/ns)");
  vvsmod1->SetMarkerStyle(21);
  vvsmod1->SetMarkerColor(2);
  vvsmod1->Draw("SAMEP");
  TLatex* tleft=new TLatex(0.2,0.82,"Side 0");
  tleft->SetNDC();
  tleft->SetTextColor(1);
  tleft->Draw();
  TLatex* tright=new TLatex(0.2,0.75,"Side 1");
  tright->SetNDC();
  tright->SetTextColor(2);
  tright->Draw();

  TCanvas* c3;
  c3=new TCanvas("c3","",900,900);
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



void ShowDriftSpeedSDD(Int_t nrun, Int_t year=2009){
  TGrid::Connect("alien:",0,0,"t");
  TString cmd=Form("gbbox find \"/alice/data/%d/OCDB/ITS/Calib/DriftSpeedSDD\" \"Run%d*.root\" > run.txt",year,nrun);
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
