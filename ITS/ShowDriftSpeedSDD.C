#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TStyle.h>
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

void ShowDriftSpeedSDD(Char_t filnam[150]="$ALICE_ROOT/ITS/Calib/DriftSpeedSDD/Run0_9999999_v0_s0.root", Int_t firstmod=0, Int_t lastmod=260){
  TFile *f= TFile::Open(filnam);
  AliCDBEntry *ent=(AliCDBEntry*)f->Get("AliCDBEntry");
  TObjArray *drspSDD = (TObjArray *)ent->GetObject();
  AliITSDriftSpeedArraySDD *vdriftarr0;
  AliITSDriftSpeedArraySDD *vdriftarr1;
  TGraph **gvdr0=new TGraph*[260];
  TGraph **gvdr1=new TGraph*[260];

  TCanvas *c0=new TCanvas("c0","Module Drift Speed",1100,500);

  TGraph *vvsmod0=new TGraph(0);
  TGraph *vvsmod1=new TGraph(0);
  vvsmod0->SetTitle("Drift Speed vs. mod. number");
  vvsmod1->SetTitle("Drift Speed vs. mod. number");
  Char_t tit[100];
  for(Int_t i=firstmod; i<lastmod; i++){
    c0->Clear();
    c0->Divide(2,1);
    Int_t i0=2*i;
    Int_t i1=1+2*i;
    vdriftarr0=(AliITSDriftSpeedArraySDD*)drspSDD->At(i0);
    vdriftarr1=(AliITSDriftSpeedArraySDD*)drspSDD->At(i1);
    
    gvdr0[i]=new TGraph(0);
    gvdr1[i]=new TGraph(0);
    gvdr0[i]->SetMarkerStyle(7);
    gvdr1[i]->SetMarkerStyle(7);
    gvdr0[i]->SetMinimum(5.);
    gvdr1[i]->SetMinimum(5.);
    gvdr0[i]->SetMaximum(9.);
    gvdr1[i]->SetMaximum(9.);
    sprintf(tit,"Mod %d\n",i);
    gvdr0[i]->SetTitle(tit);
    gvdr1[i]->SetTitle(tit);

    for(Int_t iAn=0; iAn<256; iAn++){
      Float_t vel0=0;
      if(vdriftarr0) vel0=vdriftarr0->GetDriftSpeed(1,iAn);
      Float_t vel1=0;
      if(vdriftarr1) vel1=vdriftarr1->GetDriftSpeed(1,iAn);
      gvdr0[i]->SetPoint(iAn,(Float_t)iAn,vel0);
      gvdr1[i]->SetPoint(iAn,(Float_t)iAn,vel1);
    }
    printf(" Mod. %d \t v(an 128l)= %f\n",i,vdriftarr0->GetDriftSpeed(0,128));
    printf("        \t v(an 128r)= %f\n",vdriftarr1->GetDriftSpeed(0,128));
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
    if(vdriftarr0) vel0=vdriftarr0->GetDriftSpeed(0,128);
    Float_t vel1=0;
    if(vdriftarr1) vel1=vdriftarr1->GetDriftSpeed(0,128);
    vvsmod0->SetPoint(vvsmod0->GetN(),(Float_t)i,vel0);
    vvsmod1->SetPoint(vvsmod1->GetN(),(Float_t)i,vel1);
    
    //    getchar();
  }
  TCanvas* c2;
  c2=new TCanvas("c2","",1000,700);
  vvsmod0->SetMarkerStyle(20);
  vvsmod0->Draw("AP");
  vvsmod0->GetXaxis()->SetTitle("Module Number");
  vvsmod0->GetYaxis()->SetTitle("Vdrift (#mum/ns)");
  vvsmod1->SetMarkerStyle(21);
  vvsmod1->SetMarkerColor(2);
  vvsmod1->Draw("SAMEP");
  TLatex* tleft=new TLatex(0.7,0.82,"Side 0");
  tleft->SetNDC();
  tleft->SetTextColor(1);
  tleft->Draw();
  TLatex* tright=new TLatex(0.7,0.75,"Side 1");
  tright->SetNDC();
  tright->SetTextColor(2);
  tright->Draw();

}



void ShowDriftSpeedSDD(Int_t nrun, Int_t nv){
  TGrid::Connect("alien:",0,0,"t");
  Char_t filnam[200];
  sprintf(filnam,"alien:///alice/data/2008/LHC08c/OCDB/ITS/Calib/DriftSpeedSDD/Run%d_999999999_v%d_s0.root",nrun,nv);
  printf("Open file: %s\n",filnam);
  ShowDriftSpeedSDD(filnam,0,260);
}
