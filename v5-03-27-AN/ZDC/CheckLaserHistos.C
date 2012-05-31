#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TROOT.h>
#include <TStyle.h>
#include <Riostream.h>
#include <TStopwatch.h>
#include <TDatime.h>
#include <TClassTable.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TFunction.h>
#include <TCanvas.h>
#include <TGrid.h>
#include <TFile.h>

#endif

void CheckLaserHistos(Int_t nRun=0, Int_t optPlot = 1)
{
  if(nRun==0){
    printf("\n\n YOU MUST PROVIDE A RUN NUMBER!!! \n\n");
    return;
  }
  
  TGrid::Connect("alien:",0,0,"t");
  
  char histoFName[150];
  sprintf(histoFName,"alien:///alice/data/2010/Reference/ZDC/%d_laserReference.root",nRun);
  
  TFile *file = TFile::Open(histoFName);
  file->cd();
  TH1F::AddDirectory(0);
  //
  // --- Histos for reference PMTs (high gain chains)
  TH1F *hPMRefChg = new TH1F("hPMRefChg","hPMRefChg", 100,0.,1000.);
  TH1F *hPMRefAhg = new TH1F("hPMRefAhg","hPMRefAhg", 100,0.,1000.);
  TH1F *hPMRefClg = new TH1F("hPMRefClg","hPMRefClg", 100,0.,4000.);
  TH1F *hPMRefAlg = new TH1F("hPMRefAlg","hPMRefAlg", 100,0.,4000.);
  //
  hPMRefChg = (TH1F*) file->Get("hPMRefChg");
  hPMRefAhg = (TH1F*) file->Get("hPMRefAhg");
  hPMRefClg = (TH1F*) file->Get("hPMRefClg");
  hPMRefAlg = (TH1F*) file->Get("hPMRefAlg");
  // --- Histos for detector PMTs 
  TH1F *hZNChg[5], *hZPChg[5], *hZNAhg[5], *hZPAhg[5], *hZEMhg[2];
  TH1F *hZNClg[5], *hZPClg[5], *hZNAlg[5], *hZPAlg[5], *hZEMlg[2];
  char hnamZNChg[20], hnamZPChg[20], hnamZNAhg[20], hnamZPAhg[20];
  char hnamZNClg[20], hnamZPClg[20], hnamZNAlg[20], hnamZPAlg[20];
  char hnamZEMhg[20], hnamZEMlg[20];
  for(Int_t j=0; j<5; j++){
    sprintf(hnamZNChg,"ZNChg-tow%d",j);
    sprintf(hnamZPChg,"ZPChg-tow%d",j);
    sprintf(hnamZNAhg,"ZNAhg-tow%d",j);
    sprintf(hnamZPAhg,"ZPAhg-tow%d",j);
    //
    hZNChg[j] = new TH1F(hnamZNChg, hnamZNChg, 100, 0., 1000.);
    hZPChg[j] = new TH1F(hnamZPChg, hnamZPChg, 100, 0., 1000.);
    hZNAhg[j] = new TH1F(hnamZNAhg, hnamZNAhg, 100, 0., 1000.);
    hZPAhg[j] = new TH1F(hnamZPAhg, hnamZPAhg, 100, 0., 1000.);
    //
    sprintf(hnamZNClg,"ZNClg-tow%d",j);
    sprintf(hnamZPClg,"ZPClg-tow%d",j);
    sprintf(hnamZNAlg,"ZNAlg-tow%d",j);
    sprintf(hnamZPAlg,"ZPAlg-tow%d",j);
    //
    hZNClg[j] = new TH1F(hnamZNClg, hnamZNClg, 100, 0., 4000.);
    hZPClg[j] = new TH1F(hnamZPClg, hnamZPClg, 100, 0., 4000.);
    hZNAlg[j] = new TH1F(hnamZNAlg, hnamZNAlg, 100, 0., 4000.);
    hZPAlg[j] = new TH1F(hnamZPAlg, hnamZPAlg, 100, 0., 4000.);
    //
    if(j<2){
      sprintf(hnamZEMhg,"ZEM%dhg",j);
      sprintf(hnamZEMlg,"ZEM%dlg",j);
      //
      hZEMhg[j] = new TH1F(hnamZEMhg, hnamZEMhg, 100, 0., 1000.);      
      hZEMlg[j] = new TH1F(hnamZEMlg, hnamZEMlg, 100, 0., 4000.);      
    }
    //
    hZNChg[j] = (TH1F*) file->Get(hnamZNChg);
    hZPChg[j] = (TH1F*) file->Get(hnamZPChg);
    hZNAhg[j] = (TH1F*) file->Get(hnamZNAhg);
    hZPAhg[j] = (TH1F*) file->Get(hnamZPAhg);
    //
    hZNClg[j] = (TH1F*) file->Get(hnamZNClg);
    hZPClg[j] = (TH1F*) file->Get(hnamZPClg);
    hZNAlg[j] = (TH1F*) file->Get(hnamZNAlg);
    hZPAlg[j] = (TH1F*) file->Get(hnamZPAlg);
    
    
  }

 if(optPlot==1){
  // Plot the retrieved histos
  //***********************************************************
  // #### ROOT initialization
  gROOT->Reset();
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  gStyle->SetTitleTextColor(9);
  gStyle->SetStatTextColor(4);
  gStyle->SetLineColor(1);
  gStyle->SetPalette(1);
  //***********************************************************
  TCanvas *c1 = new TCanvas("c1","ZNC",0,0,1000,400);
  c1->Divide(5,2);
  for(Int_t y=0; y<5; y++){
    c1->cd(y+1);
    hZNChg[y]->Draw();
    c1->cd(y+6);
    hZNClg[y]->Draw();
  }
  c1->Print("ZNCLaser.ps");
  //
  TCanvas *c2 = new TCanvas("c2","ZPC",300,0,1000,400);
  c2->Divide(5,2);
  for(Int_t y=0; y<5; y++){
    c2->cd(y+1);
    hZPChg[y]->Draw();
    c2->cd(y+6);
    hZPClg[y]->Draw();
  }
  c2->Print("ZPCLaser.ps");
  //
  TCanvas *c3 = new TCanvas("c3","ZEM",400,0,400,400);
  c3->Divide(2,2);
  for(Int_t y=0; y<2; y++){
    c3->cd(y+1);
    hZEMhg[y]->Draw();
    c3->cd(y+3);
    hZEMlg[y]->Draw();
  }
  c3->Print("ZEMLaser.ps");
  //
  TCanvas *c4 = new TCanvas("c4","ZNA",0,400,1000,400);
  c4->Divide(5,2);
  for(Int_t y=0; y<5; y++){
    c4->cd(y+1);
    hZNAhg[y]->Draw();
    c4->cd(y+6);
    hZNAlg[y]->Draw();
  }
  c4->Print("ZNALaser.ps");
  //
  TCanvas *c5 = new TCanvas("c5","ZPA",300,400,1000,400);
  c5->Divide(5,2);
  for(Int_t y=0; y<5; y++){
    c5->cd(y+1);
    hZPAhg[y]->Draw();
    c5->cd(y+6);
    hZPAlg[y]->Draw();
  }
  c5->Print("ZPALaser.ps");
  //
  TCanvas *c6 = new TCanvas("c6","Ref",400,0,400,400);
  c6->Divide(2,2);
  c6->cd(1);
  hPMRefChg->Draw();
  c6->cd(2);
  hPMRefAhg->Draw();
  c6->cd(3);
  hPMRefClg->Draw();
  c6->cd(4);
  hPMRefAlg->Draw();
 //
  c6->Print("RefLaser.ps");
 }
}
