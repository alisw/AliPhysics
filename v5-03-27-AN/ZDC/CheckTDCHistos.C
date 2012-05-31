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

void CheckTDCHistos(Int_t nRun=0,  Bool_t optPlot=kTRUE)
{
  if(nRun==0){
    printf("\n\n YOU MUST PROVIDE A RUN NUMBER!!! \n\n");
    return;
  }
  
  TGrid::Connect("alien:",0,0,"t");
  
  char histoFName[150];
  sprintf(histoFName,"alien:///alice/data/2011/Reference/ZDC/%d_tdcReference.root",nRun);

  TFile *file = TFile::Open(histoFName);
  file->cd();
  TH1F::AddDirectory(0);
  //
  TH1F *hTDC[6];
  for(Int_t it=0; it<6; it++){
    if(it==0)      hTDC[it] = dynamic_cast<TH1F*> (file->Get("TDCZNC"));
    else if(it==1) hTDC[it] = dynamic_cast<TH1F*> (file->Get("TDCZNA"));
    else if(it==2) hTDC[it] = dynamic_cast<TH1F*> (file->Get("TDCZPC"));
    else if(it==3) hTDC[it] = dynamic_cast<TH1F*> (file->Get("TDCZPA"));
    else if(it==4) hTDC[it] = dynamic_cast<TH1F*> (file->Get("TDCZEM1"));
    else if(it==5) hTDC[it] = dynamic_cast<TH1F*> (file->Get("TDCZEM2"));
  }
  
  
 if(optPlot){
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
  TCanvas *c6 = new TCanvas("c6","Side C correlations",0,200,1000,800);
  c6->Divide(3,2);
  for(Int_t t=0; t<6; t++){
    c6->cd(t+1); gPad->SetLogy(1);
    hTDC[t]->SetLineColor(kAzure+t);
    hTDC[t]->Draw();
  }
  char psname[16];
  sprintf(psname,"TDCrun%d.gif",nRun);
  c6->Print(psname);
 }
}
