#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGrid.h>
#include <TFile.h>
#include "AliCDBEntry.h"
#include "AliCDBGrid.h"
#include "AliCDBId.h"
#include "AliCDBLocal.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBPath.h"
#include "AliCDBRunRange.h"
#include "AliCDBStorage.h"
#include "AliDCSValue.h"
#include "AliZDC.h"
#include "AliZDCv3.h"
#include "AliZDCPedestals.h"
#include "AliZDCEnCalib.h"
#include "AliZDCTowerCalib.h"
#include "AliZDCChMap.h"
#include "AliZDCLaserCalib.h"
#include "AliZDCMBCalib.h"
#include "AliZDCDataDCS.h"

#endif


void PlotPedestalsvsTime(Int_t year=2011, Int_t firstRun=141820, 
		   Int_t lastRun=146900, Int_t ipedGainChain=0)
{


  TGrid::Connect("alien:",0,0,"t");
  gSystem->Exec(Form("gbbox find \"/alice/data/%d/OCDB/ZDC/Calib/Pedestals/\" \"Run*.root\" > calibAlienFiles.txt",year));
  FILE* listruns=fopen("calibAlienFiles.txt","r");
  
  const int kNchannels=24;
  TGraphErrors* graph[24];
  for(Int_t i=0; i<kNchannels; i++){
     graph[i] = new TGraphErrors(0);
     char name[50], title[50];
     sprintf(name,"graph%d",i); sprintf(title,"Pedestal ch.%d vs. run#",i);
     graph[i]->SetName("graph");  graph[i]->SetTitle("title");
  }

  Char_t filnam[200], filnamalien[200];
  Int_t iPoint=0;
  Int_t nrun, nrun2, nv, ns;

  while(!feof(listruns)){
    int st = fscanf(listruns,"%s\n",filnam);    
    Char_t directory[100];
    sprintf(directory,"/alice/data/%d",year);
    if(!strstr(filnam,directory)) continue;
    sscanf(filnam,"/alice/data/%d/OCDB/ZDC/Calib/Pedestals/Run%d_%d_v%d_s%d.root",&year,&nrun,&nrun2,&nv,&ns);
    if(nrun<firstRun) continue;
    if(nrun>lastRun) continue;
    sprintf(filnamalien,"alien://%s",filnam);
    printf("Opening file: %s\n",filnam);
    TFile *f = TFile::Open(filnamalien);  
    AliCDBEntry *entry = (AliCDBEntry*)f->Get("AliCDBEntry");
    AliZDCPedestals *calibdata = dynamic_cast<AliZDCPedestals*>  (entry->GetObject());
    
    for(int i=0; i<kNchannels; i++){
      if(ipedGainChain==0){
        graph[i]->SetPoint(iPoint, (Double_t)nrun, calibdata->GetMeanPed(i));
        graph[i]->SetPointError(iPoint, 0., calibdata->GetMeanPedWidth(i));
      }
      else{
        graph[i]->SetPoint(iPoint, (Double_t)nrun, calibdata->GetMeanPed(i+kNchannels));
        graph[i]->SetPointError(iPoint, 0., calibdata->GetMeanPedWidth(i+kNchannels));
      }
    }
    iPoint++;
    f->Close();
 }

 TFile *outfile=new TFile(Form("Calib%dVsTime.root",year),"recreate");
 outfile->cd();
 for(int i=0; i<kNchannels; i++) graph[i]->Write();
 outfile->Close();

 //***********************************************************
 // #### ROOT initialization
 gROOT->Reset();
 gStyle->SetCanvasColor(10);
 gStyle->SetFrameFillColor(10);
 gStyle->SetOptTitle(0);
 gStyle->SetOptStat(1111);
 gStyle->SetOptFit(0);
 gStyle->SetTitleTextColor(4);
 gStyle->SetStatTextColor(4);
 gStyle->SetStatX(0.92);
 gStyle->SetStatY(0.92);
 gStyle->SetLineColor(1);
 gStyle->SetPalette(1);
 gStyle->SetPadTopMargin(0.05);
 gStyle->SetPadRightMargin(0.05);
 gStyle->SetPadBottomMargin(0.09);
 gStyle->SetPadLeftMargin(0.09); 
 gStyle->SetTitleOffset(1.1,"Y");  
 // *************************************************************

 TCanvas *cHadPeds = new TCanvas("cHadPeds","Hadronic ZDC pedestals",0,0,1000,800);
 cHadPeds->Divide(5,4);
 for(int ic=0; ic<5; ic++){
   // *** ZNC pedestals
   cHadPeds->cd(ic+1);
   //
   TH1F *haxis1=0;
   if(ipedGainChain==0){
     if(ic==0) haxis1 = gPad->DrawFrame(firstRun-100, 80, lastRun+100, 100);
     else  haxis1 = gPad->DrawFrame(firstRun-100, 50, lastRun+100, 70);
   }
   else{
     if(ic==0) haxis1 = gPad->DrawFrame(firstRun-100, 500, lastRun+100, 800);
     else  haxis1 = gPad->DrawFrame(firstRun-100, 300, lastRun+100, 600);
   }
   haxis1->GetXaxis()->SetNoExponent();
   haxis1->SetXTitle("RUN no.");
   haxis1->SetYTitle("ZNC pedestals");
   //
   graph[ic]->SetMarkerStyle(20);
   graph[ic]->SetMarkerColor(kBlue);
   graph[ic]->Draw("P, SAME");
   // *** ZPC pedestals
   cHadPeds->cd(ic+6);
   //
   TH1F *haxis2=0;
   if(ipedGainChain==0) haxis2= gPad->DrawFrame(firstRun-100, 55, lastRun+100, 85);
   else  haxis2 = gPad->DrawFrame(firstRun-100, 400, lastRun+100, 700);
   haxis2->GetXaxis()->SetNoExponent();
   haxis2->SetXTitle("RUN no.");
   haxis2->SetYTitle("ZPC pedestals");
   //
   graph[ic+5]->SetMarkerStyle(21);
   graph[ic+5]->SetMarkerColor(kBlue+3);
   graph[ic+5]->Draw("P, SAME");
   // *** ZNA pedestals
   cHadPeds->cd(ic+11);
   //
   TH1F *haxis3=0;
   if(ipedGainChain==0) haxis3 = gPad->DrawFrame(firstRun-100, 35, lastRun+100, 85);
   else  haxis3 = gPad->DrawFrame(firstRun-100, 300, lastRun+100, 700);
   haxis3->GetXaxis()->SetNoExponent();
   haxis3->SetXTitle("RUN no.");
   haxis3->SetYTitle("ZNA pedestals");
   //
   graph[ic+12]->SetMarkerStyle(20);
   graph[ic+12]->SetMarkerColor(kRed);
   graph[ic+12]->Draw("P, SAME");
   // *** ZPA pedestals
   cHadPeds->cd(ic+16);
   //
   TH1F *haxis4=0;
   if(ipedGainChain==0) haxis4 = gPad->DrawFrame(firstRun-100, 40, lastRun+100, 80);
   else  haxis4 = gPad->DrawFrame(firstRun-100, 300, lastRun+100, 600);
   haxis4->GetXaxis()->SetNoExponent();
   haxis4->SetXTitle("RUN no.");
   haxis4->SetYTitle("ZPA pedestals");
   //
   graph[ic+17]->SetMarkerStyle(21);
   graph[ic+17]->SetMarkerColor(kRed+1);
   graph[ic+17]->Draw("P, SAME");
 }
 cHadPeds->SaveAs("ZDCPedvsTime1.gif");
 cHadPeds->SaveAs("ZDCPedvsTime1.C");

 TCanvas *cothPeds = new TCanvas("cothPeds","ZEM + Ref. pedestals",800,0,600,600);
 cothPeds->Divide(2,2);
 for(int ic=0; ic<2; ic++){
    // *** ZEM pedestals
    cothPeds->cd(ic+1);
    //
    TH1F *haxis5=0;
    if(ipedGainChain==0) haxis5 = gPad->DrawFrame(firstRun-100, 30, lastRun+20, 70);
    else  haxis5 = gPad->DrawFrame(firstRun-100, 250, lastRun+100, 550);
    haxis5->GetXaxis()->SetNoExponent();
    haxis5->SetXTitle("RUN no.");
    haxis5->SetYTitle("ZEM pedestals");
    //
    graph[ic+10]->SetMarkerStyle(22);
    graph[ic+10]->SetMarkerColor(kGreen+1);
    graph[ic+10]->Draw("P, SAME");
    // *** Ref. pedestals
    cothPeds->cd(ic+3);
    //
    TH1F *haxis6=0; 
    if(ipedGainChain==0) haxis6 = gPad->DrawFrame(firstRun-100, 50, lastRun+100, 90);
    else  haxis6 = gPad->DrawFrame(firstRun-100, 400, lastRun+100, 700);
    haxis6->GetXaxis()->SetNoExponent();
    haxis6->SetXTitle("RUN no.");
    haxis6->SetYTitle("PMRef. pedestals");
    //
    graph[ic+22]->SetMarkerStyle(23);
    graph[ic+22]->SetMarkerColor(kGreen+4);
    graph[ic+22]->Draw("P, SAME");
 }
 cothPeds->SaveAs("ZDCPedvsTime2.gif");
 cothPeds->SaveAs("ZDCPedvsTime2.C");
 
}
