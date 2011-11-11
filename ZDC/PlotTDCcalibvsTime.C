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
#include <TH1F.h>
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
#include "AliZDCTDCCalib.h"

#endif


void PlotTDCcalibvsTime(Int_t year=2011, Int_t firstRun=166000, 
		   Int_t lastRun=167000)
{


  TGrid::Connect("alien:",0,0,"t");
  gSystem->Exec(Form("gbbox find \"/alice/data/%d/OCDB/ZDC/Calib/TDCCalib/\" \"Run*.root\" > calibAlienFiles.txt",year));
  FILE* listruns=fopen("calibAlienFiles.txt","r");
  
  const int kNchannels=6;
  TGraphErrors* graph[6];
  for(Int_t i=0; i<kNchannels; i++){
     graph[i] = new TGraphErrors(0);
     char name[50], title[50];
     sprintf(name,"graph%d",i); sprintf(title,"TDC calib. coeff. %d vs. run#",i);
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
    sscanf(filnam,"/alice/data/%d/OCDB/ZDC/Calib/TDCCalib/Run%d_%d_v%d_s%d.root",&year,&nrun,&nrun2,&nv,&ns);
    if(nrun<firstRun) continue;
    if(nrun>lastRun) continue;
    sprintf(filnamalien,"alien://%s",filnam);
    printf("Opening file: %s\n",filnam);
    TFile *f = TFile::Open(filnamalien);  
    AliCDBEntry *entry = (AliCDBEntry*)f->Get("AliCDBEntry");
    AliZDCTDCCalib *calibdata = dynamic_cast<AliZDCTDCCalib*>  (entry->GetObject());
    
    for(int i=0; i<kNchannels; i++){
         graph[i]->SetPoint(iPoint, (Double_t)nrun, calibdata->GetMeanTDC(i));
        graph[i]->SetPointError(iPoint, 0., calibdata->GetWidthTDC(i));
    }
    iPoint++;
    f->Close();
 }

 TFile *outfile=new TFile(Form("TDCCalib%dVsTime.root",year),"recreate");
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

 TCanvas *cHadPeds = new TCanvas("cHadPeds","Hadronic ZDC pedestals",0,0,1200,800);
 cHadPeds->Divide(3,2);
 for(int ic=0; ic<6; ic++){
   // *** ZNC pedestals
   cHadPeds->cd(ic+1);
   //
   TH1F *haxis1 = gPad->DrawFrame(firstRun-100, -100, lastRun+100, 0);
   haxis1->GetXaxis()->SetNoExponent();
   haxis1->SetXTitle("RUN no.");
   if(ic==0) haxis1->SetYTitle("ZNC TDC calib");
   else if(ic==1)  haxis1->SetYTitle("ZNA TDC calib");
   else if(ic==2)  haxis1->SetYTitle("ZPC TDC calib");
   else if(ic==3)  haxis1->SetYTitle("ZPA TDC calib");
   else if(ic==4)  haxis1->SetYTitle("ZEM1 TDC calib");
   else if(ic==5)  haxis1->SetYTitle("ZEM2 TDC calib");
   //
   graph[ic]->SetMarkerStyle(20);
   graph[ic]->SetMarkerColor(kAzure+ic);
   graph[ic]->Draw("P, SAME");
   // 
 }
 cHadPeds->SaveAs("ZDCTDCvsTime1.gif");
 cHadPeds->SaveAs("ZDCTDCvsTime1.C");

 
}
