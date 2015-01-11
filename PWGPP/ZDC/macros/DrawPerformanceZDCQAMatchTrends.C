#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <Riostream.h>
#include <TClassTable.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TLine.h>
#include <TGrid.h>
#include <TBits.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFileMerger.h>
#include <TGridResult.h>
#include <TSystem.h>
#include <TGaxis.h>
#include <TRandom.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include "TStatToolkit.h"
#endif

TTree *tree;

void DrawPerformanceZDCQAMatchTrends(const char* inFile = "trending.root"){
  
  /*set graphic style*/
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23); 
  gStyle->SetTextFont(42);
  gStyle->SetStatColor(kWhite); 
  gStyle->SetStatBorderSize(1);
  gStyle->SetOptStat(0);  
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetLabelSize(0.03,"x");
  gStyle->SetLabelSize(0.03,"y");  
  gStyle->SetLabelOffset(0.015,"x");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetTitleOffset(1.4,"x");
  gStyle->SetTitleOffset(1.0,"y");     

  TFile *_file0 = TFile::Open(inFile);
  if(!_file0) return;
  _file0->cd();
  
  tree = (TTree*)_file0->Get("tree");
  if(!tree) return;
  
  int const entries = tree->GetEntries();
  printf("Total number of analyzed runs=%i\n", entries);
  
  Int_t offset_signals=200;
  Int_t offset_centroids=3;
  Int_t offset_tdc=6;  
  
  //Int_t runNumber=0;
  Double_t ZNC_mean=0;   Double_t ZNC_tot=0;    Double_t ZNC_avg=0;  
  Double_t ZNA_mean=0;   Double_t ZNA_tot=0;    Double_t ZNA_avg=0;
  Double_t ZPC_mean=0;   Double_t ZPC_tot=0;    Double_t ZPC_avg=0;     
  Double_t ZPA_mean=0;   Double_t ZPA_tot=0;    Double_t ZPA_avg=0;
  Double_t ZEM1_mean=0;  Double_t ZEM1_tot=0;   Double_t ZEM1_avg=0;    
  Double_t ZEM2_mean=0;  Double_t ZEM2_tot=0;   Double_t ZEM2_avg=0;   
  
  Double_t ZNCx_mean=0;  Double_t ZNCx_tot=0;   Double_t ZNCx_avg=0; 
  Double_t ZNCy_mean=0;  Double_t ZNCy_tot=0;   Double_t ZNCy_avg=0;
  Double_t ZNAx_mean=0;  Double_t ZNAx_tot=0;   Double_t ZNAx_avg=0; 
  Double_t ZNAy_mean=0;  Double_t ZNAy_tot=0;   Double_t ZNAy_avg=0;
  
  Double_t TdcSum_mean=0;  Double_t TdcSum_tot=0;   Double_t TdcSum_avg=0;
  Double_t TdcDiff_mean=0; Double_t TdcDiff_tot=0;  Double_t TdcDiff_avg=0;  

  //create plots of trending variables as a function of the run
  TGraphErrors *gr_znc = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZNC_mean_value:run","",20,kRed,1.0);
  TGraphErrors *gr_zna = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZNA_mean_value:run","",22,kBlue,1.2);  
  TGraphErrors *gr_zpc = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZPC_mean_value:run","",20,kRed,1.0);
  TGraphErrors *gr_zpa = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZPA_mean_value:run","",22,kBlue,1.2);   
  TGraphErrors *gr_zem1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZEM1_mean_value:run","",20,kRed,1.0);
  TGraphErrors *gr_zem2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZEM2_mean_value:run","",22,kBlue,1.2);   
  TGraphErrors *gr_zncXcen = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZNC_X_Centroid:run","",20,kRed,1.0);
  TGraphErrors *gr_zncYcen = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZNC_Y_Centroid:run","",22,kBlue,1.2);     
  TGraphErrors *gr_znaXcen = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZNA_X_Centroid:run","",20,kRed,1.0);
  TGraphErrors *gr_znaYcen = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZNA_Y_Centroid:run","",22,kBlue,1.2);    
  TGraphErrors *gr_zntdcsum = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZN_TDC_Sum:run","",20,kRed,1.0);
  TGraphErrors *gr_zntdcdiff = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZN_TDC_Diff:run","",22,kBlue,1.2);  
    
  for (int i=0;i<entries;i++){
   tree->GetEntry(i);
   ZNC_tot+=ZNC_mean;
   ZNA_tot+=ZNA_mean;
   ZPC_tot+=ZPC_mean;
   ZPA_tot+=ZPA_mean;
   ZEM1_tot+=ZEM1_mean;   
   ZEM2_tot+=ZEM2_mean;
   ZNCx_tot+=ZNCx_mean;
   ZNCy_tot+=ZNCy_mean;
   ZNAx_tot+=ZNAx_mean;
   ZNAy_tot+=ZNAy_mean;
   TdcSum_tot+=TdcSum_mean;
   TdcDiff_tot+=TdcDiff_mean;  
  }
  
   ZNC_avg=ZNC_tot/entries;
   ZNA_avg=ZNA_tot/entries; 
   ZPC_avg=ZPC_tot/entries; 
   ZPA_avg=ZPA_tot/entries;  
   ZEM1_avg=ZEM1_tot/entries;   
   ZEM2_avg=ZEM2_tot/entries;
   ZNCx_avg=ZNCx_tot/entries; 
   ZNCy_avg=ZNCy_tot/entries; 
   ZNAx_avg=ZNAx_tot/entries; 
   ZNAy_avg=ZNAy_tot/entries; 
   TdcSum_avg=TdcSum_tot/entries;
   TdcDiff_avg=TdcDiff_tot/entries;
  
  TCanvas *c1 = new TCanvas("c1","c1",1600,1000);
  c1->SetGrid(3); c1->Divide(1,2);
  c1->cd(1); 
  gr_znc->Draw("");
  gr_znc->SetTitle("ZNC mean values");
  gr_znc->GetYaxis()->SetTitle("(ADC ch.)");
  gr_znc->GetYaxis()->SetRangeUser(ZNC_avg-offset_signals,ZNC_avg+offset_signals);
  c1->cd(2); 
  gr_zna->Draw("");
  gr_zna->SetTitle("ZNA mean values");
  gr_zna->GetYaxis()->SetTitle("(ADC ch.)");
  gr_zna->GetYaxis()->SetRangeUser(ZNA_avg-offset_signals,ZNA_avg+offset_signals);  
  c1->Print("ZN_signals_trending.png");
  
  TCanvas *c2 = new TCanvas("c2","c2",1600,1000);
  c2->SetGrid(3); 
  c2->Divide(1,2);
  c2->cd(1); 
  gr_zpc->Draw("");
  gr_zpc->SetTitle("ZPC mean values");
  gr_zpc->GetYaxis()->SetTitle("(ADC ch.)"); 
  gr_zpc->GetYaxis()->SetRangeUser(ZPC_avg-offset_signals,ZPC_avg+offset_signals);  
  c2->cd(2); 
  gr_zpa->Draw("");
  gr_zpa->SetTitle("ZPA mean values");
  gr_zpa->GetYaxis()->SetTitle("(ADC ch.)");
  gr_zpa->GetYaxis()->SetRangeUser(ZPA_avg-offset_signals,ZPA_avg+offset_signals);    
  c2->Print("ZP_signals_trending.png");
  
  TCanvas *c3 = new TCanvas("c3","c3",1600,1000);
  c3->SetGrid(3); c3->Divide(1,2);
  c3->cd(1); 
  gr_zem1->Draw("");
  gr_zem1->SetTitle("ZEM1 mean values");  
  gr_zem1->GetYaxis()->SetTitle("(ADC ch.)");
  gr_zem1->GetYaxis()->SetRangeUser(ZEM1_avg-offset_signals,ZEM1_avg+offset_signals);   
  c3->cd(2); 
  gr_zem2->Draw("");
  gr_zem2->SetTitle("ZEM2 mean values");  
  gr_zem2->GetYaxis()->SetTitle("(ADC ch.)");
  gr_zem2->GetYaxis()->SetRangeUser(ZEM2_avg-offset_signals,ZEM2_avg+offset_signals);    
  c3->Print("ZEM_signals_trending.png");    
  
  TCanvas *c4 = new TCanvas("c4","c4",1600,1000);
  c4->SetGrid(3); c4->Divide(1,2);
  c4->cd(1); 
  gr_zncXcen->Draw("");
  gr_zncXcen->SetTitle("ZNC X coordinate");    
  gr_zncXcen->GetYaxis()->SetTitle("(cm)");
  gr_zncXcen->GetYaxis()->SetRangeUser(ZNCx_avg-offset_centroids,ZNCx_avg+offset_centroids);   
  c4->cd(2); 
  gr_zncYcen->Draw("");
  gr_zncYcen->SetTitle("ZNC Y coordinate");  
  gr_zncYcen->GetYaxis()->SetTitle("(cm)");
  gr_zncYcen->GetYaxis()->SetRangeUser(ZNCy_avg-offset_centroids,ZNCy_avg+offset_centroids);     
  c4->Print("ZNC_centroids_trending.png");   
  
  TCanvas *c5 = new TCanvas("c5","c5",1600,1000);
  c5->SetGrid(3); 
  c5->Divide(1,2);
  c5->cd(1); 
  gr_znaXcen->Draw("");
  gr_znaXcen->SetTitle("ZNA X coordinate");  
  gr_znaXcen->GetYaxis()->SetTitle("(cm)");
  gr_znaXcen->GetYaxis()->SetRangeUser(ZNAx_avg-offset_centroids,ZNAx_avg+offset_centroids);   
  c5->cd(2); 
  gr_znaYcen->Draw("");
  gr_znaYcen->SetTitle("ZNA Y coordinate");    
  gr_znaYcen->GetYaxis()->SetTitle("(cm)");
  gr_znaYcen->GetYaxis()->SetRangeUser(ZNAy_avg-offset_centroids,ZNAy_avg+offset_centroids);    
  c5->Print("ZNA_centroids_trending.png");    
  
  TCanvas *c6 = new TCanvas("c6","c6",1600,1000);
  c6->SetGrid(3); c6->Divide(1,2);
  c6->cd(1);
  gr_zntdcsum->Draw("");
  gr_zntdcsum->SetTitle("ZNC tdc + ZNA tdc"); 
  gr_zntdcsum->GetYaxis()->SetTitle("(ns)");
  gr_zntdcsum->GetYaxis()->SetRangeUser(TdcSum_avg-offset_tdc,TdcSum_avg+offset_tdc);   
  c6->cd(2); 
  gr_zntdcdiff->Draw("");
  gr_zntdcdiff->SetTitle("ZNC tdc - ZNA tdc");   
  gr_zntdcdiff->GetYaxis()->SetTitle("(ns)"); 
  gr_zntdcdiff->GetYaxis()->SetRangeUser(TdcDiff_avg-offset_tdc,TdcDiff_avg+offset_tdc);  
  c6->Print("ZN_timing_trending.png");    
    
}
