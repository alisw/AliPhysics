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

#include "TStatToolkit.h"

TTree *tree;

void DrawPerformanceZDCQAMatchTrends(const char* inFile = "prodQAhistos.root"){

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

  TFile *file0 = TFile::Open(inFile);
  if(!file0) return;
  file0->cd();

  tree = (TTree*) file0->Get("tree");
  if(!tree) return;
  //tree->Print();

  int const entries = tree->GetEntries();
  printf("Total number of analyzed runs=%i\n", entries);

Int_t runNumber=0;
Double_t ZNC_mean=0;
Double_t ZNA_mean=0;
Double_t ZPA_mean=0;
Double_t ZPC_mean=0;
Double_t ZNCuncalib_mean=0;
Double_t ZNAuncalib_mean=0;
Double_t ZPAuncalib_mean=0;
Double_t ZPCuncalib_mean=0;
Double_t ZEM1_mean=0;
Double_t ZEM2_mean=0;
Double_t ZNC_XCent=0;
Double_t ZNC_YCent=0;
Double_t ZNA_XCent=0;
Double_t ZNA_YCent=0;
Double_t ZNC_XCent_err=0;
Double_t ZNC_YCent_err=0;
Double_t ZNA_XCent_err=0;
Double_t ZNA_YCent_err=0;
Double_t ZN_TDC_Sum=0;
Double_t ZN_TDC_Diff=0;
Double_t ZN_TDC_Sum_err=0;
Double_t ZN_TDC_Diff_err=0;

tree->SetBranchAddress("run",&runNumber);
tree->SetBranchAddress("ZNC_mean_value",&ZNC_mean);
tree->SetBranchAddress("ZNA_mean_value",&ZNA_mean);
tree->SetBranchAddress("ZPC_mean_value",&ZPC_mean);
tree->SetBranchAddress("ZPA_mean_value",&ZPA_mean);
tree->SetBranchAddress("ZNC_mean_uncalib",&ZNCuncalib_mean);
tree->SetBranchAddress("ZNA_mean_uncalib",&ZNAuncalib_mean);
tree->SetBranchAddress("ZPC_mean_uncalib",&ZPCuncalib_mean);
tree->SetBranchAddress("ZPA_mean_uncalib",&ZPAuncalib_mean);
tree->SetBranchAddress("ZEM1_mean_value",&ZEM1_mean);
tree->SetBranchAddress("ZEM2_mean_value",&ZEM2_mean);
tree->SetBranchAddress("ZNC_X_Centroid",&ZNC_XCent);
tree->SetBranchAddress("ZNC_Y_Centroid",&ZNC_YCent);
tree->SetBranchAddress("ZNA_X_Centroid",&ZNA_XCent);
tree->SetBranchAddress("ZNA_Y_Centroid",&ZNA_YCent);
tree->SetBranchAddress("ZNC_X_Centroid_Err",&ZNC_XCent_err);
tree->SetBranchAddress("ZNC_Y_Centroid_Err",&ZNC_YCent_err);
tree->SetBranchAddress("ZNA_X_Centroid_Err",&ZNA_XCent_err);
tree->SetBranchAddress("ZNA_Y_Centroid_Err",&ZNA_YCent_err);
tree->SetBranchAddress("ZN_TDC_Sum",&ZN_TDC_Sum);
tree->SetBranchAddress("ZN_TDC_Diff",&ZN_TDC_Diff);
tree->SetBranchAddress("ZN_TDC_Sum_Err",&ZN_TDC_Sum_err);
tree->SetBranchAddress("ZN_TDC_Diff_Err",&ZN_TDC_Diff_err);

  Int_t offset_signals=50;
  Float_t offset_centroids=1.5;
  Int_t offset_tdc=6;

  //create plots of trending variables as a function of the run
  TGraphErrors *gr_znc = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZNC_mean_value:run","",21,kBlue+2,1.2);
  TGraphErrors *gr_zna = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZNA_mean_value:run","",20,kPink,1.2);
  TGraphErrors *gr_zpc = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZPC_mean_value:run","",21,kBlue+2,1.2);
  TGraphErrors *gr_zpa = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZPA_mean_value:run","",20,kPink,1.2);
  TGraphErrors *gr_zncUncal = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZNC_mean_uncalib:run","",21,kAzure+8,1.2);
  TGraphErrors *gr_znaUncal = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZNA_mean_uncalib:run","",20,kPink,1.2);
  TGraphErrors *gr_zpcUncal = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZPC_mean_uncalib:run","",21,kAzure+8,1.2);
  TGraphErrors *gr_zpaUncal = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZPA_mean_uncalib:run","",20,kPink,1.2);
  TGraphErrors *gr_zem1 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZEM1_mean_value:run","",20,kGreen,1.2);
  TGraphErrors *gr_zem2 = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZEM2_mean_value:run","",21,kTeal+1,1.2);
  TGraphErrors *gr_zncXcen = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZNC_X_Centroid:run","",20,kAzure+10,1.2);
  TGraphErrors *gr_zncYcen = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZNC_Y_Centroid:run","",21,kAzure-2,1.2);
  TGraphErrors *gr_znaXcen = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZNA_X_Centroid:run","",20,kPink-2,1.2);
  TGraphErrors *gr_znaYcen = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZNA_Y_Centroid:run","",21,kPink-3,1.2);
  TGraphErrors *gr_zntdcsum = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZN_TDC_Sum:run","",20,kPink,1.2);
  TGraphErrors *gr_zntdcdiff = (TGraphErrors*) TStatToolkit::MakeGraphSparse(tree,"ZN_TDC_Diff:run","",20,kBlue+3,1.2);

  Double_t ZNC_tot=0;    Double_t ZNC_avg=0;
  Double_t ZNA_tot=0;    Double_t ZNA_avg=0;
  Double_t ZPC_tot=0;    Double_t ZPC_avg=0;
  Double_t ZPA_tot=0;    Double_t ZPA_avg=0;
  Double_t ZEM1_tot=0;   Double_t ZEM1_avg=0;
  Double_t ZEM2_tot=0;   Double_t ZEM2_avg=0;
  Double_t ZNC_mean_uncal=0;   Double_t ZNC_uncal_tot=0;    Double_t ZNC_uncal_avg=0;
  Double_t ZNA_mean_uncal=0;   Double_t ZNA_uncal_tot=0;    Double_t ZNA_uncal_avg=0;
  Double_t ZPC_mean_uncal=0;   Double_t ZPC_uncal_tot=0;    Double_t ZPC_uncal_avg=0;
  Double_t ZPA_mean_uncal=0;   Double_t ZPA_uncal_tot=0;    Double_t ZPA_uncal_avg=0;

  Double_t ZNCx_mean=0;  Double_t ZNCx_tot=0;   Double_t ZNCx_avg=0;
  Double_t ZNCy_mean=0;  Double_t ZNCy_tot=0;   Double_t ZNCy_avg=0;
  Double_t ZNAx_mean=0;  Double_t ZNAx_tot=0;   Double_t ZNAx_avg=0;
  Double_t ZNAy_mean=0;  Double_t ZNAy_tot=0;   Double_t ZNAy_avg=0;

  Double_t TdcSum_mean=0;  Double_t TdcSum_tot=0;   Double_t TdcSum_avg=0;
  Double_t TdcDiff_mean=0; Double_t TdcDiff_tot=0;  Double_t TdcDiff_avg=0;

  for (int i=0; i<entries; i++){
     tree->GetEntry(i);
    printf("ZNC %f ZNA %f\n", ZNC_mean, ZNA_mean);
    printf("ZEM1 %f ZEM2 %f\n", ZEM1_mean, ZEM2_mean);
     ZNC_tot += ZNC_mean;
     ZNA_tot += ZNA_mean;
     ZPC_tot += ZPC_mean;
     ZPA_tot += ZPA_mean;
     ZNC_uncal_tot += ZNCuncalib_mean;
     ZNA_uncal_tot += ZNAuncalib_mean;
     ZPC_uncal_tot += ZPCuncalib_mean;
     ZPA_uncal_tot += ZPAuncalib_mean;
     ZEM1_tot += ZEM1_mean;
     ZEM2_tot += ZEM2_mean;
     ZNCx_tot += ZNCx_mean;
     ZNCy_tot += ZNCy_mean;
     ZNAx_tot += ZNAx_mean;
     ZNAy_tot += ZNAy_mean;
     TdcSum_tot += TdcSum_mean;
     TdcDiff_tot += TdcDiff_mean;
  }

   ZNC_avg=ZNC_tot/entries;
   ZNA_avg=ZNA_tot/entries;
   ZPC_avg=ZPC_tot/entries;
   ZPA_avg=ZPA_tot/entries;
   ZNC_uncal_avg=ZNC_uncal_tot/entries;
   ZNA_uncal_avg=ZNA_uncal_tot/entries;
   ZPC_uncal_avg=ZPC_uncal_tot/entries;
   ZPA_uncal_avg=ZPA_uncal_tot/entries;
   ZEM1_avg=ZEM1_tot/entries;
   ZEM2_avg=ZEM2_tot/entries;
   ZNCx_avg=ZNCx_tot/entries;
   ZNCy_avg=ZNCy_tot/entries;
   ZNAx_avg=ZNAx_tot/entries;
   ZNAy_avg=ZNAy_tot/entries;
   TdcSum_avg=TdcSum_tot/entries;
   TdcDiff_avg=TdcDiff_tot/entries;

  TCanvas *c1 = new TCanvas("c1","c1",0,0,1600,1000);
  c1->SetGrid(3); c1->Divide(1,2);
  c1->cd(1);
  gr_znc->Draw("");
  gr_znc->SetTitle("ZNC mean values");
  gr_znc->GetYaxis()->SetTitle("ZNC");
  float min = 0;
  if(ZNC_avg-offset_signals>0) min = ZNC_avg-offset_signals;
  else min = 0.;
  gr_znc->GetYaxis()->SetRangeUser(min, ZNC_avg+offset_signals);
  c1->cd(2);
  gr_zna->Draw("");
  gr_zna->SetTitle("ZNA mean values");
  gr_zna->GetYaxis()->SetTitle("ZNA");
  if(ZNA_avg-offset_signals>0) min = ZNA_avg-offset_signals;
  else min = 0.;
  gr_zna->GetYaxis()->SetRangeUser(min, ZNA_avg+offset_signals);
  c1->Print("ZN_signals_trending.png");

  TCanvas *c2 = new TCanvas("c2","c2",0,0,1600,1000);
  c2->SetGrid(3);
  c2->Divide(1,2);
  c2->cd(1);
  gr_zpc->Draw("");
  gr_zpc->SetTitle("ZPC mean values");
  gr_zpc->GetYaxis()->SetTitle("(ADC ch.)");
  if(ZPC_avg-offset_signals>0) min = ZPC_avg-offset_signals;
  else min = 0.;
  gr_zpc->GetYaxis()->SetRangeUser(min, ZPC_avg+offset_signals);
  c2->cd(2);
  gr_zpa->Draw("");
  gr_zpa->SetTitle("ZPA mean values");
  gr_zpa->GetYaxis()->SetTitle("(ADC ch.)");
  if(ZPA_avg-offset_signals>0) min = ZPA_avg-offset_signals;
  else min = 0.;
  gr_zpa->GetYaxis()->SetRangeUser(min, ZPA_avg+offset_signals);
  c2->Print("ZP_signals_trending.png");

  TCanvas *c1u = new TCanvas("c1u","c1u",0,0,1600,1000);
  c1u->SetGrid(3); c1u->Divide(1,2);
  c1u->cd(1);
  gr_zncUncal->Draw("");
  gr_zncUncal->SetTitle("ZNC uncalibrated average");
  gr_zncUncal->GetYaxis()->SetTitle("ZNC");
  if(ZNC_uncal_avg-offset_signals>0) min = ZNC_uncal_avg-offset_signals;
  else min = 0.;
  gr_zncUncal->GetYaxis()->SetRangeUser(min, ZNC_uncal_avg+offset_signals);
  c1u->cd(2);
  gr_znaUncal->Draw("");
  gr_znaUncal->SetTitle("ZNA uncalibrated average");
  gr_znaUncal->GetYaxis()->SetTitle("ZNA");
  if(ZNA_uncal_avg-offset_signals>0) min = ZNA_uncal_avg-offset_signals;
  else min = 0.;
  gr_znaUncal->GetYaxis()->SetRangeUser(min, ZNA_uncal_avg+offset_signals);
  c1u->Print("ZN_uncalibsignals_trending.png");

  TCanvas *c2u = new TCanvas("c2u","c2u",0,0,1600,1000);
  c2u->SetGrid(3);
  c2u->Divide(1,2);
  c2u->cd(1);
  gr_zpcUncal->Draw("");
  gr_zpcUncal->SetTitle("ZPC uncalibrated average");
  gr_zpcUncal->GetYaxis()->SetTitle("(ADC ch.)");
  if(ZPC_uncal_avg-offset_signals>0) min = ZPC_uncal_avg-offset_signals;
  else min = 0.;
  gr_zpcUncal->GetYaxis()->SetRangeUser(min, ZPC_uncal_avg+offset_signals);
  c2u->cd(2);
  gr_zpaUncal->Draw("");
  gr_zpaUncal->SetTitle("ZPA uncalibrated average");
  gr_zpaUncal->GetYaxis()->SetTitle("(ADC ch.)");
  if(ZPA_uncal_avg-offset_signals>0) min = ZPA_uncal_avg-offset_signals;
  else min = 0.;
  gr_zpaUncal->GetYaxis()->SetRangeUser(min, ZPA_uncal_avg+offset_signals);
  c2u->Print("ZP_uncalibsignals_trending.png");

  TCanvas *c3 = new TCanvas("c3","c3",1600,1000);
  c3->SetGrid(3); c3->Divide(1,2);
  c3->cd(1);
  gr_zem1->Draw("");
  gr_zem1->SetTitle("ZEM1 mean values");
  gr_zem1->GetYaxis()->SetTitle("(ADC ch.)");
  if(ZEM1_avg-offset_signals>0) min = ZEM1_avg-offset_signals;
  else min = 0.;
  gr_zem1->GetYaxis()->SetRangeUser(min, ZEM1_avg+offset_signals);
  c3->cd(2);
  gr_zem2->Draw("");
  gr_zem2->SetTitle("ZEM2 mean values");
  gr_zem2->GetYaxis()->SetTitle("(ADC ch.)");
  if(ZEM2_avg-offset_signals>0) min = ZEM2_avg-offset_signals;
  else min = 0.;
  gr_zem2->GetYaxis()->SetRangeUser(min, ZEM2_avg+offset_signals);
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
