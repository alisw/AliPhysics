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
#endif

Int_t DrawPerformanceZDCQAMatch(const char* inFile){

// Draw control histograms and generate output pictures

gSystem->Load("libSTAT");
gSystem->Load("libANALYSIS");
gSystem->Load("libANALYSISalice");
gSystem->Load("libANALYSIScalib");
gSystem->Load("libCORRFW");
gSystem->Load("libANALYSISalice");
gSystem->Load("libANALYSIScalib");
gSystem->Load("libTender");
gSystem->Load("libPWGPP");

gROOT->Reset();
gROOT->SetStyle("Plain");
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetTitleSize(0.025);
TH1::AddDirectory(kFALSE);

TFile * fin = TFile::Open(inFile);
TTree * ttree = (TTree*) fin->Get("trending");  //in
TTree * tree = new TTree("tree","tree");        //out (to be summed)

if (!ttree){
  Printf("Invalid trending tree!!!!!!!!!!!!!!");
  return 2;
}

Int_t nRuns=ttree->GetEntries();
TList list;

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
gStyle->SetLabelSize(0.02,"xyz");
gStyle->SetLabelOffset(0.03,"xyz");

TString plotDir(".");

legend = new TLegend(0.9,0.1,1.0,0.9);
legend->SetFillColor(kWhite);

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
Double_t ZNC_TDC=0;
Double_t ZNA_TDC=0;

TH1F *hZNCpmcUncalib = new TH1F("hZNCpmcUncalib","hZNCpmcUncalib",200.,0.,2000.);
TH1F *hZNApmcUncalib = new TH1F("hZNApmcUncalib","hZNApmcUncalib",200.,0.,2000.);
TH1F *hZPCpmcUncalib = new TH1F("hZPCpmcUncalib","hZPCpmcUncalib",200.,0.,2000.);
TH1F *hZPApmcUncalib = new TH1F("hZPApmcUncalib","hZPApmcUncalib",200.,0.,2000.);
TH1F *hZEM1 = new TH1F("hZEM1","hZEM1",200.,0.,2000.);
TH1F *hZEM2 = new TH1F("hZEM2","hZEM2",200.,0.,2000.);

ttree->SetBranchAddress("run",&runNumber);
ttree->SetBranchAddress("ZNC_mean_value",&ZNC_mean);
ttree->SetBranchAddress("ZNA_mean_value",&ZNA_mean);
ttree->SetBranchAddress("ZPC_mean_value",&ZPC_mean);
ttree->SetBranchAddress("ZPA_mean_value",&ZPA_mean);
ttree->SetBranchAddress("ZNCuncalib_mean",&ZNCuncalib_mean);
ttree->SetBranchAddress("ZNAuncalib_mean",&ZNAuncalib_mean);
ttree->SetBranchAddress("ZPCuncalib_mean",&ZPCuncalib_mean);
ttree->SetBranchAddress("ZPAuncalib_mean",&ZPAuncalib_mean);
ttree->SetBranchAddress("ZEM1_mean_value",&ZEM1_mean);
ttree->SetBranchAddress("ZEM2_mean_value",&ZEM2_mean);
ttree->SetBranchAddress("ZNC_X_Centroid",&ZNC_XCent);
ttree->SetBranchAddress("ZNC_Y_Centroid",&ZNC_YCent);
ttree->SetBranchAddress("ZNA_X_Centroid",&ZNA_XCent);
ttree->SetBranchAddress("ZNA_Y_Centroid",&ZNA_YCent);
ttree->SetBranchAddress("ZNC_X_Centroid_Err",&ZNC_XCent_err);
ttree->SetBranchAddress("ZNC_Y_Centroid_Err",&ZNC_YCent_err);
ttree->SetBranchAddress("ZNA_X_Centroid_Err",&ZNA_XCent_err);
ttree->SetBranchAddress("ZNA_Y_Centroid_Err",&ZNA_YCent_err);
ttree->SetBranchAddress("ZN_TDC_Sum",&ZN_TDC_Sum);
ttree->SetBranchAddress("ZN_TDC_Diff",&ZN_TDC_Diff);
ttree->SetBranchAddress("ZN_TDC_Sum_Err",&ZN_TDC_Sum_err);
ttree->SetBranchAddress("ZN_TDC_Diff_Err",&ZN_TDC_Diff_err);
//ttree->SetBranchAddress("ZNC_TDC",&ZNC_TDC);
//ttree->SetBranchAddress("ZNA_TDC",&ZNA_TDC);

TH1F *hznc = new TH1F("hznc","ZNC average signal",3,-1,1);
hznc->GetXaxis()->SetRangeUser(-1.,1.);
hznc->SetDrawOption("EP");
hznc->SetMarkerStyle(20);
hznc->SetMarkerColor(kRed);
hznc->SetLineColor(kRed);

TH1F *hzna = new TH1F("hzna","ZNA average signal",3,-1,1);
hzna->GetXaxis()->SetRangeUser(-1.,1.);
hzna->SetDrawOption("EP");
hzna->SetMarkerStyle(20);
hzna->SetMarkerColor(kRed);
hzna->SetLineColor(kRed);

TH1F *hzpc = new TH1F("hzpc","ZPC average signal",3,-1,1);
hzpc->GetXaxis()->SetRangeUser(-1.,1.);
hzpc->SetDrawOption("EP");
hzpc->SetMarkerStyle(20);
hzpc->SetMarkerColor(kRed);
hzpc->SetLineColor(kRed);

TH1F *hzpa = new TH1F("hzpa","ZPA average signal",3,-1,1);
hzpa->GetXaxis()->SetRangeUser(-1.,1.);
hzpa->SetDrawOption("EP");
hzpa->SetMarkerStyle(20);
hzpa->SetMarkerColor(kRed);
hzpa->SetLineColor(kRed);

TH1F *hzncUncalib = new TH1F("hzncUncalib","ZNC uncalibrated average signal",3,-1,1);
hzncUncalib->GetXaxis()->SetRangeUser(-1.,1.);
hzncUncalib->SetDrawOption("EP");
hzncUncalib->SetMarkerStyle(20);
hzncUncalib->SetMarkerColor(kAzure+10);
hzncUncalib->SetLineColor(kAzure+10);

TH1F *hznaUncalib = new TH1F("hznaUncalib","ZNA uncalibrated average signal",3,-1,1);
hznaUncalib->GetXaxis()->SetRangeUser(-1.,1.);
hznaUncalib->SetDrawOption("EP");
hznaUncalib->SetMarkerStyle(20);
hznaUncalib->SetMarkerColor(kAzure+10);
hznaUncalib->SetLineColor(kAzure+10);

TH1F *hzpcUncalib = new TH1F("hzpcUncalib","ZPC uncalibrated average signal",3,-1,1);
hzpcUncalib->GetXaxis()->SetRangeUser(-1.,1.);
hzpcUncalib->SetDrawOption("EP");
hzpcUncalib->SetMarkerStyle(20);
hzpcUncalib->SetMarkerColor(kAzure+10);
hzpcUncalib->SetLineColor(kAzure+10);

TH1F *hzpaUncalib = new TH1F("hzpaUncalib","ZPA uncalibrated average signal",3,-1,1);
hzpaUncalib->GetXaxis()->SetRangeUser(-1.,1.);
hzpaUncalib->SetDrawOption("EP");
hzpaUncalib->SetMarkerStyle(20);
hzpaUncalib->SetMarkerColor(kAzure+10);
hzpaUncalib->SetLineColor(kAzure+10);

TH1F *hzem1 = new TH1F("hzem1","ZEM1 average signal",3,-1,1);
hzem1->GetXaxis()->SetRangeUser(-1.,1.);
hzem1->SetDrawOption("EP");
hzem1->SetMarkerStyle(20);
hzem1->SetMarkerColor(kRed);
hzem1->SetLineColor(kRed);

TH1F *hzem2 = new TH1F("hzem2","ZEM2 average signal",3,-1,1);
hzem2->GetXaxis()->SetRangeUser(-1.,1.);
hzem2->SetDrawOption("EP");
hzem2->SetMarkerStyle(20);
hzem2->SetMarkerColor(kRed);
hzem2->SetLineColor(kRed);

TH1F *hzna_Xcentroid = new TH1F("hzna_Xcentroid","ZNA X centroid",3,-1,1);
hzna_Xcentroid->GetXaxis()->SetRangeUser(-1.,1.);
hzna_Xcentroid->SetDrawOption("EP");
hzna_Xcentroid->SetMarkerStyle(20);
hzna_Xcentroid->SetMarkerColor(kRed);
hzna_Xcentroid->SetLineColor(kRed);

TH1F *hzna_Ycentroid = new TH1F("hzna_Ycentroid","ZNA Y centroid",3,-1,1);
hzna_Ycentroid->GetXaxis()->SetRangeUser(-1.,1.);
hzna_Ycentroid->SetDrawOption("EP");
hzna_Ycentroid->SetMarkerStyle(20);
hzna_Ycentroid->SetMarkerColor(kRed);
hzna_Ycentroid->SetLineColor(kRed);

TH1F *hznc_Xcentroid = new TH1F("hznc_Xcentroid","ZNC X centroid",3,-1,1);
hznc_Xcentroid->GetXaxis()->SetRangeUser(-1.,1.);
hznc_Xcentroid->SetDrawOption("EP");
hznc_Xcentroid->SetMarkerStyle(20);
hznc_Xcentroid->SetMarkerColor(kRed);
hznc_Xcentroid->SetLineColor(kRed);

TH1F *hznc_Ycentroid = new TH1F("hznc_Ycentroid","ZNC Y centroid",3,-1,1);
hznc_Ycentroid->GetXaxis()->SetRangeUser(-1.,1.);
hznc_Ycentroid->SetDrawOption("EP");
hznc_Ycentroid->SetMarkerStyle(20);
hznc_Ycentroid->SetMarkerColor(kRed);
hznc_Ycentroid->SetLineColor(kRed);

TH1F *hzn_TDC_Sum = new TH1F("hzn_TDC_Sum","ZNC TDC + ZNA TDC",3,-1,1);
hzn_TDC_Sum->GetXaxis()->SetRangeUser(-1.,1.);
hzn_TDC_Sum->SetDrawOption("EP");
hzn_TDC_Sum->SetMarkerStyle(20);
hzn_TDC_Sum->SetMarkerColor(kRed);
hzn_TDC_Sum->SetLineColor(kRed);

TH1F *hzn_TDC_Diff = new TH1F("hzn_TDC_Diff","ZNC TDC - ZNA TDC",3,-1,1);
hzn_TDC_Diff->GetXaxis()->SetRangeUser(-1.,1.);
hzn_TDC_Diff->SetDrawOption("EP");
hzn_TDC_Diff->SetMarkerStyle(20);
hzn_TDC_Diff->SetMarkerColor(kRed);
hzn_TDC_Diff->SetLineColor(kRed);

TH1F *hznc_TDC = new TH1F("hznc_TDC","ZNC TDC",3,-1,1);
hznc_TDC->GetXaxis()->SetRangeUser(-1.,1.);
hznc_TDC->SetDrawOption("EP");
hznc_TDC->SetMarkerStyle(20);
hznc_TDC->SetMarkerColor(kRed);
hznc_TDC->SetLineColor(kRed);

TH1F *hzna_TDC = new TH1F("hzna_TDC","ZNA TDC",3,-1,1);
hzna_TDC->GetXaxis()->SetRangeUser(-1.,1.);
hzna_TDC->SetDrawOption("EP");
hzna_TDC->SetMarkerStyle(20);
hzna_TDC->SetMarkerColor(kRed);
hzna_TDC->SetLineColor(kRed);

char runlabel[10];

for (Int_t irun=0;irun<nRuns;irun++){
  ttree->GetEntry(irun);
 }

sprintf(runlabel,"%i",runNumber);

//----------------------------------------------------------------------
//spectrum vs run
//----------------------------------------------------------------------

hZNCpmcUncalib = dynamic_cast<TH1F*> (fin->Get("fhZNCpmcUncalib"));
//hZNCpmc->GetXaxis()->SetRangeUser(0.,2000);
if(hZNCpmcUncalib->GetEntries()>0. ) hZNCpmcUncalib->Scale(1./hZNCpmcUncalib->GetEntries());
hZNCpmcUncalib->SetLineColor(kRed);
hZNCpmcUncalib->SetLineWidth(2);
hZNCpmcUncalib->SetTitle("ZNC spectrum");
hZNCpmcUncalib->SetXTitle("ZNC signal ");

hZNApmcUncalib = dynamic_cast<TH1F*> (fin->Get("fhZNApmcUncalib"));
//hZNApmc->GetXaxis()->SetRangeUser(0.,2000);
if(hZNApmcUncalib->GetEntries()>0. ) hZNApmcUncalib->Scale(1./hZNApmcUncalib->GetEntries());
hZNApmcUncalib->SetLineColor(kRed);
hZNApmcUncalib->SetLineWidth(2);
hZNApmcUncalib->SetTitle("ZNA spectrum");
hZNApmcUncalib->SetXTitle("ZNA signal ");

hZPCpmcUncalib = dynamic_cast<TH1F*> (fin->Get("fhZPCpmcUncalib"));
//hZPCpmc->GetXaxis()->SetRangeUser(0.,2000);
if(hZPCpmcUncalib->GetEntries()>0. ) hZPCpmcUncalib->Scale(1./hZPCpmcUncalib->GetEntries());
hZPCpmcUncalib->SetLineColor(kRed);
hZPCpmcUncalib->SetLineWidth(2);
hZPCpmcUncalib->SetTitle("ZPC spectrum");
hZPCpmcUncalib->SetXTitle("ZPC signal ");

hZPApmcUncalib = dynamic_cast<TH1F*> (fin->Get("fhZPApmcUncalib"));
//hZPApmc->GetXaxis()->SetRangeUser(0.,2000);
if(hZPApmcUncalib->GetEntries()>0. ) hZPApmcUncalib->Scale(1./hZPApmcUncalib->GetEntries());
hZPApmcUncalib->SetLineColor(kRed);
hZPApmcUncalib->SetLineWidth(2);
hZPApmcUncalib->SetTitle("ZPA spectrum");
hZPApmcUncalib->SetXTitle("ZPA signal ");

hZEM1 = dynamic_cast<TH1F*> (fin->Get("fhZEM1Spectrum"));
//hZEM1->GetXaxis()->SetRangeUser(0.,2000);
if(hZEM1->GetEntries()>0.) hZEM1->Scale(1./hZEM1->GetEntries());
hZEM1->SetLineColor(kRed);
hZEM1->SetLineWidth(2);
hZEM1->SetTitle("ZEM1 spectrum");
hZEM1->SetXTitle("ZEM1 signal (ADC ch.)");

hZEM2 = dynamic_cast<TH1F*> (fin->Get("fhZEM2Spectrum"));
//hZEM2->GetXaxis()->SetRangeUser(0.,2000);
if(hZEM2->GetEntries()>0.) hZEM2->Scale(1./hZEM2->GetEntries());
hZEM2->SetLineColor(kRed);
hZEM2->SetLineWidth(2);
hZEM2->SetTitle("ZEM2 spectrum");
hZEM2->SetXTitle("ZEM2 signal (ADC ch.)");

//----------------------------------------------------------------------
//variables vs run
//----------------------------------------------------------------------
hzna->SetBinContent(2,ZNA_mean);
hzna->GetXaxis()->SetBinLabel(2,runlabel);
hzna->GetXaxis()->SetLabelSize(0.05);

hzpc->SetBinContent(2,ZPC_mean);
hzpc->GetXaxis()->SetBinLabel(2,runlabel);
hzpc->GetXaxis()->SetLabelSize(0.05);

hznc->SetBinContent(2,ZNC_mean);
hznc->GetXaxis()->SetBinLabel(2,runlabel);
hznc->GetXaxis()->SetLabelSize(0.05);

hzpa->SetBinContent(2,ZPA_mean);
hzpa->GetXaxis()->SetBinLabel(2,runlabel);
hzpa->GetXaxis()->SetLabelSize(0.05);

hznaUncalib->SetBinContent(2,ZNAuncalib_mean);
hznaUncalib->GetXaxis()->SetBinLabel(2,runlabel);
hznaUncalib->GetXaxis()->SetLabelSize(0.05);

hzpcUncalib->SetBinContent(2,ZPCuncalib_mean);
hzpcUncalib->GetXaxis()->SetBinLabel(2,runlabel);
hzpcUncalib->GetXaxis()->SetLabelSize(0.05);

hzncUncalib->SetBinContent(2,ZNCuncalib_mean);
hzncUncalib->GetXaxis()->SetBinLabel(2,runlabel);
hzncUncalib->GetXaxis()->SetLabelSize(0.05);

hzpaUncalib->SetBinContent(2,ZPAuncalib_mean);
hzpaUncalib->GetXaxis()->SetBinLabel(2,runlabel);
hzpaUncalib->GetXaxis()->SetLabelSize(0.05);

hzem1->SetBinContent(2,ZEM1_mean);
hzem1->GetXaxis()->SetBinLabel(2,runlabel);
hzem1->GetXaxis()->SetLabelSize(0.05);

hzem2->SetBinContent(2,ZEM2_mean);
hzem2->GetXaxis()->SetBinLabel(2,runlabel);
hzem2->GetXaxis()->SetLabelSize(0.05);

hzna_Xcentroid->SetBinContent(2,ZNA_XCent);
hzna_Xcentroid->SetBinError(2,ZNA_XCent_err);
hzna_Xcentroid->GetXaxis()->SetBinLabel(2,runlabel);
hzna_Xcentroid->GetXaxis()->SetLabelSize(0.05);
hzna_Xcentroid->GetYaxis()->SetTitle("(cm)");

hzna_Ycentroid->SetBinContent(2,ZNA_YCent);
hzna_Ycentroid->SetBinError(2,ZNA_YCent_err);
hzna_Ycentroid->GetXaxis()->SetBinLabel(2,runlabel);
hzna_Ycentroid->GetXaxis()->SetLabelSize(0.05);
hzna_Ycentroid->GetYaxis()->SetTitle("(cm)");

hznc_Xcentroid->SetBinContent(2,ZNC_XCent);
hznc_Xcentroid->SetBinError(2,ZNC_XCent_err);
hznc_Xcentroid->GetXaxis()->SetBinLabel(2,runlabel);
hznc_Xcentroid->GetXaxis()->SetLabelSize(0.05);
hznc_Xcentroid->GetYaxis()->SetTitle("(cm)");

hznc_Ycentroid->SetBinContent(2,ZNC_YCent);
hznc_Ycentroid->SetBinError(2,ZNC_YCent_err);
hznc_Ycentroid->GetXaxis()->SetBinLabel(2,runlabel);
hznc_Ycentroid->GetXaxis()->SetLabelSize(0.05);
hznc_Ycentroid->GetYaxis()->SetTitle("(cm)");

hzn_TDC_Sum->SetBinContent(2,ZN_TDC_Sum);
hzn_TDC_Sum->SetBinError(2,ZN_TDC_Sum_err);
hzn_TDC_Sum->GetXaxis()->SetBinLabel(2,runlabel);
hzn_TDC_Sum->GetXaxis()->SetLabelSize(0.05);
hzn_TDC_Sum->GetYaxis()->SetTitle("(ns)");

hzn_TDC_Diff->SetBinContent(2,ZN_TDC_Diff);
hzn_TDC_Diff->SetBinError(2,ZN_TDC_Diff_err);
hzn_TDC_Diff->GetXaxis()->SetBinLabel(2,runlabel);
hzn_TDC_Diff->GetXaxis()->SetLabelSize(0.05);
hzn_TDC_Diff->GetYaxis()->SetTitle("(ns)");

hznc_TDC->SetBinContent(2,ZNC_TDC);
//hznc_TDC->SetBinError(2,ZNC_TDC_err);
hznc_TDC->GetXaxis()->SetBinLabel(2,runlabel);
hznc_TDC->GetXaxis()->SetLabelSize(0.05);
hznc_TDC->GetYaxis()->SetTitle("(ns)");

hzna_TDC->SetBinContent(2,ZNA_TDC);
//hzna_TDC->SetBinError(2,ZNA_TDC_err);
hzna_TDC->GetXaxis()->SetBinLabel(2,runlabel);
hzna_TDC->GetXaxis()->SetLabelSize(0.05);
hzna_TDC->GetYaxis()->SetTitle("(ns)");

//----------------------------------------------------------------------
//spectra
//----------------------------------------------------------------------

TCanvas* cZNC_Spectra_Uncal = new TCanvas("cZNC_Spectra_Uncal","cZNC_Spectra_Uncal",0,0,1200,900);
gPad->SetLogy();
hZNCpmcUncalib->Draw();
cZNC_Spectra_Uncal->Print(Form("%s/cZNC_Spectra_Uncal.png",plotDir.Data()));

TCanvas* cZNA_Spectra_Uncal = new TCanvas("cZNA_Spectra_Uncal","cZNA_Spectra_Uncal",0,0,1200,900);
gPad->SetLogy();
hZNApmcUncalib->Draw();
cZNA_Spectra_Uncal->Print(Form("%s/cZNA_Spectra_Uncal.png",plotDir.Data()));

TCanvas* cZPC_Spectra_Uncal = new TCanvas("cZPC_Spectra_Uncal","cZPC_Spectra_Uncal",0,0,1200,900);
gPad->SetLogy();
hZPCpmcUncalib->Draw();
cZPC_Spectra_Uncal->Print(Form("%s/cZPC_Spectra_Uncal.png",plotDir.Data()));

TCanvas* cZPA_Spectra_Uncal = new TCanvas("cZPA_Spectra_Uncal","cZPA_Spectra_Uncal",0,0,1200,900);
gPad->SetLogy();
hZPApmcUncalib->Draw();
cZPA_Spectra_Uncal->Print(Form("%s/cZPA_Spectra_Uncal.png",plotDir.Data()));

TCanvas* cZEM1_Spectra = new TCanvas("cZEM1_Spectra","cZEM1_Spectra",0,0,1200,900);
gPad->SetLogy();
hZEM1->Draw();
cZEM1_Spectra->Print(Form("%s/cZEM1_Spectra.png",plotDir.Data()));

TCanvas* cZEM2_Spectra = new TCanvas("cZEM2_Spectra","cZEM2_Spectra",0,0,1200,900);
gPad->SetLogy();
hZEM2->Draw();
cZEM2_Spectra->Print(Form("%s/cZEM2_Spectra.png",plotDir.Data()));


//---------------------------------------------------------------------------------------------------
//variables
//---------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------
//means
//---------------------------------------------------------------------------------------------------
TCanvas* cZNC_Mean_Values = new TCanvas("cZNC_Mean_Values","cZNC_Mean_Values", 0,0,750,900);
hznc->Draw("ep");
cZNC_Mean_Values->Print(Form("%s/cZNC_Mean_Values.png",plotDir.Data()));

TCanvas* cZNA_Mean_Values = new TCanvas("cZNA_Mean_Values","cZNA_Mean_Values", 0,0,750,900);
hzna->Draw("ep");
cZNA_Mean_Values->Print(Form("%s/cZNA_Mean_Values.png",plotDir.Data()));

TCanvas* cZPC_Mean_Values = new TCanvas("cZPC_Mean_Values","cZPC_Mean_Values", 0,0,750,900);
hzpc->Draw("ep");
cZPC_Mean_Values->Print(Form("%s/cZPC_Mean_Values.png",plotDir.Data()));

TCanvas* cZPA_Mean_Values = new TCanvas("cZPA_Mean_Values","cZPA_Mean_Values", 0,0,750,900);
hzpa->Draw("ep");
cZPA_Mean_Values->Print(Form("%s/cZPA_Mean_Values.png",plotDir.Data()));

TCanvas* cZNC_Mean_Uncalib = new TCanvas("cZNC_Mean_Uncalib","cZNC_Mean_Uncalib", 0,0,750,900);
hzncUncalib->Draw("ep");
cZNC_Mean_Uncalib->Print(Form("%s/cZNC_Mean_Uncalib.png",plotDir.Data()));

TCanvas* cZNA_Mean_Uncalib = new TCanvas("cZNA_Mean_Uncalib","cZNA_Mean_Uncalib", 0,0,750,900);
hznaUncalib->Draw("ep");
cZNA_Mean_Uncalib->Print(Form("%s/cZNA_Mean_Uncalib.png",plotDir.Data()));

TCanvas* cZPC_Mean_Uncalib = new TCanvas("cZPC_Mean_Uncalib","cZPC_Mean_Uncalib", 0,0,750,900);
hzpcUncalib->Draw("ep");
cZPC_Mean_Uncalib->Print(Form("%s/cZPC_Mean_Uncalib.png",plotDir.Data()));

TCanvas* cZPA_Mean_Uncalib = new TCanvas("cZPA_Mean_Uncalib","cZPA_Mean_Uncalib", 0,0,750,900);
hzpaUncalib->Draw("ep");
cZPA_Mean_Uncalib->Print(Form("%s/cZPA_Mean_Uncalib.png",plotDir.Data()));

TCanvas* cZEM1_Mean_Values = new TCanvas("cZEM1_Mean_Values","cZEM1_Mean_Values", 0,0,750,900);
hzem1->Draw("ep");
cZEM1_Mean_Values->Print(Form("%s/cZEM1_Mean_Values.png",plotDir.Data()));

TCanvas* cZEM2_Mean_Values = new TCanvas("cZEM2_Mean_Values","cZEM2_Mean_Values", 0,0,750,900);
hzem2->Draw("ep");
cZEM2_Mean_Values->Print(Form("%s/cZEM2_Mean_Values.png",plotDir.Data()));

//---------------------------------------------------------------------------------------------------
//centroids
//---------------------------------------------------------------------------------------------------
TCanvas* cZNA_X_centroid = new TCanvas("cZNA_X_centroid","cZNA_X_centroid", 0,0,750,900);
hzna_Xcentroid->Draw();
cZNA_X_centroid->Print(Form("%s/cZNA_X_centroid.png",plotDir.Data()));

TCanvas* cZNA_Y_centroid = new TCanvas("cZNA_Y_centroid","cZNA_Y_centroid", 0,0,750,900);
hzna_Ycentroid->Draw();
cZNA_Y_centroid->Print(Form("%s/cZNA_Y_centroid.png",plotDir.Data()));

TCanvas* cZNC_X_centroid = new TCanvas("cZNC_X_centroid","cZNC_X_centroid", 0,0,750,900);
hznc_Xcentroid->Draw();
cZNC_X_centroid->Print(Form("%s/cZNC_X_centroid.png",plotDir.Data()));

TCanvas* cZNC_Y_centroid = new TCanvas("cZNC_Y_centroid","cZNC_Y_centroid", 0,0,750,900);
hznc_Ycentroid->Draw();
cZNC_Y_centroid->Print(Form("%s/cZNC_Y_centroid.png",plotDir.Data()));

//---------------------------------------------------------------------------------
//timing
//---------------------------------------------------------------------------------
TCanvas* cTimingSum = new TCanvas("cTimingSum","cTimingSum",0,0,750,900);
hzn_TDC_Sum->Draw();
cTimingSum->Print(Form("%s/cTimingSum.png",plotDir.Data()));

TCanvas* cTimingDiff = new TCanvas("cTimingDiff","cTimingDiff",0,0,750,900);
hzn_TDC_Diff->Draw();
cTimingDiff->Print(Form("%s/cTimingDiff.png",plotDir.Data()));

//----------------------------------------------------------------------
//out
//----------------------------------------------------------------------
tree->Branch("run",&runNumber,"runNumber/I");
tree->Branch("ZNC_mean_value",&ZNC_mean,"ZNC_mean/D");
tree->Branch("ZNA_mean_value",&ZNA_mean,"ZNA_mean/D");
tree->Branch("ZPC_mean_value",&ZPC_mean,"ZPC_mean/D");
tree->Branch("ZPA_mean_value",&ZPA_mean,"ZPA_mean/D");
tree->Branch("ZNC_mean_uncalib",&ZNCuncalib_mean,"ZNCuncalib_mean/D");
tree->Branch("ZNA_mean_uncalib",&ZNAuncalib_mean,"ZNAuncalib_mean/D");
tree->Branch("ZPC_mean_uncalib",&ZPCuncalib_mean,"ZPCuncalib_mean/D");
tree->Branch("ZPA_mean_uncalib",&ZPAuncalib_mean,"ZPAuncalib_mean/D");
tree->Branch("ZEM1_mean_value",&ZEM1_mean,"ZEM1_mean/D");
tree->Branch("ZEM2_mean_value",&ZEM2_mean,"ZEM2_mean/D");
tree->Branch("ZNC_X_Centroid",&ZNC_XCent,"ZNC_XCent/D");
tree->Branch("ZNC_Y_Centroid",&ZNC_YCent,"ZNC_YCent/D");
tree->Branch("ZNA_X_Centroid",&ZNA_XCent,"ZNA_XCent/D");
tree->Branch("ZNA_Y_Centroid",&ZNA_YCent,"ZNA_YCent/D");
tree->Branch("ZNC_X_Centroid_Err",&ZNC_XCent_err,"ZNC_XCent_err/D");
tree->Branch("ZNC_Y_Centroid_Err",&ZNC_YCent_err,"ZNC_YCent_err/D");
tree->Branch("ZNA_X_Centroid_Err",&ZNA_XCent_err,"ZNA_XCent_err/D");
tree->Branch("ZNA_Y_Centroid_Err",&ZNA_YCent_err,"ZNA_YCent_err/D");
tree->Branch("ZN_TDC_Sum",&ZN_TDC_Sum,"ZN_TDC_Sum/D");
tree->Branch("ZN_TDC_Diff",&ZN_TDC_Diff,"ZN_TDC_Diff/D");
tree->Branch("ZN_TDC_Sum_Err",&ZN_TDC_Sum_err,"ZN_TDC_Sum_err/D");
tree->Branch("ZN_TDC_Diff_Err",&ZN_TDC_Diff_err,"ZN_TDC_Diff_err/D");
tree->Fill();

list.Add(cZNC_Spectra_Uncal);
list.Add(cZNA_Spectra_Uncal);
list.Add(cZPC_Spectra_Uncal);
list.Add(cZPA_Spectra_Uncal);
list.Add(cZEM1_Spectra);
list.Add(cZEM2_Spectra);
list.Add(cTimingSum);
list.Add(cTimingDiff);
list.Add(cZNA_X_centroid);
list.Add(cZNA_Y_centroid);
list.Add(cZNC_X_centroid);
list.Add(cZNC_Y_centroid);

TFile *fout;
fout = TFile::Open("prodQAhistos.root", "update");
if(!fout) fout = new TFile("prodQAhistos.root");
fout->cd();
list.Write();
tree->Write();
fout->Close();

return 0;

}
