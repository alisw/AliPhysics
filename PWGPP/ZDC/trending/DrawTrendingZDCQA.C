/******************************************************************************************************************************************
Contact person: Marco Leoncino (leoncino@to.infn.it)
Macro to draw the ZDC QA trending plots by accessing the std tree. To be mainly used with the automatic scripts to fill the QA repository.
Launch with aliroot -x -b -q "DrawTrendingZDCQA.C" The macro produces one png file for trending variables and a file with the histograms.
******************************************************************************************************************************************/
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

Int_t DrawTrendingZDCQA(TString mergedTrendFile = "trending.root", // trending tree file 
                        Bool_t displayAll = kFALSE)                // set to kTRUE to display trending for expert plots
{
  //
  //reads merged trending.root file and draws trending plots from tree
  //
  if (!mergedTrendFile) {
    Printf("Cannot open merged trend file with ZDC QA");
    return 1;
  }
  
  /*set graphic style*/
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTextFont(42);
  gStyle->SetStatColor(kWhite); 
  gStyle->SetStatBorderSize(1);
  gStyle->SetOptStat(0);  
  TGaxis::SetMaxDigits(3);
 
  char  outfilename[200]= "ProductionQA.hist.root";
  TString plotDir(".");
  //TString plotDir(Form("PlotsTrending"));
  //gSystem->Exec(Form("mkdir %s",plotDir.Data()));
  
  legend = new TLegend(0.9,0.1,1.0,0.9);
  legend->SetFillColor(kWhite);

  Int_t runNumber=0;
  Double_t ZNC_mean=0;
  Double_t ZNA_mean=0;
  Double_t ZPA_mean=0;
  Double_t ZPC_mean=0;
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

  TFile * fin = TFile::Open(mergedTrendFile.Data());
  TTree * ttree = (TTree*) fin->Get("trending");
  if (!ttree){
    Printf("Invalid trending tree.");
    return 2;
  }
 
  Int_t nRuns=ttree->GetEntries();
  printf("Number of runs to be analyzed:%d\n",nRuns);
  TList list;
 
  TH1D *hZNCpmcSpectrum = new TH1D("hZNCpmcSpectrum","ZNC common tower mean ADC channel vs Run",200.,0.,2000.);
  TH1D *hZNApmcSpectrum = new TH1D("hZNApmcSpectrum","ZNA common tower mean ADC channel vs Run",200.,0.,2000.);
  TH1D *hZPCpmcSpectrum = new TH1D("hZPCpmcSpectrum","ZPC common tower mean ADC channel vs Run",200.,0.,2000.);  
  TH1D *hZPApmcSpectrum = new TH1D("hZPApmcSpectrum","ZPA common tower mean ADC channel vs Run",200.,0.,2000.);
  TH1D *hZEM1Spectrum = new TH1D("hZEM1Spectrum","ZEM1 mean ADC channel vs Run",200.,0.,2000.);  
  TH1D *hZEM2Spectrum = new TH1D("hZEM2Spectrum","ZEM2 mean ADC channel vs Run",200.,0.,2000.);  
  
  ttree->SetBranchAddress("run",&runNumber);
  ttree->SetBranchAddress("ZNC_mean_value",&ZNC_mean);
  ttree->SetBranchAddress("ZNA_mean_value",&ZNA_mean);
  ttree->SetBranchAddress("ZPC_mean_value",&ZPC_mean);
  ttree->SetBranchAddress("ZPA_mean_value",&ZPA_mean); 
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
        
  TH1D *hznc = new TH1D("hznc","ZNC common tower mean ADC channel vs Run",nRuns,0.,nRuns);
  hznc->SetDrawOption("E");
  hznc->SetDrawOption("E");
  
  TH1D *hzna = new TH1D("hzna","ZNA common tower mean ADC channel vs Run",nRuns,0.,nRuns);
  hzna->SetDrawOption("E");
  hzna->SetMarkerStyle(20);  
  
  TH1D *hzpc = new TH1D("hzpc","ZPC common tower mean ADC channel vs Run",nRuns,0.,nRuns);
  hzpc->SetDrawOption("E");
  hzpc->SetMarkerStyle(20);
  
  TH1D *hzpa = new TH1D("hzpa","ZPA common tower mean ADC channel vs Run",nRuns,0.,nRuns);
  hzpa->SetDrawOption("E");
  hzpa->SetMarkerStyle(20);  
  
  TH1D *hzem1 = new TH1D("hzem1","ZEM1 mean ADC channel vs Run",nRuns,0.,nRuns);
  hzem1->SetDrawOption("E");
  hzem1->SetMarkerStyle(20);  
  
  TH1D *hzem2 = new TH1D("hzem2","ZEM2 mean ADC channel vs Run",nRuns,0.,nRuns);
  hzem2->SetDrawOption("E");
  hzem2->SetMarkerStyle(20);
  
  TH1D *hzna_Xcentroid = new TH1D("hzna_Xcentroid","ZNA X centroid vs Run",nRuns,0.,nRuns);
  hzna_Xcentroid->SetDrawOption("E");
  hzna_Xcentroid->SetMarkerStyle(20);

  TH1D *hzna_Ycentroid = new TH1D("hzna_Ycentroid","ZNA Y centroid vs Run",nRuns,0.,nRuns);
  hzna_Ycentroid->SetDrawOption("E");
  hzna_Ycentroid->SetMarkerStyle(20);
  
  TH1D *hznc_Xcentroid = new TH1D("hznc_Xcentroid","ZNC X centroid vs Run",nRuns,0.,nRuns);
  hznc_Xcentroid->SetDrawOption("E");
  hznc_Xcentroid->SetMarkerStyle(20);

  TH1D *hznc_Ycentroid = new TH1D("hznc_Ycentroid","ZNC Y centroid vs Run",nRuns,0.,nRuns);
  hznc_Ycentroid->SetDrawOption("E");
  hznc_Ycentroid->SetMarkerStyle(20);
  
  TH1D *hzn_TDC_Sum = new TH1D("hzn_TDC_Sum","ZNC TDC + ZNA TDC",nRuns,0.,nRuns);
  hzn_TDC_Sum->SetDrawOption("E");
  hzn_TDC_Sum->SetMarkerStyle(20);
  
  TH1D *hzn_TDC_Diff = new TH1D("hzn_TDC_Diff","ZNC TDC - ZNA TDC",nRuns,0.,nRuns);
  hzn_TDC_Diff->SetDrawOption("E");
  hzn_TDC_Diff->SetMarkerStyle(20);   
  
  char runlabel[10];
  int  icol = 0;

  for (Int_t irun=0;irun<nRuns;irun++){
    
    ttree->GetEntry(irun);
    sprintf(runlabel,"%i",runNumber);
    icol = 99*gRandom->Rndm();
    
    //----------------------------------------------------------------------
    //spectra vs run
    //----------------------------------------------------------------------
    
    hZNCpmcSpectrum = dynamic_cast<TH1D*> (fin->Get("fhZNCpmc"));
    hZNCpmcSpectrum->Scale(1./hZNCpmcSpectrum->GetEntries());
    hZNCpmcSpectrum->SetLineColor(icol);
    hZNCpmcSpectrum->SetTitle("ZNC spectra as a function of Run Number");
    hZNCpmcSpectrum->SetXTitle("ZNC signal (ADC ch.)");
    
    hZNApmcSpectrum = dynamic_cast<TH1D*> (fin->Get("fhZNApmc"));
    hZNApmcSpectrum->Scale(1./hZNApmcSpectrum->GetEntries());
    hZNApmcSpectrum->SetLineColor(icol);
    hZNApmcSpectrum->SetTitle("ZNA spectra as a function of Run Number");
    hZNApmcSpectrum->SetXTitle("ZNA signal (ADC ch.)");
    
    hZPCpmcSpectrum = dynamic_cast<TH1D*> (fin->Get("fhZPCpmc"));
    hZPCpmcSpectrum->Scale(1./hZPCpmcSpectrum->GetEntries());
    hZPCpmcSpectrum->SetLineColor(icol);
    hZPCpmcSpectrum->SetTitle("ZPC spectra as a function of Run Number");
    hZPCpmcSpectrum->SetXTitle("ZPC signal (ADC ch.)");
    
    hZPApmcSpectrum = dynamic_cast<TH1D*> (fin->Get("fhZPApmc"));
    hZPApmcSpectrum->Scale(1./hZPApmcSpectrum->GetEntries());
    hZPApmcSpectrum->SetLineColor(icol);
    hZPApmcSpectrum->SetTitle("ZPA spectra as a function of Run Number");
    hZPApmcSpectrum->SetXTitle("ZPA signal (ADC ch.)");    
    
    hZEM1Spectrum = dynamic_cast<TH1D*> (fin->Get("fhZEM1Spectrum"));
    hZEM1Spectrum->Scale(1./hZEM1Spectrum->GetEntries());
    hZEM1Spectrum->SetLineColor(icol);
    hZEM1Spectrum->SetTitle("ZEM1 spectra as a function of Run Number");
    hZEM1Spectrum->SetXTitle("ZEM1 signal (ADC ch.)");    
    
    hZEM2Spectrum = dynamic_cast<TH1D*> (fin->Get("fhZEM2Spectrum"));
    hZEM2Spectrum->Scale(1./hZEM2Spectrum->GetEntries());
    hZEM2Spectrum->SetLineColor(icol);
    hZEM2Spectrum->SetTitle("ZEM2 spectra as a function of Run Number");
    hZEM2Spectrum->SetXTitle("ZEM2 signal (ADC ch.)");    
     
    legend->AddEntry(hZPApmcSpectrum, runlabel, "l");    
      
    //----------------------------------------------------------------------    
    //variables vs run
    //----------------------------------------------------------------------
    
    hznc->SetBinContent(irun+1, ZNC_mean);
    hznc->GetXaxis()->SetBinLabel(irun+1,runlabel);  

    hzna->SetBinContent(irun+1, ZNA_mean);
    hzna->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
    hzpc->SetBinContent(irun+1, ZPC_mean);
    hzpc->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
    hzpa->SetBinContent(irun+1, ZPA_mean);
    hzpa->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
    hzem1->SetBinContent(irun+1, ZEM1_mean);
    hzem1->GetXaxis()->SetBinLabel(irun+1,runlabel);
    
    hzem2->SetBinContent(irun+1, ZEM2_mean);
    hzem2->GetXaxis()->SetBinLabel(irun+1,runlabel); 
    
    hzna_Xcentroid->SetBinContent(irun+1, ZNA_XCent);
    hzna_Xcentroid->SetBinError(irun+1, ZNA_XCent_err);
    hzna_Xcentroid->GetXaxis()->SetBinLabel(irun+1,runlabel);    
    
    hzna_Ycentroid->SetBinContent(irun+1, ZNA_YCent);
    hzna_Ycentroid->SetBinError(irun+1, ZNA_YCent_err); 
    hzna_Ycentroid->GetXaxis()->SetBinLabel(irun+1,runlabel); 
    
    hznc_Xcentroid->SetBinContent(irun+1, ZNC_XCent);
    hznc_Xcentroid->SetBinError(irun+1, ZNC_XCent_err);
    hznc_Xcentroid->GetXaxis()->SetBinLabel(irun+1,runlabel); 
    
    hznc_Ycentroid->SetBinContent(irun+1, ZNC_YCent);
    hznc_Ycentroid->SetBinError(irun+1, ZNC_YCent_err);
    hznc_Ycentroid->GetXaxis()->SetBinLabel(irun+1,runlabel); 
    
    hzn_TDC_Sum->SetBinContent(irun+1, ZN_TDC_Sum);
    hzn_TDC_Sum->SetBinError(irun+1, ZN_TDC_Sum_err);
    hzn_TDC_Sum->GetXaxis()->SetBinLabel(irun+1,runlabel); 
    
    hzn_TDC_Diff->SetBinContent(irun+1, ZN_TDC_Diff);
    hzn_TDC_Diff->SetBinError(irun+1, ZN_TDC_Diff_err); 
    hzn_TDC_Diff->GetXaxis()->SetBinLabel(irun+1,runlabel);     
  }

//----------------------------------------------------------------------
//spectra vs run
//----------------------------------------------------------------------  
  
TCanvas* cZNC_Spectra = new TCanvas("cZNC_Spectra","cZNC_Spectra",50,50,1200,900);
hZNCpmcSpectrum->Draw();
cZNC_Spectra->Print(Form("%s/cZNC_Spectra.png",plotDir.Data()));
gPad->SetLogy();
legend->Draw("same");

TCanvas* cZNA_Spectra = new TCanvas("cZNA_Spectra","cZNA_Spectra",50,50,1200,900);
hZNApmcSpectrum->Draw();
cZNA_Spectra->Print(Form("%s/cZNA_Spectra.png",plotDir.Data()));
gPad->SetLogy();
legend->Draw("same");

TCanvas* cZPC_Spectra = new TCanvas("cZPC_Spectra","cZPC_Spectra",50,50,1200,900);
hZPCpmcSpectrum->Draw();
cZPC_Spectra->Print(Form("%s/cZPC_Spectra.png",plotDir.Data()));
gPad->SetLogy();
legend->Draw("same");

TCanvas* cZPA_Spectra = new TCanvas("cZPA_Spectra","cZPA_Spectra",50,50,1200,900);
hZPApmcSpectrum->Draw();
cZPA_Spectra->Print(Form("%s/cZPA_Spectra.png",plotDir.Data()));
gPad->SetLogy();
legend->Draw("same");

TCanvas* cZEM1_Spectra = new TCanvas("cZEM1_Spectra","cZEM1_Spectra",50,50,1200,900);
hZEM1Spectrum->Draw();
cZEM1_Spectra->Print(Form("%s/cZEM1_Spectra.png",plotDir.Data()));
gPad->SetLogy();
legend->Draw("same");

TCanvas* cZEM2_Spectra = new TCanvas("cZEM2_Spectra","cZEM2_Spectra",50,50,1200,900);
hZEM2Spectrum->Draw();
cZEM2_Spectra->Print(Form("%s/cZEM2_Spectra.png",plotDir.Data()));
gPad->SetLogy();
legend->Draw("same");

//----------------------------------------------------------------------
//variables vs run
//---------------------------------------------------------------------- 
TCanvas* cZNC_mean_values = new TCanvas("cZNC_mean_values","cZNC_mean_values", 50,50,750,550);
hznc->Draw();
cZNC_mean_values->Print(Form("%s/cZNC_mean_values.png",plotDir.Data()));
  
TCanvas* cZNA_mean_values = new TCanvas("cZNA_mean_values","cZNA_mean_values", 50,50,750,550);
hzna->Draw();
cZNA_mean_values->Print(Form("%s/cZNA_mean_values.png",plotDir.Data()));

TCanvas* cZPC_mean_values = new TCanvas("cZPC_mean_values","cZPC_mean_values", 50,50,750,550);
hzpc->Draw();
cZPC_mean_values->Print(Form("%s/cZNP_mean_values.png",plotDir.Data()));

TCanvas* cZPA_mean_values = new TCanvas("cZPA_mean_values","cZPA_mean_values", 50,50,750,550);
hzpa->Draw();
cZPA_mean_values->Print(Form("%s/cZPA_mean_values.png",plotDir.Data()));

TCanvas* cZNA_X_centroid = new TCanvas("cZNA_X_centroid","cZNA_X_centroid", 50,50,750,550);
hzna_Xcentroid->Draw();
cZNA_X_centroid->Print(Form("%s/cZNA_X_centroid.png",plotDir.Data()));

TCanvas* cZNA_Y_centroid = new TCanvas("cZNA_Y_centroid","cZNA_Y_centroid", 50,50,750,550);
hzna_Ycentroid->Draw();
cZNA_Y_centroid->Print(Form("%s/cZNA_Y_centroid.png",plotDir.Data()));

TCanvas* cZNC_X_centroid = new TCanvas("cZNC_X_centroid","cZNC_X_centroid", 50,50,750,550);
hznc_Xcentroid->Draw();
cZNC_X_centroid->Print(Form("%s/cZNC_X_centroid.png",plotDir.Data()));

TCanvas* cZNC_Y_centroid = new TCanvas("cZNC_Y_centroid","cZNC_Y_centroid", 50,50,750,550);
hznc_Ycentroid->Draw();
cZNC_Y_centroid->Print(Form("%s/cZNC_Y_centroid.png",plotDir.Data()));

list.Add(cZNC_Spectra);
list.Add(cZNA_Spectra);
list.Add(cZPC_Spectra);
list.Add(cZPA_Spectra);
list.Add(cZEM1_Spectra);
list.Add(cZEM2_Spectra);
list.Add(hznc);
list.Add(hzna);
list.Add(hzpc);  
list.Add(hzpa);
list.Add(hzem1);  
list.Add(hzem2);
list.Add(hzna_Xcentroid);   
list.Add(hzna_Ycentroid);
list.Add(hznc_Xcentroid);   
list.Add(hznc_Ycentroid);
list.Add(hzn_TDC_Sum);   
list.Add(hzn_TDC_Diff);

TFile * fout=new TFile(outfilename,"recreate");
fout->cd();
list.Write();
fout->Close();

return 0;  
    
}
