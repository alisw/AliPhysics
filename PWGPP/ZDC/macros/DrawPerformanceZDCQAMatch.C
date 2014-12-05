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
gSystem->Load("libANALYSISalice.so");
gSystem->Load("libANALYSIScalib.so");
gSystem->Load("libTENDER.so");
gSystem->Load("libPWGPP.so");

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
  Printf("Invalid trending tree.");
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

char  outfilename[200]= "trending.root";
TString plotDir(".");

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

TH1D *hZNCpmcSpectrum = new TH1D("hZNCpmcSpectrum","hZNCpmcSpectrum",200.,0.,2000.);
TH1D *hZNApmcSpectrum = new TH1D("hZNApmcSpectrum","hZNApmcSpectrum",200.,0.,2000.);
TH1D *hZPCpmcSpectrum = new TH1D("hZPCpmcSpectrum","hZPCpmcSpectrum",200.,0.,2000.);  
TH1D *hZPApmcSpectrum = new TH1D("hZPApmcSpectrum","hZPApmcSpectrum",200.,0.,2000.);
TH1D *hZEM1Spectrum = new TH1D("hZEM1Spectrum","hZEM1Spectrum",200.,0.,2000.);  
TH1D *hZEM2Spectrum = new TH1D("hZEM2Spectrum","hZEM2Spectrum",200.,0.,2000.);  

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

TH1D *hznc = new TH1D("hznc","ZNC TC mean signal",3,-1,1);
hznc->GetXaxis()->SetRangeUser(-1.,1.);
hznc->SetDrawOption("EP");
hznc->SetMarkerStyle(20);
hznc->SetMarkerColor(kRed);
hznc->SetLineColor(kRed);   

TH1D *hzna = new TH1D("hzna","ZNA TC mean signal",3,-1,1);
hzna->GetXaxis()->SetRangeUser(-1.,1.);  
hzna->SetDrawOption("EP");
hzna->SetMarkerStyle(20);
hzna->SetMarkerColor(kRed);
hzna->SetLineColor(kRed);   

TH1D *hzpc = new TH1D("hzpc","ZPC TC mean signal",3,-1,1);
hzpc->GetXaxis()->SetRangeUser(-1.,1.); 
hzpc->SetDrawOption("EP");
hzpc->SetMarkerStyle(20);
hzpc->SetMarkerColor(kRed);
hzpc->SetLineColor(kRed);   

TH1D *hzpa = new TH1D("hzpa","ZPA TC mean signal",3,-1,1);
hzpa->GetXaxis()->SetRangeUser(-1.,1.);   
hzpa->SetDrawOption("EP");
hzpa->SetMarkerStyle(20);
hzpa->SetMarkerColor(kRed);
hzpa->SetLineColor(kRed);   

TH1D *hzem1 = new TH1D("hzem1","ZEM1 mean signal",3,-1,1);
hzem1->GetXaxis()->SetRangeUser(-1.,1.);   
hzem1->SetDrawOption("EP");
hzem1->SetMarkerStyle(20);
hzem1->SetMarkerColor(kRed);  
hzem1->SetLineColor(kRed); 

TH1D *hzem2 = new TH1D("hzem2","ZEM2 mean signal",3,-1,1);
hzem2->GetXaxis()->SetRangeUser(-1.,1.);     
hzem2->SetDrawOption("EP");
hzem2->SetMarkerStyle(20);
hzem2->SetMarkerColor(kRed);
hzem2->SetLineColor(kRed);   

TH1D *hzna_Xcentroid = new TH1D("hzna_Xcentroid","ZNA X centroid",3,-1,1);
hzna_Xcentroid->GetXaxis()->SetRangeUser(-1.,1.);       
hzna_Xcentroid->SetDrawOption("EP");
hzna_Xcentroid->SetMarkerStyle(20);
hzna_Xcentroid->SetMarkerColor(kRed);
hzna_Xcentroid->SetLineColor(kRed);     

TH1D *hzna_Ycentroid = new TH1D("hzna_Ycentroid","ZNA Y centroid",3,-1,1);
hzna_Ycentroid->GetXaxis()->SetRangeUser(-1.,1.);    
hzna_Ycentroid->SetDrawOption("EP");
hzna_Ycentroid->SetMarkerStyle(20);
hzna_Ycentroid->SetMarkerColor(kRed);
hzna_Ycentroid->SetLineColor(kRed);      

TH1D *hznc_Xcentroid = new TH1D("hznc_Xcentroid","ZNC X centroid",3,-1,1);
hznc_Xcentroid->GetXaxis()->SetRangeUser(-1.,1.);     
hznc_Xcentroid->SetDrawOption("EP");
hznc_Xcentroid->SetMarkerStyle(20);
hznc_Xcentroid->SetMarkerColor(kRed);
hznc_Xcentroid->SetLineColor(kRed); 

TH1D *hznc_Ycentroid = new TH1D("hznc_Ycentroid","ZNC Y centroid",3,-1,1);
hznc_Ycentroid->GetXaxis()->SetRangeUser(-1.,1.);     
hznc_Ycentroid->SetDrawOption("EP");
hznc_Ycentroid->SetMarkerStyle(20);
hznc_Ycentroid->SetMarkerColor(kRed);
hznc_Ycentroid->SetLineColor(kRed); 

TH1D *hzn_TDC_Sum = new TH1D("hzn_TDC_Sum","ZNC TDC + ZNA TDC",3,-1,1);
hzn_TDC_Sum->GetXaxis()->SetRangeUser(-1.,1.);   
hzn_TDC_Sum->SetDrawOption("EP");
hzn_TDC_Sum->SetMarkerStyle(20);
hzn_TDC_Sum->SetMarkerColor(kRed);
hzn_TDC_Sum->SetLineColor(kRed);   

TH1D *hzn_TDC_Diff = new TH1D("hzn_TDC_Diff","ZNC TDC - ZNA TDC",3,-1,1);
hzn_TDC_Diff->GetXaxis()->SetRangeUser(-1.,1.);    
hzn_TDC_Diff->SetDrawOption("EP");
hzn_TDC_Diff->SetMarkerStyle(20);
hzn_TDC_Diff->SetMarkerColor(kRed);
hzn_TDC_Diff->SetLineColor(kRed);   
  
char runlabel[10];

for (Int_t irun=0;irun<nRuns;irun++){  
  ttree->GetEntry(irun);
 }
    
sprintf(runlabel,"%i",runNumber);
    
//----------------------------------------------------------------------
//spectrum vs run
//----------------------------------------------------------------------

hZNCpmcSpectrum = dynamic_cast<TH1D*> (fin->Get("fhZNCpmc"));
hZNCpmcSpectrum->GetXaxis()->SetRangeUser(0.,2000);   
hZNCpmcSpectrum->Scale(1./hZNCpmcSpectrum->GetEntries());
hZNCpmcSpectrum->SetLineColor(kRed);
hZNCpmcSpectrum->SetLineWidth(2);    
hZNCpmcSpectrum->SetTitle("ZNC spectrum");
hZNCpmcSpectrum->SetXTitle("ZNC signal (ADC ch.)");

hZNApmcSpectrum = dynamic_cast<TH1D*> (fin->Get("fhZNApmc"));
hZNApmcSpectrum->GetXaxis()->SetRangeUser(0.,2000);     
hZNApmcSpectrum->Scale(1./hZNApmcSpectrum->GetEntries());
hZNApmcSpectrum->SetLineColor(kRed);
hZNApmcSpectrum->SetLineWidth(2);        
hZNApmcSpectrum->SetTitle("ZNA spectrum");
hZNApmcSpectrum->SetXTitle("ZNA signal (ADC ch.)");

hZPCpmcSpectrum = dynamic_cast<TH1D*> (fin->Get("fhZPCpmc"));
hZPCpmcSpectrum->GetXaxis()->SetRangeUser(0.,2000);       
hZPCpmcSpectrum->Scale(1./hZPCpmcSpectrum->GetEntries());
hZPCpmcSpectrum->SetLineColor(kRed);
hZPCpmcSpectrum->SetLineWidth(2); 
hZPCpmcSpectrum->SetTitle("ZPC spectrum");
hZPCpmcSpectrum->SetXTitle("ZPC signal (ADC ch.)");

hZPApmcSpectrum = dynamic_cast<TH1D*> (fin->Get("fhZPApmc"));
hZPApmcSpectrum->GetXaxis()->SetRangeUser(0.,2000);        
hZPApmcSpectrum->Scale(1./hZPApmcSpectrum->GetEntries());
hZPApmcSpectrum->SetLineColor(kRed);
hZPApmcSpectrum->SetLineWidth(2); 
hZPApmcSpectrum->SetTitle("ZPA spectrum");
hZPApmcSpectrum->SetXTitle("ZPA signal (ADC ch.)");    

hZEM1Spectrum = dynamic_cast<TH1D*> (fin->Get("fhZEM1Spectrum"));
hZEM1Spectrum->GetXaxis()->SetRangeUser(0.,2000);           
hZEM1Spectrum->Scale(1./hZEM1Spectrum->GetEntries());
hZEM1Spectrum->SetLineColor(kRed);
hZEM1Spectrum->SetLineWidth(2);     
hZEM1Spectrum->SetTitle("ZEM1 spectrum");
hZEM1Spectrum->SetXTitle("ZEM1 signal (ADC ch.)");    

hZEM2Spectrum = dynamic_cast<TH1D*> (fin->Get("fhZEM2Spectrum"));
hZEM2Spectrum->GetXaxis()->SetRangeUser(0.,2000);      
hZEM2Spectrum->Scale(1./hZEM2Spectrum->GetEntries());
hZEM2Spectrum->SetLineColor(kRed);
hZEM2Spectrum->SetLineWidth(2);      
hZEM2Spectrum->SetTitle("ZEM2 spectrum");
hZEM2Spectrum->SetXTitle("ZEM2 signal (ADC ch.)");    
          
//----------------------------------------------------------------------    
//variables vs run
//----------------------------------------------------------------------
hzna->SetBinContent(2,ZNA_mean);
hzna->GetXaxis()->SetBinLabel(2,runlabel);
hzna->GetXaxis()->SetLabelSize(0.05);

hzpc->SetBinContent(2,ZPC_mean);
hzpc->GetXaxis()->SetBinLabel(2,runlabel);
hzpc->GetXaxis()->SetLabelSize(0.05);    

hzpa->SetBinContent(2,ZPA_mean);
hzpa->GetXaxis()->SetBinLabel(2,runlabel);
hzpa->GetXaxis()->SetLabelSize(0.05);      

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
    
//----------------------------------------------------------------------
//spectra
//----------------------------------------------------------------------  
  
TCanvas* cZNC_Spectra = new TCanvas("cZNC_Spectra","cZNC_Spectra",50,50,1200,900); gPad->SetLogy();
hZNCpmcSpectrum->Draw();
cZNC_Spectra->Print(Form("%s/cZNC_Spectra.png",plotDir.Data()));

TCanvas* cZNA_Spectra = new TCanvas("cZNA_Spectra","cZNA_Spectra",50,50,1200,900); gPad->SetLogy();
hZNApmcSpectrum->Draw();
cZNA_Spectra->Print(Form("%s/cZNA_Spectra.png",plotDir.Data()));

TCanvas* cZPC_Spectra = new TCanvas("cZPC_Spectra","cZPC_Spectra",50,50,1200,900); gPad->SetLogy();
hZPCpmcSpectrum->Draw();
cZPC_Spectra->Print(Form("%s/cZPC_Spectra.png",plotDir.Data()));

TCanvas* cZPA_Spectra = new TCanvas("cZPA_Spectra","cZPA_Spectra",50,50,1200,900); gPad->SetLogy();
hZPApmcSpectrum->Draw();
cZPA_Spectra->Print(Form("%s/cZPA_Spectra.png",plotDir.Data()));

TCanvas* cZEM1_Spectra = new TCanvas("cZEM1_Spectra","cZEM1_Spectra",50,50,1200,900); gPad->SetLogy();
hZEM1Spectrum->Draw();
cZEM1_Spectra->Print(Form("%s/cZEM1_Spectra.png",plotDir.Data()));

TCanvas* cZEM2_Spectra = new TCanvas("cZEM2_Spectra","cZEM2_Spectra",50,50,1200,900); gPad->SetLogy();
hZEM2Spectrum->Draw();
cZEM2_Spectra->Print(Form("%s/cZEM2_Spectra.png",plotDir.Data()));


//---------------------------------------------------------------------------------------------------
//variables
//---------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------
//means
//---------------------------------------------------------------------------------------------------
TCanvas* cZNC_Mean_Values = new TCanvas("cZNC_Mean_Values","cZNC_Mean_Values", 50,50,750,550);
hznc->Draw("ep");
hznc->SetBinContent(2,ZNC_mean);
hznc->GetXaxis()->SetBinLabel(2,runlabel);  
hznc->GetXaxis()->SetLabelSize(0.05);
cZNC_Mean_Values->Print(Form("%s/cZNC_Mean_Values.png",plotDir.Data()));
  
TCanvas* cZNA_Mean_Values = new TCanvas("cZNA_Mean_Values","cZNA_Mean_Values", 50,50,750,550);
hzna->Draw("ep");
cZNA_Mean_Values->Print(Form("%s/cZNA_Mean_Values.png",plotDir.Data()));

TCanvas* cZPC_Mean_Values = new TCanvas("cZPC_Mean_Values","cZPC_Mean_Values", 50,50,750,550);
hzpc->Draw("ep");
cZPC_Mean_Values->Print(Form("%s/cZPC_Mean_Values.png",plotDir.Data()));

TCanvas* cZPA_Mean_Values = new TCanvas("cZPA_Mean_Values","cZPA_Mean_Values", 50,50,750,550);
hzpa->Draw("ep");
cZPA_Mean_Values->Print(Form("%s/cZPA_Mean_Values.png",plotDir.Data()));

TCanvas* cZEM1_Mean_Values = new TCanvas("cZEM1_Mean_Values","cZEM1_Mean_Values", 50,50,750,550);
hzem1->Draw("ep");
cZEM1_Mean_Values->Print(Form("%s/cZEM1_Mean_Values.png",plotDir.Data()));

TCanvas* cZEM2_Mean_Values = new TCanvas("cZEM2_Mean_Values","cZEM2_Mean_Values", 50,50,750,550);
hzem2->Draw("ep");
cZEM2_Mean_Values->Print(Form("%s/cZEM2_Mean_Values.png",plotDir.Data()));

//---------------------------------------------------------------------------------------------------
//centroids
//---------------------------------------------------------------------------------------------------
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

//---------------------------------------------------------------------------------
//timing
//---------------------------------------------------------------------------------
TCanvas* cTimingSum = new TCanvas("cTimingSum","cTimingSum",50,50,750,550);
hzn_TDC_Sum->Draw();
cTimingSum->Print(Form("%s/cTimingSum.png",plotDir.Data()));

TCanvas* cTimingDiff = new TCanvas("cTimingDiff","cTimingDiff",50,50,750,550);
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
tree.Write();
fout->Close();

return 0;  
    
}
