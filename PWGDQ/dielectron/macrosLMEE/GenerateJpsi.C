#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TLine.h"
#include "TError.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TMinuit.h"
#include "TTree.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "AliAnalysisTaskJpsi.h"

void GenerateJpsi(TString resolution_file = "resolution_133_LHC18g3_Fast_cent_woSDD_deltaXvsP_ROOT5.root", TString param_file ="pp502TeV_22_01_2018.root", TString param_std ="5TeV_default", TString param_low ="5TeV_default", TString param_high ="5TeV_default", Double_t ptmin=0.2, Double_t ptmax=10.)
{

  //gROOT->ProcessLine(".L AliAnalysisTaskJPsi.cxx+g");
  AliAnalysisTaskJpsi *analysis = new AliAnalysisTaskJpsi();  
  
  // 
  // Resolution
  // 
  TFile fRes(resolution_file.Data());
  if (fRes.IsOpen() && ((TObjArray*)fRes.Get("ptSlices"))!=0x0) { // Smearing Run 1: momentum and opening angle
    analysis->SetResolutionP((TObjArray*)fRes.Get("ptSlices"), kFALSE);
  }
  else { // New Smearing Run 2
    analysis->ReadResoFile(&fRes);
  }

  //
  // Parametrizations
  //
  TFile fParam(param_file.Data());
  TDirectoryFile *d = (TDirectoryFile *) fParam.Get(param_std.Data());
  if(d) {
    TF1 *fjpsi = (TF1 *) d->Get("443_pt");
    analysis->SetScaleFunction(fjpsi);
  } else printf("No central param\n");
  TDirectoryFile *d_low = (TDirectoryFile *) fParam.Get(param_low.Data());
  if(d_low) {
    TF1 *fjpsi_low = (TF1 *) d_low->Get("443_pt");
    analysis->SetScaleFunctionLow(fjpsi_low);
  } else printf("No clow param\n");
  TDirectoryFile *d_high = (TDirectoryFile *) fParam.Get(param_high.Data());
  if(d_high) {
    TF1 *fjpsi_high = (TF1 *) d_high->Get("443_pt");
    analysis->SetScaleFunctionHigh(fjpsi_high);
  } else printf("No high param\n");

  //
  // Settings
  //
  analysis->SetPtRange(ptmin,ptmax);
  analysis->SetEtaRange(-0.8,0.8);

  
  analysis->Init();
  

  //
  // Input tree
  //
  TFile *fTree = new TFile("./input/rawTree_jpsi.root","READ");
  TTree *tree = static_cast<TTree*> (fTree->Get("JPsi_ee"));
  analysis->ConnectTree(tree);

  //
  // Analysis
  //
  analysis->DetermineWeights();
  analysis->Fill();
  
  //
  // Write in File
  //
  TFile *fOut = new TFile("./output/Jpsi_cocktail.root","RECREATE");
  fOut->cd();
  TList *l = analysis->GetOutputList();
  l->Write(Form("%s",param_std.Data()),TObject::kSingleKey);
  //l = analysis->GetOutputListLow();
  //l->Write(Form("%s",param_low.Data()),TObject::kSingleKey);
  //l = analysis->GetOutputListHigh();
  //l->Write(Form("%s",param_high.Data()),TObject::kSingleKey);
  fOut->Close();

}
