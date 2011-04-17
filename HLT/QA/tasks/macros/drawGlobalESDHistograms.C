// $Id$
/*
 * Plotting macro for comparing offline- and HLT- ESD trees from
 * HLT-OFFLINE-GLOBAL-comparison.root produced using $ALICE_ROOT/HLT/QA/tasks/AliAnalysisTaskHLT.*
 *
 * Usage: aliroot drawGlobalESDHistograms.C'("HLT-OFFLINE-GLOBAL-comparison.root")'
 *
 * or aliroot drawGlobalESDHistograms.C++ in compiled mode
 *
 * It saves the canvas with the output histograms in a png and a ROOT file.
 *
 * @ingroup alihlt_qa
 * @author Camilla.Stokkevag@student.uib.no, Kalliopi.Kanaki@ift.uib.no 
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TList.h"
#include "TCanvas.h"
#include "TText.h"
#include "TPaveStats.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include <iostream>
#include <cstdlib>
using std::endl;
#endif

// --------------------- forward declerations --------------//

void plot(TH1F *h1, TH1F *h2);
void printStats(TH1F *h1);
void printStats(TH1F *h1, TH1F *h2);

//==========================================================//

void drawGlobalESDHistograms(const char* filename="HLT-OFFLINE-GLOBAL-comparison.root"){
 
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat("emr");
 gStyle->SetTitleX(gStyle->GetPadLeftMargin());

 TFile *file = TFile::Open(filename); 
 if(!file || file->IsZombie()) {
    printf("file %s does not exist or there is an error opening it\n", filename);
    return;
 }

 TList *list = (TList*)file->Get("global_histograms");
 if(!list){
    printf("No list %s contained in your input file\n", list->GetName()); 
    return; 
 }
 
 TText *hText = (TText*)list->FindObject("text");
 if(!hText) printf("No hText\n");

 TString folder = "GlobalTask_";
 folder += hText->GetTitle();
 folder.ReplaceAll(" ",""); 
 folder.ReplaceAll(",","_");
 gSystem->Exec("mkdir "+folder); // create a folder whose name contains run number and date of run

 TCanvas *c1 = new TCanvas("c1","track properties HLT vs. OFF",1200,700);
 c1->Divide(4,2);

 TH1F *h1 = NULL;
 TH1F *h2 = NULL;


 h1 = (TH1F*)list->FindObject("fMomentum_hlt"); if(!h1) { printf("Empty histogram fMomentum_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fMomentum_off"); if(!h2) { printf("Empty histogram fMomentum_off\n"); return; }
 h1->SetXTitle("p_{t} (GeV/c)"); 

 c1->cd(1);
 plot(h1,h2);

 TLegend *leg1 = new TLegend(0.6,0.2,0.8,0.5);
 leg1->SetFillColor(10);
 leg1->SetLineColor(10);
 leg1->AddEntry(h1,"HLT", "l");
 leg1->AddEntry(h2,"OFF", "l");
 leg1->Draw("same");

//------------------------------------------------- 

 h1 = (TH1F*)list->FindObject("fNcluster_hlt"); if(!h1) { printf("Empty histogram fNcluster_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fNcluster_off"); if(!h2) { printf("Empty histogram fNcluster_off\n"); return; }
 h1->SetXTitle("TPC clusters per track");

 c1->cd(2);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)list->FindObject("fEta_hlt"); if(!h1) { printf("Empty histogram fEta_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fEta_off"); if(!h2) { printf("Empty histogram fEta_off\n"); return; }
 h1->SetXTitle("#eta");

 c1->cd(3);
 plot(h1,h2);

//-------------------------------------------------

 h1 = (TH1F*)list->FindObject("fPhi_hlt"); if(!h1) { printf("Empty histogram fPhi_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fPhi_off"); if(!h2) { printf("Empty histogram fPhi_off\n"); return; }
 h1->SetXTitle("#phi (deg)");

 c1->cd(4);
 plot(h1,h2);

//-------------------------------------------------

 h1 = (TH1F*)list->FindObject("fDCAr_hlt"); if(!h1) { printf("Empty histogram fDCAr_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fDCAr_off"); if(!h2) { printf("Empty histogram fDCAr_off\n"); return; }
 h1->SetXTitle("DCAr (cm)");
 
 c1->cd(5);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)list->FindObject("fDCAz_hlt"); if(!h1) { printf("Empty histogram fDCAz_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fDCAz_off"); if(!h2) { printf("Empty histogram fDCAz_off\n"); return; }
 h1->SetXTitle("DCAz (cm)");
 
 c1->cd(6);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)list->FindObject("fCharge_hlt"); if(!h1) { printf("Empty histogram fCharge_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fCharge_off"); if(!h2) { printf("Empty histogram fCharge_off\n"); return; }

 c1->cd(7);
 plot(h1,h2);

//------------------------------------------------- 
 
 h1 = (TH1F*)list->FindObject("fNITScluster_hlt"); if(!h1) { printf("Empty histogram fNITScluster_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fNITScluster_off"); if(!h2) { printf("Empty histogram fNITScluster_off\n"); return; }
 h1->SetXTitle("ITS clusters per track");

 c1->cd(8);
 plot(h1,h2);

//============= TRACK PROPERTIES WITH CUTS ===============//

 TCanvas *c4 = new TCanvas("c4","track properties HLT vs. OFF",1200,700);
 c4->Divide(4,2);

 h1 = (TH1F*)list->FindObject("fMomentum_hlt");    if(!h1) { printf("Empty histogram fMomentum_hlt\n");    return; }
 h2 = (TH1F*)list->FindObject("fMomentum_hltcut"); if(!h2) { printf("Empty histogram fMomentum_hltcut\n"); return; }
 h1->SetXTitle("p_{t} (GeV/c)"); 

 c4->cd(1);
 plot(h1,h2);

 TLegend *leg2 = new TLegend(0.6,0.2,0.8,0.5);
 leg2->SetFillColor(10);
 leg2->SetLineColor(10);
 leg2->AddEntry(h1,"HLT", "l");
 leg2->AddEntry(h2,"HLT with cuts", "l");
 leg2->Draw("same");

//------------------------------------------------- 
 
 h1 = (TH1F*)list->FindObject("fNcluster_hlt");    if(!h1) { printf("Empty histogram fNcluster_hlt\n");    return; }
 h2 = (TH1F*)list->FindObject("fNcluster_hltcut"); if(!h2) { printf("Empty histogram fNcluster_hltcut\n"); return; }
 h1->SetXTitle("TPC clusters per track");

 c4->cd(2);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)list->FindObject("fEta_hlt");    if(!h1) { printf("Empty histogram fEta_hlt\n");    return; }
 h2 = (TH1F*)list->FindObject("fEta_hltcut"); if(!h2) { printf("Empty histogram fEta_hltcut\n"); return; }
 h1->SetXTitle("#eta");

 c4->cd(3);
 plot(h1,h2);

//-------------------------------------------------

 h1 = (TH1F*)list->FindObject("fPhi_hlt");    if(!h1) { printf("Empty histogram fPhi_hlt\n");    return; }
 h2 = (TH1F*)list->FindObject("fPhi_hltcut"); if(!h2) { printf("Empty histogram fPhi_hltcut\n"); return; }
 h1->SetXTitle("#phi (deg)");

 c4->cd(4);
 plot(h1,h2);

//-------------------------------------------------

 h1 = (TH1F*)list->FindObject("fDCAr_hlt");    if(!h1) { printf("Empty histogram fDCAr_hlt\n");    return; }
 h2 = (TH1F*)list->FindObject("fDCAr_hltcut"); if(!h2) { printf("Empty histogram fDCAr_hltcut\n"); return; }
 h1->SetXTitle("DCAr (cm)");
 
 c4->cd(5);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)list->FindObject("fDCAz_hlt");    if(!h1) { printf("Empty histogram fDCAz_hlt\n");    return; }
 h2 = (TH1F*)list->FindObject("fDCAz_hltcut"); if(!h2) { printf("Empty histogram fDCAz_hltcut\n"); return; }
 h1->SetXTitle("DCAz (cm)");
 
 c4->cd(6);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)list->FindObject("fCharge_hlt");    if(!h1) { printf("Empty histogram fCharge_hlt\n");    return; }
 h2 = (TH1F*)list->FindObject("fCharge_hltcut"); if(!h2) { printf("Empty histogram fCharge_hltcut\n"); return; }
 h1->SetXTitle("polarity"); 

 c4->cd(7);
 plot(h1,h2);

//------------------------------------------------- 
 
 h1 = (TH1F*)list->FindObject("fNITScluster_hlt");    if(!h1) { printf("Empty histogram fNITScluster_hlt\n");    return; }
 h2 = (TH1F*)list->FindObject("fNITScluster_hltcut"); if(!h2) { printf("Empty histogram fNITScluster_hltcut\n"); return; }
 h1->SetXTitle("ITS clusters per track");

 c4->cd(8);
 plot(h1,h2);

//============= EVENT PROPERTIES ===============//

 TCanvas *c2 = new TCanvas("c2","vertex event properties",1200,700);
 c2->Divide(3,2);

 h1 = (TH1F*)list->FindObject("fXvertex_hlt"); if(!h1) { printf("Empty histogram fXvertex_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fXvertex_off"); if(!h2) { printf("Empty histogram fXvertex_off\n"); return; }
 h1->SetXTitle("x (cm)");

 c2->cd(1);
 plot(h1,h2);
 leg1->Draw("same");

//------------------------------------------------- 

 h1 = (TH1F*)list->FindObject("fYvertex_hlt"); if(!h1) { printf("Empty histogram fYvertex_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fYvertex_off"); if(!h2) { printf("Empty histogram fYvertex_off\n"); return; }
 h1->SetXTitle("y (cm)");
 
 c2->cd(2);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)list->FindObject("fZvertex_hlt"); if(!h1) { printf("Empty histogram fZvertex_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fZvertex_off"); if(!h2) { printf("Empty histogram fZvertex_off\n"); return; }
 h1->SetXTitle("z (cm)");

 c2->cd(3);
 plot(h1,h2);

//------------------------------------------------- 
 
 h1 = (TH1F*)list->FindObject("fSPDXvertex_hlt"); if(!h1) { printf("Empty histogram fSPDXvertex_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fSPDXvertex_off"); if(!h2) { printf("Empty histogram fSPDXvertex_off\n"); return; }
 h1->SetXTitle("x (cm)");

 c2->cd(4);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)list->FindObject("fSPDYvertex_hlt"); if(!h1) { printf("Empty histogram fSPDYvertex_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fSPDYvertex_off"); if(!h2) { printf("Empty histogram fSPDYvertex_off\n"); return; }
 h1->SetXTitle("y (cm)");
 
 c2->cd(5);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)list->FindObject("fSPDZvertex_hlt"); if(!h1) { printf("Empty histogram fSPDZvertex_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fSPDZvertex_off"); if(!h2) { printf("Empty histogram fSPDZvertex_off\n"); return; }
 h1->SetXTitle("z (cm)");

 c2->cd(6);
 plot(h1,h2);

//------------------------------------------------- 

 TCanvas *c3 = new TCanvas("c3","general event properties",1200,500);
 c3->Divide(3,1);

 h1 = (TH1F*)list->FindObject("fMult_hlt"); if(!h1) { printf("Empty histogram fMult_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fMult_off"); if(!h2) { printf("Empty histogram fMult_off\n"); return; }

 c3->cd(1);
 plot(h1,h2);
 leg1->Draw("same");

//------------------------------------------------- 

 h1 = (TH1F*)list->FindObject("fNcont_hlt"); if(!h1) { printf("Empty histogram fNcont_hlt\n"); return; }
 h2 = (TH1F*)list->FindObject("fNcont_off"); if(!h2) { printf("Empty histogram fNcont_off\n"); return; }

 c3->cd(2);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)list->FindObject("fV0cent"); if(!h1) { printf("Empty histogram fV0cent\n"); return; }
 c3->cd(3);
 h1->Draw();
 printStats(h1);

//------------------------------------------------- 

 c1->SaveAs(folder+"/track_properties.png");  
 c1->SaveAs(folder+"/track_properties.root");  
 c2->SaveAs(folder+"/vertex_event_properties.png");  
 c2->SaveAs(folder+"/vertex_event_properties.root");  
 c3->SaveAs(folder+"/general_event_properties.png");  
 c3->SaveAs(folder+"/general_event_properties.root");  
 c4->SaveAs(folder+"/HLT_track_properties_cuts.png");  
 c4->SaveAs(folder+"/HLT_track_properties_cuts.root");  
 return;	
}

void printStats(TH1F* h1, TH1F* h2){  
  
  gPad->Update();
  TPaveStats *st1 = (TPaveStats*)h1->FindObject("stats"); if(!st1) { printf("TPaveStats st1 is 0x0\n"); return; }	
  st1->SetLineColor(0);
 
  gPad->Update();
  TPaveStats *st2 = (TPaveStats*)h2->FindObject("stats"); if(!st2) { printf("TPaveStats st2 is 0x0\n"); return; }
  st2->SetY2NDC(st1->GetY1NDC()-0.05);
  st2->SetY1NDC(st2->GetY2NDC()-TMath::Abs(st1->GetY1NDC()-st1->GetY2NDC()));
  st2->SetLineColor(0);
  st2->SetTextColor(h2->GetLineColor());
  st2->Draw();  
  return;
}

void printStats(TH1F* h1){    
  gPad->Update();
  TPaveStats *st1 = (TPaveStats*)h1->FindObject("stats"); if(!st1) { printf("TPaveStats st1 is 0x0\n"); return; }	
  st1->SetLineColor(0);
  return;
}

void plot(TH1F *h1, TH1F *h2){ 
  //Y axis
  if(h1->GetMaximum() > h2->GetMaximum()) h2->SetMaximum(1.1*h1->GetMaximum());
  else h1->SetMaximum(1.1*h2->GetMaximum());
  
  h1->SetMinimum(0);
  h2->SetMinimum(0);
  h2->SetLineColor(2);
 
  // X axis  
  double xmin, xmax;  
  if(h1->GetBinLowEdge(1) > h2->GetBinLowEdge(1)) xmin = h1->GetBinLowEdge(1);
  else xmin = h2->GetBinLowEdge(1);
  if(h1->GetBinLowEdge(h1->GetNbinsX()+1) > h2->GetBinLowEdge(h1->GetNbinsX()+1)) xmax = h1->GetBinLowEdge(h1->GetNbinsX()+1);
  else xmax = h2->GetBinLowEdge(h2->GetNbinsX()+1);
  
  h2->SetAxisRange(xmin, xmax, "X");  
  
  h1->Draw();
  h2->Draw("sames");
  printStats(h1,h2);
  return;
}
