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
#include "TPaveText.h"
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

void printStats(TH1F *h1, TH1F *h2);
void plot(TH1F *h1, TH1F *h2);

//==========================================================//

void drawGlobalESDHistograms(const char* filename="HLT-OFFLINE-GLOBAL-comparison.root"){
 
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat("emr");
 gStyle->SetTitleX(gStyle->GetPadLeftMargin());

 TFile *f1 = TFile::Open(filename); 
 if(!f1 || f1->IsZombie()) {
    printf("file %s does not exist or there is an error opening it\n", filename);
    return;
 }

 TList *l1 = (TList*)f1->Get("global_histograms");
 if(!l1){
    printf("No list %s contained in your input file\n", l1->GetName()); 
    return; 
 }
 
 TCanvas *c1 = new TCanvas("c1","HLT vs. offline",1200,700);
 c1->Divide(3,3);
 TH1F *h1 = NULL;
 TH1F *h2 = NULL;

 h1 = (TH1F*)l1->FindObject("fNcluster_hlt"); if(!h1) { printf("Empty histogram fNcluster_hlt\n"); return; }
 h2 = (TH1F*)l1->FindObject("fNcluster_off"); if(!h2) { printf("Empty histogram fNcluster_off\n"); return; }
 h1->SetTitle("TPC cluster distribution");
 h1->GetXaxis()->SetTitle("TPC clusters per track");

 TLegend *leg1 = new TLegend(0.7,0.6,0.88,0.77);
 leg1->SetFillColor(10);
 leg1->SetLineColor(10);
 leg1->AddEntry(h1,"HLT", "l");
 leg1->AddEntry(h2,"OFF", "l");

 c1->cd(1);
 plot(h1,h2);
 leg1->Draw("same");

//-------------------------------------------------

 h1 = (TH1F*)l1->FindObject("fDCA_hlt"); if(!h1) { printf("Empty histogram fDCA_hlt\n"); return; }
 h2 = (TH1F*)l1->FindObject("fDCA_off"); if(!h2) { printf("Empty histogram fDCA_off\n"); return; }
 h1->SetTitle("DCA between track and vertex on XY plane");
 h1->SetXTitle("DCAr (cm)");
 
 c1->cd(2);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)l1->FindObject("fMult_hlt"); if(!h1) { printf("Empty histogram fMult_hlt\n"); return; }
 h2 = (TH1F*)l1->FindObject("fMult_off"); if(!h2) { printf("Empty histogram fMult_off\n"); return; }
 h1->SetTitle("track multiplicity");

 c1->cd(3);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)l1->FindObject("fCharge_hlt"); if(!h1) { printf("Empty histogram fCharge_hlt\n"); return; }
 h2 = (TH1F*)l1->FindObject("fCharge_off"); if(!h2) { printf("Empty histogram fCharge_off\n"); return; }
 h1->SetXTitle("polarity"); 
 h1->SetTitle("charge distribution");

 c1->cd(4);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)l1->FindObject("fMomentum_hlt"); if(!h1) { printf("Empty histogram fMomentum_hlt\n"); return; }
 h2 = (TH1F*)l1->FindObject("fMomentum_off"); if(!h2) { printf("Empty histogram fMomentum_off\n"); return; }
 h1->SetXTitle("p_{t} (GeV/c)"); 
 h1->SetTitle("transverse momentum");

 c1->cd(5);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)l1->FindObject("fEta_hlt"); if(!h1) { printf("Empty histogram fEta_hlt\n"); return; }
 h2 = (TH1F*)l1->FindObject("fEta_off"); if(!h2) { printf("Empty histogram fEta_off\n"); return; }
 h1->SetTitle("pseudorapidity");
 h1->SetXTitle("#eta");

 c1->cd(6);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)l1->FindObject("fXvertex_hlt"); if(!h1) { printf("Empty histogram fXvertex_hlt\n"); return; }
 h2 = (TH1F*)l1->FindObject("fXvertex_off"); if(!h2) { printf("Empty histogram fXvertex_off\n"); return; }
 h1->SetXTitle("x (cm)");
 h1->SetTitle("x of primary vertex");

 c1->cd(7);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)l1->FindObject("fYvertex_hlt"); if(!h1) { printf("Empty histogram fYvertex_hlt\n"); return; }
 h2 = (TH1F*)l1->FindObject("fYvertex_off"); if(!h2) { printf("Empty histogram fYvertex_off\n"); return; }
 h1->SetXTitle("y (cm)");
 h1->SetTitle("y of primary vertex");
 
 c1->cd(8);
 plot(h1,h2);

//------------------------------------------------- 

 h1 = (TH1F*)l1->FindObject("fZvertex_hlt"); if(!h1) { printf("Empty histogram fZvertex_hlt\n"); return; }
 h2 = (TH1F*)l1->FindObject("fZvertex_off"); if(!h2) { printf("Empty histogram fZvertex_off\n"); return; }
 h1->SetXTitle("z (cm)");
 h1->SetTitle("z of primary vertex");

 c1->cd(9);
 plot(h1,h2);

//------------------------------------------------- 

 c1->SaveAs("HLT-offline.png");  
 c1->SaveAs("HLT-offline.root");  
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
  st2->SetFillStyle(0);
  st2->Draw();  
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
  printStats(h1,h2);
  
  h1->Draw();
  h2->Draw("sames");
  return;
}
