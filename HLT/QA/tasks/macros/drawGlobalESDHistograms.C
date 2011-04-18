// $Id$
/*
 * Plotting macro for comparing offline- and HLT- ESD trees from
 * HLT-OFFLINE-GLOBAL-comparison.root produced using $ALICE_ROOT/HLT/QA/tasks/AliAnalysisTaskHLT.*
 *
 * Usage: aliroot drawGlobalESDHistograms.C'("HLT-OFFLINE-GLOBAL-comparison.root")'
 *
 * or aliroot drawGlobalESDHistograms.C++ in compiled mode
 *
 * It saves all canvases with the output histograms in a png and a ROOT file.
 * The second argument of the macro will produce individual files for all pads,
 * in case it is turned to kTRUE.
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
#include "TPaveText.h"
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

void drawGlobalESDHistograms(const char* filename="HLT-OFFLINE-GLOBAL-comparison.root", bool option=kFALSE){
 
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
 
 TText *runInfo = (TText*)list->FindObject("text");
 if(!runInfo) printf("No runInfo string\n");

 TText *cuts = (TText*)list->FindObject("cuts");
 if(!cuts) printf("No cuts string\n");

 TString folder = "GlobalTask_";
 folder += runInfo->GetTitle();
 folder.ReplaceAll(" ",""); 
 folder.ReplaceAll(",","_");
 gSystem->Exec("mkdir "+folder); // create a folder whose name contains run number and date of run

 const static int sizeTrack=8;
 const char *trackHLT[sizeTrack]={"fMomentum_hlt","fNcluster_hlt","fEta_hlt","fPhi_hlt","fDCAr_hlt","fDCAz_hlt","fCharge_hlt","fNITScluster_hlt"};
 const char *trackOFF[sizeTrack]={"fMomentum_off","fNcluster_off","fEta_off","fPhi_off","fDCAr_off","fDCAz_off","fCharge_off","fNITScluster_off"};
 const char *trackHLTcut[sizeTrack]={"fMomentum_hltcut","fNcluster_hltcut","fEta_hltcut","fPhi_hltcut","fDCAr_hltcut","fDCAz_hltcut","fCharge_hltcut","fNITScluster_hltcut"};

 const static int sizeEvent=9;
 const char *eventHLT[sizeEvent]={"fXvertex_hlt","fYvertex_hlt","fZvertex_hlt","fSPDXvertex_hlt","fSPDYvertex_hlt","fSPDZvertex_hlt","fMult_hlt","fNcont_hlt","fV0cent"};
 const char *eventOFF[sizeEvent]={"fXvertex_off","fYvertex_off","fZvertex_off","fSPDXvertex_off","fSPDYvertex_off","fSPDZvertex_off","fMult_off","fNcont_off","fV0cent"};

 TCanvas *c1 = new TCanvas("c1","track properties HLT vs. OFF",1250,700);
 c1->Divide(4,2);

 TH1F *h1 = NULL;
 TH1F *h2 = NULL;
 
 for(int i=0;i<sizeTrack;i++){
     c1->cd(i+1);
     h1 = (TH1F*)list->FindObject(trackHLT[i]); if(!h1) { return; }
     h2 = (TH1F*)list->FindObject(trackOFF[i]); if(!h2) { return; }
     plot(h1,h2);
     if(i==0){
        TLegend *leg1 = new TLegend(0.6,0.2,0.8,0.5);
        leg1->SetFillColor(10);
        leg1->SetLineColor(10);
        leg1->AddEntry(h1,"HLT", "l");
        leg1->AddEntry(h2,"OFF", "l");
        leg1->Draw("same");
     }
 }

 TCanvas *c4 = new TCanvas("c4","HLT track properties with and w/o cuts",1250,700);
 c4->Divide(4,2);
 
 for(int i=0;i<sizeTrack;i++){
     c4->cd(i+1);
     h1 = (TH1F*)list->FindObject(trackHLT[i]);    if(!h1) { return; }
     h2 = (TH1F*)list->FindObject(trackHLTcut[i]); if(!h2) { return; }
     plot(h1,h2);
     if(i==0){
  	TPaveText *pave = new TPaveText(2.1,24000,8.3,31600);
  	pave->SetFillColor(kWhite);
  	pave->SetLineColor(kWhite);
  	pave->SetShadowColor(kWhite);
  	TString tmp=cuts->GetTitle();
  	pave->SetTextColor(2);
  	pave->AddText(tmp);
  	pave->SetTextFont(42);
  	pave->SetTextSize(0.04);
  	pave->Draw();
  	c4->Update();
     }
 }

 TCanvas *c2 = new TCanvas("c2","vertex event properties",1200,700);
 c2->Divide(3,2);
 
 for(int i=0;i<6;i++){
     c2->cd(i+1);
     h1 = (TH1F*)list->FindObject(eventHLT[i]); if(!h1) { return; }
     h2 = (TH1F*)list->FindObject(eventOFF[i]); if(!h2) { return; }
     plot(h1,h2);
     if(i==0){
        TLegend *leg1 = new TLegend(0.6,0.2,0.8,0.5);
        leg1->SetFillColor(10);
        leg1->SetLineColor(10);
        leg1->AddEntry(h1,"HLT", "l");
        leg1->AddEntry(h2,"OFF", "l");
        leg1->Draw("same");
     }
 }

 TCanvas *c3 = new TCanvas("c3","general event properties",1200,500);
 c3->Divide(3,1);

 for(int i=6;i<9;i++){
     c3->cd(i-5);
     h2 = (TH1F*)list->FindObject(eventOFF[i]); if(!h2) { return; }
     h1 = (TH1F*)list->FindObject(eventHLT[i]); if(!h1) { return; }
     plot(h1,h2);
     if(i==6){
        TLegend *leg1 = new TLegend(0.6,0.2,0.8,0.5);
        leg1->SetFillColor(10);
        leg1->SetLineColor(10);
        leg1->AddEntry(h1,"HLT", "l");
        leg1->AddEntry(h2,"OFF", "l");
        leg1->Draw("same");
     }
 }
 
 c1->SaveAs(folder+"/track_properties.png");  
 c1->SaveAs(folder+"/track_properties.root");  
 c2->SaveAs(folder+"/vertex_event_properties.png");  
 c2->SaveAs(folder+"/vertex_event_properties.root");  
 c3->SaveAs(folder+"/general_event_properties.png");  
 c3->SaveAs(folder+"/general_event_properties.root");  
 c4->SaveAs(folder+"/HLT_track_properties_cuts.png");  
 c4->SaveAs(folder+"/HLT_track_properties_cuts.root");  
 
 if(option==kTRUE){ 
    TPad *pad = NULL; 
    for(int i=1; i<=sizeTrack; i++){
       pad = (TPad*)c1->GetListOfPrimitives()->FindObject(Form("c1_%d",i));
       if(!pad){
  	  printf("Empty pad %d in canvas %s.\n", i, c1->GetName());
  	  continue;	    
       }
       pad->SaveAs(Form(folder+"/c1_%s_off.png",trackHLT[i-1]));
    }

    for(int i=1; i<=sizeTrack; i++){
       pad = (TPad*)c4->GetListOfPrimitives()->FindObject(Form("c4_%d",i));
       if(!pad){
  	  printf("Empty pad %d in canvas %s.\n", i, c4->GetName());
  	  continue;	    
       }
       pad->SaveAs(Form(folder+"/c4_%s_cuts.png",trackHLT[i-1]));
    }
    
    for(int i=1; i<7; i++){
       pad = (TPad*)c2->GetListOfPrimitives()->FindObject(Form("c2_%d",i));
       if(!pad){
  	  printf("Empty pad %d in canvas %s.\n", i, c2->GetName());
  	  continue;	    
       }
       pad->SaveAs(Form(folder+"/c2_%s_off.png",eventHLT[i-1]));
    }

    for(int i=6; i<9; i++){
       pad = (TPad*)c3->GetListOfPrimitives()->FindObject(Form("c3_%d",i-5));
       if(!pad){
  	  printf("Empty pad %d in canvas %s.\n", i-5, c3->GetName());
  	  continue;	    
       }
       pad->SaveAs(Form(folder+"/c3_%s_off.png",eventHLT[i]));
    }
 }
 
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
