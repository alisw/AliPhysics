//-*- Mode: C++ -*-

// Make Hists for trigger efficiency studies
// Author: Jochen Thaeder <jochen@thaeder.de> 

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TLegend.h"
#include "TH1.h"
#include "TString.h"
#include "TMath.h"
#include "TLine.h"

// ----------------------------------------------------------------------
// --                      Static Variablen                            --
// ----------------------------------------------------------------------
TObjArray* fObjArray = NULL;

Int_t fgImgIdx = 0;

const Int_t fgkNTrigger = 10;

const Double_t fgkTriggerPt[fgkNTrigger] = {1.0, 2.0, 2.5, 3.0, 5.0, 7.0, 10.0, 1000., 1001., 1002.};

Char_t *fgkTrigger[fgkNTrigger] = {
  "p_{t}> 1.0",              "p_{t}> 2.0",              "p_{t}> 2.5",
  "p_{t}> 3.0",              "p_{t}> 5.0",              "p_{t}> 7",
  "p_{t}> 10.",
  "S1",                      "S2",                      "S3"
};

const Int_t fgkNSelectionCuts  = 5;
Char_t *fgkSelectionCuts[fgkNSelectionCuts] = { 
  "All Events", 
  "AliPhysSel", 
  "AliPhysSel - PrimVertex",
  "AliPhysSel - PrimVertex - Track (OFF)",
  "AliPhysSel - PrimVertex - Track (HLT)",
};

Char_t *fgkSelectionCutsMC[fgkNSelectionCuts] = { 
  "AllEvents", 
  "AliPhysSel", 
  "AliPhysSel - PrimVertex - Charged Primary",
  "AliPhysSel - PrimVertex - Track (MC)",
  "AliPhysSel - PrimVertex - Track (MC)",
};

Char_t *fgkSelectionCutsRatioMC[fgkNSelectionCuts] = { 
  "AllEvents", 
  "AliPhysSel", 
  "AliPhysSel - PrimVertex - Charged Primary",
  "AliPhysSel - PrimVertex - Track (OFF)/ Track (MC)",
  "AliPhysSel - PrimVertex - Track (HLT)/ Track (MC)",
};

// ----------------------------------------------------------------------
// --                       Main Function                              --
// ----------------------------------------------------------------------

void readTrigger(const Char_t* folder = "../..");

// ----------------------------------------------------------------------
// --                      Canvas Functions                            --
// ----------------------------------------------------------------------

void CanvasTrigger();
void CanvasTriggerFactors();
void CanvasTriggerFactorsBackup();

void CreateCanvasCuts();

void CreateCanvasDistributions( Int_t maxLayer = 99, Int_t draw = 2, Int_t data = 2, Int_t minLayer=0);
// draw 0 : pt  || draw 1 : mult || draw 2 : both
// data 0 : OFF || data 1 : MC   || data 2 : both

// ----------------------------------------------------------------------
// --                      Helper Functions                            --
// ----------------------------------------------------------------------

void SetHist(TH1F* hist, Int_t hType, Char_t* dataType, Char_t* histType, Float_t minY, Float_t maxY );

void DivideHist(TH1F* h0, TH1F* h1, TH1F* h2);
void FillReduxHistograms( TH1F* hN, TH1F* hRedux);
void FillReduxHistogramsPt( TH1F* hN, TH1F* hRedux, TH1F* hReduxW);

void FillRatioHistograms( TH1F* hN, TH1F* hF, TH1F* hRatio);
void FillRatioHistogramsPt( TH1F* hN, TH1F* hF, TH1F* hRatio);

// ----------------------------------------------------------------------
// --                       Draw Functions                             --
// ----------------------------------------------------------------------

void DrawHistogramTrigger( TCanvas* canvas, Int_t idx, TH1* hist, Int_t iLayer, Bool_t bLogY);
void DrawHistogram( TCanvas* canvas, Int_t idx, TH1* hist, Int_t iLayer, 
		    Bool_t bScale, Bool_t bLogY, Bool_t bLogX);

//#######################################################################
void readTrigger(const Char_t* folder) {

  gStyle->SetPalette(1);
  gROOT->SetStyle("Plain");

  TString path(folder);
  path += "/jthaeder_trigger.root";
  
  TFile *f = new TFile(path.Data());

  // ----------------------------------------------------------------------------

  fObjArray = static_cast<TObjArray*>(f->Get("jthaeder_trigger"));
  
  // == SHOW ===================
  CanvasTriggerFactors();
  CreateCanvasDistributions(  7, 0, 2, 0 );

  // CreateCanvasDistributions( 10, 1, 0, 7 );
  // == SHOW ===================

  // == BACKUP =================
#if 1
  CanvasTriggerFactorsBackup();
  CanvasTrigger();
  CreateCanvasCuts();
#endif
  // == BACKUP =================
}

// ----------------------------------------------------------------------
// --                      Canvas Functions                            --
// ----------------------------------------------------------------------

//#######################################################################
void CanvasTrigger() {
  
  TCanvas *cTrigger = new TCanvas("cTrigger", "Trigger Events", 10, 10, 1400, 800);
  cTrigger->Divide(2,3);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TLegend* leg = new TLegend(0.60, 0.12, 0.89, 0.37, "Event Selection");
  leg->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggered")),    "OFF Trig - All Events","LP");
  leg->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggeredSel")), "OFF Trig - AliPhysSel - OFF PrimVertex","LP");
  leg->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggered")),    "HLT Trig - All Events","LP");
  leg->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggeredSel")), "HLT Trig - AliPhysSel - HLT PrimVertex","LP");
  leg->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggered")),     "MC Trig - All Events","LP");
  leg->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggeredSel")),  "MC Trig - AliPhysSel - MC PrimVertex","LP");
  leg->SetFillColor(kWhite);

  TLegend* legW = new TLegend(0.60, 0.12, 0.89, 0.37, "Event Selection");
  legW->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggered")),    "OFF Trig - All Events","LP");
  legW->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggeredSel")), "OFF Trig - AliPhysSel - OFF PrimVertex","LP");
  legW->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggered")),    "HLT Trig - All Events","LP");
  legW->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggeredSel")), "HLT Trig - AliPhysSel - HLT PrimVertex","LP");
  legW->SetFillColor(kWhite);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  
  TH1F *hOFF = static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggered"));

  hOFF->SetTitle("Triggered Events");
  hOFF->SetAxisRange(0,fgkNTrigger);
  hOFF->SetMinimum(1);

  hOFF->GetXaxis()->SetBinLabel(1, "Total N events");  
  hOFF->GetXaxis()->SetBinLabel(fgkNTrigger+2, "Total N events");  
  for ( Int_t ii=2 ; ii <= fgkNTrigger+1 ; ii++ ) {
    hOFF->GetXaxis()->SetBinLabel(ii, fgkTrigger[ii-2]);
    hOFF->GetXaxis()->SetBinLabel(ii+fgkNTrigger+1, fgkTrigger[ii-2]);
  }

  DrawHistogramTrigger(cTrigger, 1, static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggered")),    0, kTRUE);
  DrawHistogramTrigger(cTrigger, 1, static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggeredSel")), 1, kTRUE);
  DrawHistogramTrigger(cTrigger, 1, static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggered")),    2, kTRUE);
  DrawHistogramTrigger(cTrigger, 1, static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggeredSel")), 3, kTRUE);
  DrawHistogramTrigger(cTrigger, 1, static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggered")),     4, kTRUE);
  DrawHistogramTrigger(cTrigger, 1, static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggeredSel")),  5, kTRUE);
  leg->Draw();

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TH1F* h1Clone = static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggered")->Clone());
  h1Clone->SetAxisRange(fgkNTrigger+1,2*(fgkNTrigger+1));
  h1Clone->GetYaxis()->SetTitle("N Events * N Tracks");

  TString title(h1Clone->GetTitle());
  title += " - weighted nTracks";
  h1Clone->SetTitle(title);    

  DrawHistogramTrigger(cTrigger, 2, h1Clone, 0, kTRUE);
  DrawHistogramTrigger(cTrigger, 2, static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggeredSel")), 1, kTRUE);
  DrawHistogramTrigger(cTrigger, 2, static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggered")),    2, kTRUE);
  DrawHistogramTrigger(cTrigger, 2, static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggeredSel")), 3, kTRUE);
  legW->Draw();

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  // ---------------------------------------------------------------
  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TLegend* legF = new TLegend(0.65, 0.12, 0.89, 0.32, "Event Selection");
  legF->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToOFF")),
		 "HLT Trig vs OFF - All Events","LP");
  legF->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToOFFSel")),
		 "HLT Trig vs OFF - AliPhysSel - OFF PrimVertex","LP");

  legF->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToMC")),
		 "HLT Trig vs MC - All Events","LP");
  legF->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToMCSel")),
		 "HLT Trig vs MC - AliPhysSel - OFF PrimVertex","LP");

  legF->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNOFFFakeToMC")),
		 "OFF Trig vs MC - AllEvents","LP");
  legF->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNOFFFakeToMCSel")),
		 "OFF Trig vs MC - AliPhysSel - OFF PrimVertex","LP");

  legF->SetFillColor(kWhite);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  
  TH1F* hFake = static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToOFF"));
  hFake->SetAxisRange(1,fgkNTrigger);
  hFake->SetAxisRange(1, 1e5,"Y");
  hFake->GetYaxis()->SetTitle("N Events #left(HLT_{Triggered} && ! [OFF|MC]_{Triggered}#right)");

  for ( Int_t ii=2 ; ii <= fgkNTrigger+1 ; ii++ ) {
    hFake->GetXaxis()->SetBinLabel(ii, fgkTrigger[ii-2]);
    hFake->GetXaxis()->SetBinLabel(ii+fgkNTrigger+1, fgkTrigger[ii-2]);
  }

  DrawHistogramTrigger(cTrigger, 3, hFake, 0, kTRUE);
  DrawHistogramTrigger(cTrigger, 3, static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToOFFSel")), 1, kTRUE);
  DrawHistogramTrigger(cTrigger, 3, static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToMC")),     2, kTRUE);
  DrawHistogramTrigger(cTrigger, 3, static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToMCSel")),  3, kTRUE);
  DrawHistogramTrigger(cTrigger, 3, static_cast<TH1F*>(fObjArray->FindObject("fNOFFFakeToMC")),     4, kTRUE);
  DrawHistogramTrigger(cTrigger, 3, static_cast<TH1F*>(fObjArray->FindObject("fNOFFFakeToMCSel")),  5, kTRUE);
  legF->Draw();

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TH1F* h2Clone = static_cast<TH1F*>(hFake->Clone());
  h2Clone->SetAxisRange(fgkNTrigger+2,2*(fgkNTrigger+1));
  h2Clone->SetAxisRange(1, 1e6, "Y");
  h2Clone->GetYaxis()->SetTitle("N Events * N Tracks");

  TString titleF(h2Clone->GetTitle());
  titleF += " - weighted nTracks";
  h2Clone->SetTitle(titleF);    

  DrawHistogramTrigger(cTrigger, 4, h2Clone, 0, kTRUE);
  DrawHistogramTrigger(cTrigger, 4, static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToOFFSel")), 1, kTRUE);
  DrawHistogramTrigger(cTrigger, 4, static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToMC")),     2, kTRUE);
  DrawHistogramTrigger(cTrigger, 4, static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToMCSel")),  3, kTRUE);
  DrawHistogramTrigger(cTrigger, 4, static_cast<TH1F*>(fObjArray->FindObject("fNOFFFakeToMC")),     4, kTRUE);
  DrawHistogramTrigger(cTrigger, 4, static_cast<TH1F*>(fObjArray->FindObject("fNOFFFakeToMCSel")),  5, kTRUE);

  legF->Draw();

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  // ---------------------------------------------------------------
  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TLegend* legM = new TLegend(0.65, 0.12, 0.89, 0.32, "Event Selection");
  legM->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToOFF")),
		 "HLT Trig vs OFF - All Events","LP");
  legM->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToOFFSel")),
		 "HLT Trig vs OFF - AliPhysSel - OFF PrimVertex","LP");

  legM->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToMC")),
		 "HLT Trig vs MC - All Events","LP");
  legM->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToMCSel")),
		 "HLT Trig vs MC - AliPhysSel - OFF PrimVertex","LP");

  legM->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNOFFMissToMC")),
		 "OFF Trig vs MC - All Events","LP");
  legM->AddEntry(static_cast<TH1F*>(fObjArray->FindObject("fNOFFMissToMCSel")),
		 "OFF Trig vs MC - AliPhysSel - OFF PrimVertex","LP");

  legM->SetFillColor(kWhite);

  // --- --- --- --- --- --- -

  TH1F* hMiss = static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToOFF"));
  hMiss->SetAxisRange(1,fgkNTrigger);
  hMiss->SetAxisRange(1, 1e4, "Y");
  hMiss->GetYaxis()->SetTitle("N Events #left(!HLT_{Triggered} && [OFF|MC]_{Triggered}#right)");

  for ( Int_t ii=2 ; ii <= fgkNTrigger+1 ; ii++ ) {
    hMiss->GetXaxis()->SetBinLabel(ii, fgkTrigger[ii-2]);
    hMiss->GetXaxis()->SetBinLabel(ii+fgkNTrigger+1, fgkTrigger[ii-2]);
  }

  DrawHistogramTrigger(cTrigger, 5, hMiss, 0, kTRUE);
  DrawHistogramTrigger(cTrigger, 5, static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToOFFSel")), 1, kTRUE);
  DrawHistogramTrigger(cTrigger, 5, static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToMC")),     2, kTRUE);
  DrawHistogramTrigger(cTrigger, 5, static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToMCSel")),  3, kTRUE);
  DrawHistogramTrigger(cTrigger, 5, static_cast<TH1F*>(fObjArray->FindObject("fNOFFMissToMC")),     4, kTRUE);
  DrawHistogramTrigger(cTrigger, 5, static_cast<TH1F*>(fObjArray->FindObject("fNOFFMissToMCSel")),  5, kTRUE);

  legM->Draw();

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TH1F* h3Clone = static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToOFF")->Clone());
  h3Clone->SetAxisRange(fgkNTrigger+2,2*(fgkNTrigger+1));
  h3Clone->SetAxisRange(1, 1e6, "Y");
  h3Clone->GetYaxis()->SetTitle("N Events * N Tracks");

  TString titleM(h3Clone->GetTitle());
  titleM += " - weighted nTracks";
  h3Clone->SetTitle(titleM);    

  DrawHistogramTrigger(cTrigger, 6, h3Clone, 0, kTRUE);
  DrawHistogramTrigger(cTrigger, 6, static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToOFFSel")), 1, kTRUE);
  DrawHistogramTrigger(cTrigger, 6, static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToMC")),     2, kTRUE);
  DrawHistogramTrigger(cTrigger, 6, static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToMCSel")),  3, kTRUE);
  DrawHistogramTrigger(cTrigger, 6, static_cast<TH1F*>(fObjArray->FindObject("fNOFFMissToMC")),     4, kTRUE);
  DrawHistogramTrigger(cTrigger, 6, static_cast<TH1F*>(fObjArray->FindObject("fNOFFMissToMCSel")),  5, kTRUE);

  legM->Draw();

  cTrigger->SaveAs(Form("TriggerHists_%d_B.png",++fgImgIdx));
}
// === === === === === === === === === === === === === === === ===
// ===============================================================
// === === === === === === === === === === === === === === === ===

//#######################################################################
void CanvasTriggerFactorsBackup() {

  TCanvas *cTriggerFB = new TCanvas("cTriggerFB", "Reduction - Purity - Efficiency", 10, 10, 1400, 800);
  cTriggerFB->Divide(3,3);
 
  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  Float_t maxLine = 0.9;
  Float_t minLine = 0.8;

  TLine* line1 = new TLine(0.,maxLine,fgkNTrigger,maxLine);
  line1->SetLineStyle(2);
  line1->SetLineColor(kGreen);

  TLine* line2 = new TLine(0.,minLine,fgkNTrigger,minLine);
  line2->SetLineStyle(2);
  line2->SetLineColor(kRed);


  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  TH1F* hOFFReduxB   = new TH1F("hOFFReduxB", 
			       "Reduction Factors;;Reduction #left(#frac{N Events_{Triggered}}{N Events_{Total}}#right)", 
			       2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));
  TH1F* hOFFReduxBS  = new TH1F("hOFFReduxBS", "",  2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));
  TH1F* hHLTReduxB   = new TH1F("hHLTReduxB", "",  2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));
  TH1F* hHLTReduxBS  = new TH1F("hHLTReduxBS", "", 2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));
  TH1F* hMCReduxB    = new TH1F("hMCReduxB", "",  2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));
  TH1F* hMCReduxBS   = new TH1F("hMCReduxBS", "", 2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));

  FillReduxHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggered")),    hOFFReduxB);
  FillReduxHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggeredSel")), hOFFReduxBS);
  FillReduxHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggered")),    hHLTReduxB);
  FillReduxHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggeredSel")), hHLTReduxBS);
  FillReduxHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggered")),     hMCReduxB);
  FillReduxHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggeredSel")),  hMCReduxBS);

  TLegend* legR = new TLegend(0.65, 0.12, 0.89, 0.37, "Event Selection");
  legR->AddEntry(hOFFReduxB,  "OFF Trig - All Events","LP");
  legR->AddEntry(hOFFReduxBS, "OFF Trig - AliPhysSel - OFF PrimVertex","LP");
  legR->AddEntry(hHLTReduxB, "HLT Trig - All Events","LP");
  legR->AddEntry(hHLTReduxBS,"HLT Trig - AliPhysSel - OFF PrimVertex","LP");
  legR->AddEntry(hMCReduxB, "MC Trig - All Events","LP");
  legR->AddEntry(hMCReduxBS,"MC Trig - AliPhysSel - OFF PrimVertex","LP");
  legR->SetFillColor(kWhite);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  hOFFReduxB->SetAxisRange(0,fgkNTrigger-1);

  for ( Int_t ii=1 ; ii <= fgkNTrigger ; ii++ ) {
    hOFFReduxB->GetXaxis()->SetBinLabel(ii, fgkTrigger[ii-1]);
    hOFFReduxB->GetXaxis()->SetBinLabel(ii+fgkNTrigger, fgkTrigger[ii-1]);
  }

  DrawHistogramTrigger(cTriggerFB, 1, hOFFReduxB,  0, kTRUE);
  DrawHistogramTrigger(cTriggerFB, 1, hOFFReduxBS, 1, kTRUE);
  DrawHistogramTrigger(cTriggerFB, 1, hHLTReduxB,  2, kTRUE);
  DrawHistogramTrigger(cTriggerFB, 1, hHLTReduxBS, 3, kTRUE);
  DrawHistogramTrigger(cTriggerFB, 1, hMCReduxB,   4, kTRUE);
  DrawHistogramTrigger(cTriggerFB, 1, hMCReduxBS,  5, kTRUE);
  legR->Draw();

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TH1F* hOFFRClone = static_cast<TH1F*>(hOFFReduxB->Clone());
  hOFFRClone->SetAxisRange(fgkNTrigger,2*(fgkNTrigger));

  TString titleR(hOFFRClone->GetTitle());
  titleR += " - weighted nTracks";
  hOFFRClone->SetTitle(titleR);    

  DrawHistogramTrigger(cTriggerFB, 2, hOFFRClone, 0, kTRUE);
  DrawHistogramTrigger(cTriggerFB, 2, hOFFReduxBS, 1, kTRUE);
  DrawHistogramTrigger(cTriggerFB, 2, hHLTReduxB,  2, kTRUE);
  DrawHistogramTrigger(cTriggerFB, 2, hHLTReduxBS, 3, kTRUE);
  DrawHistogramTrigger(cTriggerFB, 2, hMCReduxB,   4, kTRUE);
  DrawHistogramTrigger(cTriggerFB, 2, hMCReduxBS,  5, kTRUE);
  legR->Draw();

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  // ---------------------------------------------------------------
  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TH1F* hOFFPurityHLTB  = new TH1F("hOFFPurityHLTB", "HLT Trigger Purity vs OFF", 2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));
  hOFFPurityHLTB->GetYaxis()->SetTitle("Purity #left(#frac{HLT_{Triggered} - HLT_{Fake}}{HLT_{Triggered}}#right)");
  TH1F* hOFFPurityHLTBS = new TH1F("hOFFPurityHLTBS", "", 2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));
  
  TH1F* hMCPurityHLTB  = new TH1F("hMCPurityHLTB", "HLT Trigger Purity vs MC", 2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));
  hMCPurityHLTB->GetYaxis()->SetTitle("Purity #left(#frac{HLT_{Triggered} - HLT_{Fake}}{HLT_{Triggered}}#right)");
  TH1F* hMCPurityHLTBS = new TH1F("hMCPurityHLTBS", "", 2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));

  TH1F* hMCPurityOFFB  = new TH1F("hMCPurityOFFB", "OFF Trigger Purity vs MC", 2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));
  hMCPurityOFFB->GetYaxis()->SetTitle("Purity #left(#frac{OFF_{Triggered} - OFF_{Fake}}{OFF_{Triggered}}#right)");
  TH1F* hMCPurityOFFBS = new TH1F("hMCPurityOFFBS", "", 2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));

  FillRatioHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggered")),    
		      static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToOFF")), hOFFPurityHLTB);
  FillRatioHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggeredSel")),
		      static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToOFFSel")), hOFFPurityHLTBS);

  FillRatioHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggered")),    
		      static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToMC")), hMCPurityHLTB);
  FillRatioHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggeredSel")),
		      static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToMCSel")), hMCPurityHLTBS);

  FillRatioHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggered")),    
		      static_cast<TH1F*>(fObjArray->FindObject("fNOFFFakeToMC")), hMCPurityOFFB);
  FillRatioHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggeredSel")),
		      static_cast<TH1F*>(fObjArray->FindObject("fNOFFFakeToMCSel")), hMCPurityOFFBS);

  TLegend* legPH = new TLegend(0.65, 0.72, 0.89, 0.87, "Event Selection");
  legPH->AddEntry(hOFFPurityHLTB, "HLT Trig - All Events","LP");
  legPH->AddEntry(hOFFPurityHLTBS,"HLT Trig - AliPhysSel - OFF PrimVertex","LP");
  legPH->SetFillColor(kWhite);

  TLegend* legPO = new TLegend(0.65, 0.72, 0.89, 0.87, "Event Selection");
  legPO->AddEntry(hMCPurityOFFB, "OFF Trig - All Events","LP");
  legPO->AddEntry(hMCPurityOFFBS,"OFF Trig - AliPhysSel - OFF PrimVertex","LP");
  legPO->SetFillColor(kWhite);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  hOFFPurityHLTB->SetAxisRange(0,fgkNTrigger-1);
  hOFFPurityHLTB->SetAxisRange(0., 1.2, "Y");

  for ( Int_t ii=1 ; ii <= fgkNTrigger ; ii++ ){
    hOFFPurityHLTB->GetXaxis()->SetBinLabel(ii, fgkTrigger[ii-1]);
    hOFFPurityHLTB->GetXaxis()->SetBinLabel(ii+fgkNTrigger, fgkTrigger[ii-1]);
  }

  DrawHistogramTrigger(cTriggerFB, 4, hOFFPurityHLTB,  0, kFALSE);
  DrawHistogramTrigger(cTriggerFB, 4, hOFFPurityHLTBS, 1, kFALSE);
  line1->Draw();
  line2->Draw();
  legPH->Draw(); 

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  hMCPurityHLTB->SetAxisRange(0,fgkNTrigger-1);
  hMCPurityHLTB->SetAxisRange(0., 1.2, "Y");

  for ( Int_t ii=1 ; ii <= fgkNTrigger ; ii++ ){
    hMCPurityHLTB->GetXaxis()->SetBinLabel(ii, fgkTrigger[ii-1]);
    hMCPurityHLTB->GetXaxis()->SetBinLabel(ii+fgkNTrigger, fgkTrigger[ii-1]);
  }

  DrawHistogramTrigger(cTriggerFB, 5, hMCPurityHLTB,  0, kFALSE);
  DrawHistogramTrigger(cTriggerFB, 5, hMCPurityHLTBS, 1, kFALSE);
  line1->Draw();
  line2->Draw();
  legPH->Draw(); 

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  hMCPurityOFFB->SetAxisRange(0,fgkNTrigger-1);
  hMCPurityOFFB->SetAxisRange(0., 1.2, "Y");

  for ( Int_t ii=1 ; ii <= fgkNTrigger ; ii++ ){
    hMCPurityOFFB->GetXaxis()->SetBinLabel(ii, fgkTrigger[ii-1]);
    hMCPurityOFFB->GetXaxis()->SetBinLabel(ii+fgkNTrigger, fgkTrigger[ii-1]);
  }

  DrawHistogramTrigger(cTriggerFB, 6, hMCPurityOFFB,  0, kFALSE);
  DrawHistogramTrigger(cTriggerFB, 6, hMCPurityOFFBS, 1, kFALSE);
  line1->Draw();
  line2->Draw();
  legPH->Draw(); 

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  // ---------------------------------------------------------------
  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TH1F* hOFFEffHLTB  = new TH1F("hOFFEffHLTB", "HLT Trigger Efficiency vs OFF", 2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));
  hOFFEffHLTB->GetYaxis()->SetTitle("Efficiency #left(#frac{OFF_{Triggered} - HLT_{Miss}}{OFF_{Triggered}}#right)");
  TH1F* hOFFEffHLTBS = new TH1F("hOFFEffHLTBS", "", 2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));

  TH1F* hMCEffHLTB  = new TH1F("hMCEffHLTB", "HLT Trigger Efficiency vs MC", 2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));
  hMCEffHLTB->GetYaxis()->SetTitle("Efficiency #left(#frac{MC_{Triggered} - HLT_{Miss}}{MC_{Triggered}}#right)");
  TH1F* hMCEffHLTBS = new TH1F("hMCEffHLTBS", "", 2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));

  TH1F* hMCEffOFFB  = new TH1F("hMCEffOFFB", "OFF Trigger Efficiency vs MC", 2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));
  hMCEffOFFB->GetYaxis()->SetTitle("Efficiency #left(#frac{MC_{Triggered} - OFF_{Miss}}{MC_{Triggered}}#right)");
  TH1F* hMCEffOFFBS = new TH1F("hMCEffOFFBS", "", 2*fgkNTrigger, 0., static_cast<Double_t>(2*fgkNTrigger));

  
  FillRatioHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggered")),    
		      static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToOFF")), hOFFEffHLTB);
  FillRatioHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggeredSel")),
		      static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToOFFSel")), hOFFEffHLTBS);
 
  FillRatioHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggered")),    
		      static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToMC")), hMCEffHLTB);
  FillRatioHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggeredSel")),
		      static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToMCSel")), hMCEffHLTBS);

  FillRatioHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggered")),    
		      static_cast<TH1F*>(fObjArray->FindObject("fNOFFMissToMC")), hMCEffOFFB);
  FillRatioHistograms(static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggeredSel")),
		      static_cast<TH1F*>(fObjArray->FindObject("fNOFFMissToMCSel")), hMCEffOFFBS);


  TLegend* legEH = new TLegend(0.65, 0.70, 0.89, 0.85, "Event Selection");
  legEH->AddEntry(hOFFEffHLTB, "HLT Trig - All Events","LP");
  legEH->AddEntry(hOFFEffHLTBS,"HLT Trig - AliPhysSel - OFF PrimVertex","LP");
  legEH->SetFillColor(kWhite);

  TLegend* legEO = new TLegend(0.65, 0.70, 0.89, 0.85, "Event Selection");
  legEO->AddEntry(hMCEffOFFB, "OFF Trig - All Events","LP");
  legEO->AddEntry(hMCEffOFFBS,"OFF Trig - AliPhysSel - OFF PrimVertex","LP");
  legEO->SetFillColor(kWhite);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  hOFFEffHLTB->SetAxisRange(0,fgkNTrigger-1);
  hOFFEffHLTB->SetAxisRange(0., 1.2, "Y");
 
  for ( Int_t ii=1 ; ii <= fgkNTrigger ; ii++ ){
    hOFFEffHLTB->GetXaxis()->SetBinLabel(ii, fgkTrigger[ii-1]);
    hOFFEffHLTB->GetXaxis()->SetBinLabel(ii+fgkNTrigger, fgkTrigger[ii-1]);
  }

  DrawHistogramTrigger(cTriggerFB, 7, hOFFEffHLTB, 0, kFALSE);
  DrawHistogramTrigger(cTriggerFB, 7, hOFFEffHLTBS, 1, kFALSE);
  line1->Draw();
  line2->Draw();
  legEH->Draw(); 

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  hMCEffHLTB->SetAxisRange(0,fgkNTrigger-1);
  hMCEffHLTB->SetAxisRange(0., 1.2, "Y");

  for ( Int_t ii=1 ; ii <= fgkNTrigger ; ii++ ){
    hMCEffHLTB->GetXaxis()->SetBinLabel(ii, fgkTrigger[ii-1]);
    hMCEffHLTB->GetXaxis()->SetBinLabel(ii+fgkNTrigger, fgkTrigger[ii-1]);
  }

  DrawHistogramTrigger(cTriggerFB, 8, hMCEffHLTB,  0, kFALSE);
  DrawHistogramTrigger(cTriggerFB, 8, hMCEffHLTBS, 1, kFALSE);
  line1->Draw();
  line2->Draw();
  legEH->Draw(); 

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  hMCEffOFFB->SetAxisRange(0,fgkNTrigger-1);
  hMCEffOFFB->SetAxisRange(0., 1.2, "Y");

  for ( Int_t ii=1 ; ii <= fgkNTrigger ; ii++ ){
    hMCEffOFFB->GetXaxis()->SetBinLabel(ii, fgkTrigger[ii-1]);
    hMCEffOFFB->GetXaxis()->SetBinLabel(ii+fgkNTrigger, fgkTrigger[ii-1]);
  }

  DrawHistogramTrigger(cTriggerFB, 9, hMCEffOFFB,  0, kFALSE);
  DrawHistogramTrigger(cTriggerFB, 9, hMCEffOFFBS, 1, kFALSE);
  line1->Draw();
  line2->Draw();
  legEO->Draw(); 

  cTriggerFB->SaveAs(Form("TriggerHists_%d_B.png",++fgImgIdx));

  for ( Int_t ii=1 ; ii <= fgkNTrigger ; ii++ ){
    if (fgkTriggerPt[ii-1] > 20. )
      continue;

    printf("=================================\n");
    printf(" Trigger Pt %.2f GeV/c \n", fgkTriggerPt[ii-1]);
    printf("=================================\n");
    printf(" - HLT Redux Factor  : %.2f - %.2f Hz (1500 Hz) - %.2f Hz (800 Hz) \n", 
	   hHLTReduxB->GetBinContent(ii), 1500/hHLTReduxB->GetBinContent(ii), 800/hHLTReduxB->GetBinContent(ii) );
    printf(" - OFF Redux Factor  : %.2f - %.2f Hz (1500 Hz) - %.2f Hz (800 Hz) \n", 
	   hOFFReduxB->GetBinContent(ii), 1500/hOFFReduxB->GetBinContent(ii), 800/hOFFReduxB->GetBinContent(ii) );
    printf(" - MC  Redux Factor  : %.2f - %.2f Hz (1500 Hz) - %.2f Hz (800 Hz) \n", 
	   hMCReduxB->GetBinContent(ii), 1500/hMCReduxB->GetBinContent(ii), 800/hMCReduxB->GetBinContent(ii) );
    printf("---------------------------------\n");
    printf(" - HLT Purity vs OFF : %.3f | HLT Eff vs OFF : %.3f \n", hOFFPurityHLTB->GetBinContent(ii),hOFFEffHLTB->GetBinContent(ii));
    printf(" - HLT Purity vs MC  : %.3f | HLT Eff vs MC  : %.3f \n", hMCPurityHLTB->GetBinContent(ii), hMCEffHLTB->GetBinContent(ii) );
    printf(" - OFF Purity vs MC  : %.3f | OFF Eff vs MC  : %.3f \n", hMCPurityOFFB->GetBinContent(ii), hMCEffOFFB->GetBinContent(ii) );
  }
 
  for ( Int_t ii=1 ; ii <= fgkNTrigger ; ii++ ){
    if (fgkTriggerPt[ii-1] < 20. )
      continue;

    printf("=================================\n");
    printf(" Trigger Scenario %.0f \n", fgkTriggerPt[ii-1]-1000+1);
    printf("=================================\n");
    printf(" - HLT Redux Factor  : %.2f - %.2f Hz (1500 Hz) - %.2f Hz (800 Hz) \n", 
	   hHLTReduxB->GetBinContent(ii), 1500/hHLTReduxB->GetBinContent(ii), 800/hHLTReduxB->GetBinContent(ii) );
    printf(" - OFF Redux Factor  : %.2f - %.2f Hz (1500 Hz) - %.2f Hz (800 Hz) \n", 
	   hOFFReduxB->GetBinContent(ii), 1500/hOFFReduxB->GetBinContent(ii), 800/hOFFReduxB->GetBinContent(ii) );
    printf(" - MC  Redux Factor  : %.2f - %.2f Hz (1500 Hz) - %.2f Hz (800 Hz) \n", 
	   hMCReduxB->GetBinContent(ii), 1500/hMCReduxB->GetBinContent(ii), 800/hMCReduxB->GetBinContent(ii) );
    printf("---------------------------------\n");
    printf(" - HLT Purity vs OFF : %.3f | HLT Eff vs OFF : %.3f \n", hOFFPurityHLTB->GetBinContent(ii),hOFFEffHLTB->GetBinContent(ii));
    printf(" - HLT Purity vs MC  : %.3f | HLT Eff vs MC  : %.3f \n", hMCPurityHLTB->GetBinContent(ii), hMCEffHLTB->GetBinContent(ii) );
    printf(" - OFF Purity vs MC  : %.3f | OFF Eff vs MC  : %.3f \n", hMCPurityOFFB->GetBinContent(ii), hMCEffOFFB->GetBinContent(ii) );
  }


  return;
}

//#######################################################################
void CanvasTriggerFactors() {

  TCanvas *cTriggerR = new TCanvas("cTriggerR", "Reduction Factors", 10, 10, 1400, 800);
  cTriggerR->Divide(1,2);

  TCanvas *cTriggerEP = new TCanvas("cTriggerEP", "Purity - Efficiency", 10, 10, 1400, 800);
  cTriggerEP->Divide(3,2);
 
  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  Float_t ptMin = -0.25;
  Float_t ptMax = 12.75;
  Int_t nBins   = 26;

  Float_t maxLine = 0.9;
  Float_t minLine = 0.8;

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TLine* line1 = new TLine(0.,maxLine,ptMax,maxLine);
  line1->SetLineStyle(2);
  line1->SetLineColor(kGreen);

  TLine* line2 = new TLine(0.,minLine,ptMax,minLine);
  line2->SetLineStyle(2);
  line2->SetLineColor(kRed);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TH1F* hOFFRedux   = new TH1F("hOFFRedux", "Reduction Factors", nBins, ptMin, ptMax); 
  hOFFRedux->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  hOFFRedux->GetYaxis()->SetTitle("Reduction #left(#frac{N Events_{Total}}{N Events_{Triggered}}#right)");

  TH1F* hOFFReduxS  = new TH1F("hOFFReduxS", "", nBins, ptMin, ptMax); 
  TH1F* hHLTRedux   = new TH1F("hHLTRedux",  "", nBins, ptMin, ptMax); 
  TH1F* hHLTReduxS  = new TH1F("hHLTReduxS", "", nBins, ptMin, ptMax); 
  TH1F* hMCRedux    = new TH1F("hMCRedux",   "", nBins, ptMin, ptMax); 
  TH1F* hMCReduxS   = new TH1F("hMCReduxS",  "", nBins, ptMin, ptMax); 

  TH1F* hOFFReduxW   = new TH1F("hOFFReduxW", "Reduction Factors - weighted nTracks", nBins, ptMin, ptMax); 
  hOFFReduxW->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  hOFFReduxW->GetYaxis()->SetTitle("Reduction #left(#frac{N Events_{Total}}{N Events_{Triggered}}#right)");

  TH1F* hOFFReduxWS  = new TH1F("hOFFReduxWS", "", nBins, ptMin, ptMax); 
  TH1F* hHLTReduxW   = new TH1F("hHLTReduxW",  "", nBins, ptMin, ptMax); 
  TH1F* hHLTReduxWS  = new TH1F("hHLTReduxWS", "", nBins, ptMin, ptMax); 
  TH1F* hMCReduxW    = new TH1F("hMCReduxW",   "", nBins, ptMin, ptMax); 
  TH1F* hMCReduxWS   = new TH1F("hMCReduxWS",  "", nBins, ptMin, ptMax); 

  FillReduxHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggered")),    hOFFRedux,  hOFFReduxW);
  FillReduxHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggeredSel")), hOFFReduxS, hOFFReduxWS);
  FillReduxHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggered")),    hHLTRedux,  hHLTReduxW);
  FillReduxHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggeredSel")), hHLTReduxS, hHLTReduxWS);
  FillReduxHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggered")),     hMCRedux,   hMCReduxW);
  FillReduxHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggeredSel")),  hMCReduxS,  hMCReduxWS);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TLegend* legR = new TLegend(0.65, 0.12, 0.89, 0.37, "Event Selection");
  legR->AddEntry(hOFFRedux, "OFF Trig - All Events","LP");
  legR->AddEntry(hOFFReduxS,"OFF Trig - AliPhysSel - OFF PrimVertex","LP");
  legR->AddEntry(hHLTRedux, "HLT Trig - All Events","LP");
  legR->AddEntry(hHLTReduxS,"HLT Trig - AliPhysSel - OFF PrimVertex","LP");
  legR->AddEntry(hMCRedux,  "MC Trig - All Events","LP");
  legR->AddEntry(hMCReduxS, "MC Trig - AliPhysSel - OFF PrimVertex","LP");
  legR->SetFillColor(kWhite);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  DrawHistogramTrigger(cTriggerR, 1, hOFFRedux,  0, kTRUE);
  DrawHistogramTrigger(cTriggerR, 1, hOFFReduxS, 1, kTRUE);
  DrawHistogramTrigger(cTriggerR, 1, hHLTRedux,  2, kTRUE);
  DrawHistogramTrigger(cTriggerR, 1, hHLTReduxS, 3, kTRUE);
  DrawHistogramTrigger(cTriggerR, 1, hMCRedux,   4, kTRUE);
  DrawHistogramTrigger(cTriggerR, 1, hMCReduxS,  5, kTRUE);
  legR->Draw();

  DrawHistogramTrigger(cTriggerR, 2, hOFFReduxW,  0, kTRUE);
  DrawHistogramTrigger(cTriggerR, 2, hOFFReduxWS, 1, kTRUE);
  DrawHistogramTrigger(cTriggerR, 2, hHLTReduxW,  2, kTRUE);
  DrawHistogramTrigger(cTriggerR, 2, hHLTReduxWS, 3, kTRUE);
  DrawHistogramTrigger(cTriggerR, 2, hMCReduxW,   4, kTRUE);
  DrawHistogramTrigger(cTriggerR, 2, hMCReduxWS,  5, kTRUE);
  legR->Draw();

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  // ---------------------------------------------------------------
  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TH1F* hOFFPurityHLT  = new TH1F("hOFFPurityHLT", "HLT Trigger Purity vs OFF", nBins, ptMin, ptMax);
  hOFFPurityHLT->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  hOFFPurityHLT->GetYaxis()->SetTitle("Purity #left(#frac{HLT_{Triggered} - HLT_{Fake}}{HLT_{Triggered}}#right)");
  hOFFPurityHLT->SetMinimum(0.0);
  hOFFPurityHLT->SetMaximum(1.2);
  TH1F* hOFFPurityHLTS = new TH1F("hOFFPurityHLTS", "", nBins, ptMin, ptMax);

  TH1F* hMCPurityHLT  = new TH1F("hMCPurityHLT", "HLT Trigger Purity vs MC", nBins, ptMin, ptMax);
  hMCPurityHLT->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  hMCPurityHLT->GetYaxis()->SetTitle("Purity #left(#frac{HLT_{Triggered} - HLT_{Fake}}{HLT_{Triggered}}#right)");
  hMCPurityHLT->SetMinimum(0.0);
  hMCPurityHLT->SetMaximum(1.2);
  TH1F* hMCPurityHLTS = new TH1F("hMCPurityHLTS", "", nBins, ptMin, ptMax);

  TH1F* hMCPurityOFF  = new TH1F("hMCPurityOFF", "OFF Trigger Purity vs MC", nBins, ptMin, ptMax);
  hMCPurityOFF->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  hMCPurityOFF->GetYaxis()->SetTitle("Purity #left(#frac{OFF_{Triggered} - OFF_{Fake}}{OFF_{Triggered}}#right)");
  hMCPurityOFF->SetMinimum(0.0);
  hMCPurityOFF->SetMaximum(1.2);
  TH1F* hMCPurityOFFS = new TH1F("hMCPurityOFFS", "", nBins, ptMin, ptMax);


  FillRatioHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggered")),    
			static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToOFF")), hOFFPurityHLT);

  FillRatioHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggeredSel")),
			static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToOFFSel")), hOFFPurityHLTS);

  FillRatioHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggered")),    
			static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToMC")), hMCPurityHLT);

  FillRatioHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNHLTTriggeredSel")),
			static_cast<TH1F*>(fObjArray->FindObject("fNHLTFakeToMCSel")), hMCPurityHLTS);

  FillRatioHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggered")),    
			static_cast<TH1F*>(fObjArray->FindObject("fNOFFFakeToMC")), hMCPurityOFF);

  FillRatioHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggeredSel")),
			static_cast<TH1F*>(fObjArray->FindObject("fNOFFFakeToMCSel")), hMCPurityOFFS);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TLegend* legPH = new TLegend(0.65, 0.72, 0.89, 0.87, "Event Selection");
  legPH->AddEntry(hOFFPurityHLT, "HLT Trig - All Events","LP");
  legPH->AddEntry(hOFFPurityHLTS,"HLT Trig - AliPhysSel - OFF PrimVertex","LP");
  legPH->SetFillColor(kWhite);

  TLegend* legPO = new TLegend(0.65, 0.72, 0.89, 0.87, "Event Selection");
  legPO->AddEntry(hMCPurityOFF, "OFF Trig - All Events","LP");
  legPO->AddEntry(hMCPurityOFFS,"OFF Trig - AliPhysSel - OFF PrimVertex","LP");
  legPO->SetFillColor(kWhite);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  DrawHistogramTrigger(cTriggerEP, 1, hOFFPurityHLT,  0, kFALSE);
  DrawHistogramTrigger(cTriggerEP, 1, hOFFPurityHLTS, 1, kFALSE);
  line1->Draw();
  line2->Draw();
  legPH->Draw(); 

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  DrawHistogramTrigger(cTriggerEP, 2, hMCPurityHLT,  0, kFALSE);
  DrawHistogramTrigger(cTriggerEP, 2, hMCPurityHLTS, 1, kFALSE);
  line1->Draw();
  line2->Draw();
  legPH->Draw(); 

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  DrawHistogramTrigger(cTriggerEP, 3, hMCPurityOFF,  0, kFALSE);
  DrawHistogramTrigger(cTriggerEP, 3, hMCPurityOFFS, 1, kFALSE);
  line1->Draw();
  line2->Draw();
  legPO->Draw(); 

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  // ---------------------------------------------------------------
  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TH1F* hOFFEffHLT  = new TH1F("hOFFEffHLT", "HLT Trigger Efficiency vs OFF", nBins, ptMin, ptMax);
  hOFFEffHLT->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  hOFFEffHLT->GetYaxis()->SetTitle("Efficiency #left(#frac{OFF_{Triggered} - HLT_{Miss}}{OFF_{Triggered}}#right)");
  hOFFEffHLT->SetMinimum(0.0);
  hOFFEffHLT->SetMaximum(1.2);
  TH1F* hOFFEffHLTS = new TH1F("hOFFEffHLTS", "", nBins, ptMin, ptMax);

  TH1F* hMCEffHLT  = new TH1F("hMCEffHLT","HLT Trigger Efficiency vs MC", nBins, ptMin, ptMax);
  hMCEffHLT->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  hMCEffHLT->GetYaxis()->SetTitle("Efficiency #left(#frac{MC_{Triggered} - HLT_{Miss}}{MC_{Triggered}}#right)");
  hMCEffHLT->SetMinimum(0.0);
  hMCEffHLT->SetMaximum(1.2);
  TH1F* hMCEffHLTS = new TH1F("hMCEffHLTS", "", nBins, ptMin, ptMax);

  TH1F* hMCEffOFF  = new TH1F("hMCEffOFF","OFF Trigger Efficiency vs MC", nBins, ptMin, ptMax);
  hMCEffOFF->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  hMCEffOFF->GetYaxis()->SetTitle("Efficiency #left(#frac{MC_{Triggered} - OFF_{Miss}}{MC_{Triggered}}#right)");
  hMCEffOFF->SetMinimum(0.0);
  hMCEffOFF->SetMaximum(1.2);
  TH1F* hMCEffOFFS = new TH1F("hMCEffOFFS", "", nBins, ptMin, ptMax);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  FillRatioHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggered")),    
		      static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToOFF")), hOFFEffHLT);

  FillRatioHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNOFFTriggeredSel")),
		      static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToOFFSel")), hOFFEffHLTS);

  FillRatioHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggered")),    
		      static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToMC")), hMCEffHLT);

  FillRatioHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggeredSel")),
		      static_cast<TH1F*>(fObjArray->FindObject("fNHLTMissToMCSel")), hMCEffHLTS);

  FillRatioHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggered")),    
		      static_cast<TH1F*>(fObjArray->FindObject("fNOFFMissToMC")), hMCEffOFF);

  FillRatioHistogramsPt(static_cast<TH1F*>(fObjArray->FindObject("fNMCTriggeredSel")),
		      static_cast<TH1F*>(fObjArray->FindObject("fNOFFMissToMCSel")), hMCEffOFFS);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  TLegend* legEH = new TLegend(0.65, 0.70, 0.89, 0.85, "Event Selection");
  legEH->AddEntry(hOFFEffHLT, "HLT Trig - All Events","LP");
  legEH->AddEntry(hOFFEffHLTS,"HLT Trig - AliPhysSel - OFF PrimVertex","LP");
  legEH->SetFillColor(kWhite);

  TLegend* legEO = new TLegend(0.65, 0.70, 0.89, 0.85, "Event Selection");
  legEO->AddEntry(hMCEffOFF, "OFF Trig - All Events","LP");
  legEO->AddEntry(hMCEffOFFS,"OFF Trig - AliPhysSel - OFF PrimVertex","LP");
  legEO->SetFillColor(kWhite);

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  DrawHistogramTrigger(cTriggerEP, 4, hOFFEffHLT, 0, kFALSE);
  DrawHistogramTrigger(cTriggerEP, 4, hOFFEffHLTS, 1, kFALSE);
  line1->Draw();
  line2->Draw();
  legEH->Draw(); 

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  DrawHistogramTrigger(cTriggerEP, 5, hMCEffHLT,  0, kFALSE);
  DrawHistogramTrigger(cTriggerEP, 5, hMCEffHLTS, 1, kFALSE);
  line1->Draw();
  line2->Draw();
  legEH->Draw(); 

  // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

  DrawHistogramTrigger(cTriggerEP, 6, hMCEffOFF,  0, kFALSE);
  DrawHistogramTrigger(cTriggerEP, 6, hMCEffOFFS, 1, kFALSE);
  line1->Draw();
  line2->Draw();
  legEO->Draw(); 

  cTriggerR->SaveAs(Form("TriggerHists_%d.png",++fgImgIdx));
  cTriggerEP->SaveAs(Form("TriggerHists_%d.png",++fgImgIdx));

  return;
}

//#######################################################################
void CreateCanvasCuts() {

  Char_t * type[] = {"MC","OFF","HLT","HLT/MC","HLT/OFF","OFF/MC",
		     "MC","OFF","HLT","HLT/MC","HLT/OFF","OFF/MC" };
  Char_t * histType[] = {"Pt","Mult"};

  Float_t ptXMin = 0.3;
  Float_t ptXMax = 20.;
  Float_t ptYMin = 1.;
  Float_t ptYMax = 1e6;

  Float_t multXMin = 3.;
  Float_t multXMax = 50.;
  Float_t multYMin = 1.;
  Float_t multYMax = 1e5;

  Int_t   maxPads = 12;
  Int_t   iMax = 2;

  // ----------------------------------------------------------------------------

  TCanvas *cEvCuts = new TCanvas("cEvCuts", "Pt Event Selection - Cut Studies", 10, 10, 1400, 800);
  cEvCuts->Divide(3,2);

  TCanvas *cTrCuts = new TCanvas("cTrCuts", "Pt Track Selection - Cut Studies", 10, 10, 1400, 800);
  cTrCuts->Divide(3,2);

  TCanvas *cEvCuts2 = new TCanvas("cEvCuts2", "Mult Event Selection - Cut Studies", 10, 10, 1400, 800);
  cEvCuts2->Divide(3,2);

  TCanvas *cTrCuts2 = new TCanvas("cTrCuts2", "Mult Track Selection - Cut Studies", 10, 10, 1400, 800);
  cTrCuts2->Divide(3,2);

  // ----------------------------------------------------------------------------

  TLegend* legEvCuts[maxPads];
  TLegend* legTrCuts[maxPads];

  for ( Int_t idx=0; idx < maxPads; idx++ ) {

    if ( idx < maxPads/2) { // -- Pt
    
      if ( idx >= 3 && idx < 6) { // ratio
	legEvCuts[idx] = new TLegend(0.14, 0.75, 0.55, 0.88, 
				     Form("Event Selection - P_{t} Cut Studies (%s)", type[idx])); 

	legTrCuts[idx] = new TLegend(0.14, 0.75, 0.55, 0.88, 
				     Form("Track Selection - P_{t} Cut Studies (%s)", type[idx])); 
      }
      else {
	legEvCuts[idx] = new TLegend(0.14, 0.16, 0.55, 0.29, 
				     Form("Event Selection - P_{t} Cut Studies (%s)", type[idx])); 
  
	legTrCuts[idx] = new TLegend(0.14, 0.16, 0.55, 0.29, 
				     Form("Track Selection - P_{t} Cut Studies (%s)", type[idx])); 
      }
      
      legEvCuts[idx]->SetFillColor(kWhite);
      legTrCuts[idx]->SetFillColor(kWhite);
    }
    else { // -- Multiplicity
      legEvCuts[idx] = new TLegend(0.48, 0.75, 0.89, 0.88, 
				   Form("Event Selection - Multiplicity Cut Studies (%s)", type[idx])); 
      legEvCuts[idx]->SetFillColor(kWhite);
      
      legTrCuts[idx] = new TLegend(0.48, 0.75, 0.89, 0.88, 
				   Form("Track Selection - Multiplicity Cut Studies (%s)", type[idx])); 
      legTrCuts[idx]->SetFillColor(kWhite);
    }
  }

  // ----------------------------------------------------------------------------

  for ( Int_t hType = 0; hType < 2; hType++ ) {
    
    TCanvas* can = NULL;
    TCanvas* can2 = NULL;
    if ( hType == 0 ) {
      can  = cEvCuts;
      can2 = cEvCuts2;
    }
    else if ( hType == 1 ) {
      can  = cTrCuts;
      can2  = cTrCuts2;
    }

    for ( Int_t ipad=1; ipad <= maxPads; ipad++ ) {

      Int_t pad = ipad;

      // -- Legends
      TLegend* leg = NULL;
      if ( hType == 0 )
	leg  = legEvCuts[pad-1];
      else if ( hType == 1 )
	leg  = legTrCuts[pad-1];

      // -- Pt / Multiplicity
      Bool_t xLog = kTRUE;
      Int_t hIdx = 0;
      Float_t xMin, xMax, yMin, yMax;
      if ( ipad <= maxPads/2 ) { // -- Pt
	hIdx = 0; 
	xLog = kTRUE;
	xMin = ptXMin;
	xMax = ptXMax;
	yMin = ptYMin;
	yMax = ptYMax;
      }
      else { // -- Multiplicity
	hIdx = 1; 
	xLog = kFALSE;
	pad = ipad - maxPads/2;
	can = can2;
	xMin = multXMin;
	xMax = multXMax;
	yMin = multYMin;
	yMax = multYMax;
      }
     
      // ----------------------------------------------------------------------------

      for ( Int_t idx = 0; idx <= iMax; idx++ ) {

	Int_t iIdx = idx;
	if ( hType == 1 ) // -- Track selection
	  iIdx = idx+2;

	TH1F* hO  = static_cast<TH1F*>(fObjArray->FindObject(Form("fHistOFF%s_%d", histType[hIdx],iIdx)));
	TH1F* hH = static_cast<TH1F*>(fObjArray->FindObject(Form("fHistHLT%s_%d", histType[hIdx], iIdx)));
	TH1F* hM = static_cast<TH1F*>(fObjArray->FindObject(Form("fHistMC%s_%d", histType[hIdx], iIdx)));
	
	if ( hType == 0 ) {  // Event selection
	  hO->SetTitle(Form("OFF %s Distribution - Event Selection", histType[hIdx]));
	  hH->SetTitle(Form("HLT %s Distribution - Event Selection", histType[hIdx]));
	  hM->SetTitle(Form("MC %s Distribution - Event Selection", histType[hIdx]));
	}
	else {
	  hO->SetTitle(Form("OFF %s Distribution - Track Selection", histType[hIdx]));
	  hH->SetTitle(Form("HLT %s Distribution - Track Selection", histType[hIdx]));
	  hM->SetTitle(Form("MC %s Distribution - Track Selection", histType[hIdx]));
	}
	
	hO->SetAxisRange(xMin, xMax);
	hH->SetAxisRange(xMin, xMax);
	hM->SetAxisRange(xMin, xMax);
	hO->SetAxisRange(yMin, yMax, "Y");
	hH->SetAxisRange(yMin, yMax, "Y");
	hM->SetAxisRange(yMin, yMax, "Y");

       	if (pad == 1 || pad == 7) {
	  DrawHistogram(can, pad, hM, idx, kFALSE, kTRUE, xLog);
	  leg->AddEntry(hM, fgkSelectionCutsMC[iIdx],"L");
	  leg->Draw();
	}
	else if (pad == 2 || pad == 8 ) {
	  DrawHistogram(can, pad, hO, idx, kFALSE, kTRUE, xLog);
	  leg->AddEntry(hO, fgkSelectionCuts[iIdx],"L");
	  leg->Draw();
	}
	else if (pad == 3 || pad == 9 ) {
	  DrawHistogram(can, pad, hH, idx, kFALSE, kTRUE, xLog);
	  leg->AddEntry(hH, fgkSelectionCuts[iIdx],"L");
	  leg->Draw();
	}
	else if (pad == 4 || pad == 10 ) {

	  TH1F* hClone0= static_cast<TH1F*>(hH->Clone());
	  DivideHist(hClone0, hH, hM);

	  if (hIdx == 0)
	    SetHist(hClone0, hType, type[pad-1], histType[hIdx], 0., 3.);
	  else
	    SetHist(hClone0, hType, type[pad-1], histType[hIdx], 0., 10.);

	  DrawHistogram(can, pad, hClone0, idx, kFALSE, kFALSE, xLog);
	  leg->AddEntry(hClone0, fgkSelectionCutsRatioMC[iIdx],"L");
	  leg->Draw();
	}
	else if (pad == 5 || pad == 11 ) {

	  TH1F* hClone1= static_cast<TH1F*>(hH->Clone());
	  DivideHist(hClone1, hH, hO);
	  
	  if (hIdx == 0)
	    SetHist(hClone1, hType, type[pad-1], histType[hIdx], 0., 3.);
	  else
	    SetHist(hClone1, hType, type[pad-1], histType[hIdx], 0., 10.);

	  DrawHistogram(can, pad, hClone1, idx, kFALSE, kFALSE, xLog);
	  leg->AddEntry(hClone1, fgkSelectionCuts[iIdx],"L");
	  leg->Draw();
	}
	else if (pad == 6 || pad == 12 ) {

	  TH1F* hClone2= static_cast<TH1F*>(hO->Clone());
	  DivideHist(hClone2, hO, hM);

	  if (hIdx == 0)
	    SetHist(hClone2, hType, type[pad-1], histType[hIdx], 0., 3.);
	  else
	    SetHist(hClone2, hType, type[pad-1], histType[hIdx], 0., 10.);
	      
	  DrawHistogram(can, pad, hClone2, idx, kFALSE, kFALSE, xLog);
	  leg->AddEntry(hClone2, fgkSelectionCutsRatioMC[iIdx],"L");
	  leg->Draw();
	}
      } // for ( Int_t idx = 0; idx <=2; idx++ ) {
    } // for ( Int_t pad=1; pad <=4; pad++ ) {
  } // for ( Int_t hType = 0; hType < 4; hType++ ) {
  

  if (cEvCuts)  cEvCuts->SaveAs(Form("TriggerHists_Cuts_%d_B.png",++fgImgIdx));
  if (cEvCuts2) cEvCuts2->SaveAs(Form("TriggerHists_Cuts_%d_B.png",++fgImgIdx));
  if (cTrCuts)  cTrCuts->SaveAs(Form("TriggerHists_Cuts_%d_B.png",++fgImgIdx));
  if (cTrCuts2) cTrCuts2->SaveAs(Form("TriggerHists_Cuts_%d_B.png",++fgImgIdx));

  return;
}

//#######################################################################
void CreateCanvasDistributions( Int_t maxLayer, Int_t draw, Int_t data, Int_t minLayer ) {

  Float_t maxLine = 0.9;
  Float_t minLine = 0.8;

  Int_t compHist = 3;
  Float_t ptMin = 0.3;
  Float_t ptMax = 40.;

  Char_t * type[] = {"MC","OFF","HLT"};
  Char_t * histType[] = {"Pt","Mult"};

  // ----------------------------------------------------------------------------

  TCanvas *cSpectra = NULL;
  TCanvas *cSpectraEff = NULL;
  TCanvas *cMult = NULL;
  TCanvas *cMultEff = NULL;

  TLegend* legPt[6];
  TLegend* legPtEff[6];
  TLegend* legMult[6];
  TLegend* legMultEff[6];

  // ----------------------------------------------------------------------------
  Int_t canvasRows = 1;

  if ( data == 2 )
    canvasRows = 2;

  if ( draw == 0 || draw == 2 ) {
    cSpectra = new TCanvas("cSpectra", "Triggered - P_{t} Distributions", 10, 10, 1400, 800);
    cSpectra->Divide(3,canvasRows);

   cSpectraEff = new TCanvas("cSpectraEff", "Triggered - P_{t} Efficiencies", 10, 10, 1400, 800);
   cSpectraEff->Divide(3,canvasRows);
  }

  if ( draw == 1 || draw == 2 ) {
    cMult = new TCanvas("cMult", "Triggered - Multiplicity Distributions", 10, 10, 1400, 800);
    cMult->Divide(3,canvasRows);

    cMultEff = new TCanvas("cMultEff", "Triggered - Multiplicity Efficiencies", 10, 10, 1400, 800);
    cMultEff->Divide(3,canvasRows);
  }

  // ----------------------------------------------------------------------------
  for ( Int_t hType = 0; hType < 2; hType++ ) {

    if ( hType != draw && draw != 2)
      continue;

    for ( Int_t pad=1; pad <=6; pad++ ) {

      if ( data == 0 && pad > 3 ) 
	continue;

      if ( data == 1 && pad < 4 ) 
	continue;
      
      TCanvas* can = NULL;
      TCanvas* canEff = NULL;

      TLegend* leg = NULL;
      TLegend* legEff = NULL;

      Bool_t xLog = kTRUE;
      Char_t* trigType = "";
      Char_t* dataType = "OFF";

      Int_t iIdx = 0;
      Int_t layer = 0;

      if ( pad == 1 || pad == 4 ) {
	iIdx = 0;
	trigType = "MC";
      }
      else if ( pad == 2 || pad == 5 ) {
	iIdx = 1;
	trigType = "OFF";
      }
      else if ( pad == 3 || pad == 6 ) {
	iIdx = 2;
	trigType = "HLT";
      }
      else
	continue;

      if ( pad > 3 )
	dataType = "MC";
      
      if ( hType == 0 ) {
	can  = cSpectra;
	canEff  = cSpectraEff;
	leg  = legPt[pad-1];
	legEff  = legPtEff[pad-1];
	xLog = kTRUE;
      }
      else if ( hType == 1 ) {
	can  = cMult;
	canEff  = cMultEff;
	leg  = legMult[pad-1];
	legEff  = legMultEff[pad-1];
	xLog = kFALSE;
      }

      TH1F * hcp = static_cast<TH1F*>(fObjArray->FindObject(Form("fHist%s%s_%d", dataType,
								 histType[hType],compHist)));
      
      if ( hType == 0 )
      	hcp->SetAxisRange(ptMin, ptMax);

      hcp->SetMinimum(10);
       
      DrawHistogram(can, pad, hcp, layer, kFALSE, kFALSE, xLog);
      
      leg = new TLegend(0.58, 0.70, 0.88, 0.88, Form("Trigger Selection: %s Tracks",type[iIdx])); 
      leg->AddEntry(hcp, "Track Selection","L");
      leg->SetFillColor(kWhite);

      legEff = new TLegend(0.15, 0.70, 0.45, 0.88, Form("Trigger Selection: %s Tracks",type[iIdx])); 
      legEff->SetFillColor(kWhite);

      // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
      ++layer;

      for ( Int_t idx= 0; idx < fgkNTrigger; idx++ ) { 

	if ( idx >= maxLayer )
	  continue;

	if ( idx < minLayer )
	  continue;
  
	Int_t iter = idx;

	TH1F *hist = static_cast<TH1F*>(fObjArray->FindObject(Form("fHist%s_%sTriggered_%s_%d", 
								   dataType, histType[hType], trigType, iter)));
    
	if ( hType == 0 )
	  hist->SetAxisRange(ptMin, ptMax);
	DrawHistogram(can, pad, hist, layer, kFALSE, kTRUE, xLog);

	leg->AddEntry(hist, fgkTrigger[iter],"L");
	leg->Draw();

	// ------------------------------------------------

	TH1F *eff = static_cast<TH1F*>(hist->Clone()); 
	DivideHist(eff,hist,hcp);
	if ( hType == 0 )
	  eff->SetAxisRange(ptMin, ptMax);

	eff->SetAxisRange(0., 1.4, "Y");
	eff->SetTitle(Form("%s Trigger Efficiency - applied on %s data", type[iIdx], dataType));
	eff->GetYaxis()->SetTitle(Form("%s Trigger Efficiecy", type[iIdx]));

	DrawHistogram(canEff, pad, eff, layer-1, kFALSE, kFALSE, xLog);
	legEff->AddEntry(eff, fgkTrigger[iter],"L");
	legEff->Draw();

	// ------------------------------------------------

	// -- Draw line at maxLine
	if ( idx == minLayer ){

	  Float_t max = 50.;
	  if ( hType == 0)
	    max = ptMax;

	  TLine* line = NULL;
	  line = new TLine(0.,maxLine,max,maxLine);
	  line->SetLineStyle(2);
	  line->SetLineColor(kGreen);
	  line->Draw();

	  TLine* line2 = NULL;
	  line2 = new TLine(0.,minLine,max,minLine);
	  line2->SetLineStyle(2);
	  line2->SetLineColor(kRed);
	  line2->Draw();
	}

	// ------------------------------------------------

	++layer;
      } // for ( Int_t idx= 0; idx < 3; idx++ ) {    

    } // for ( Int_t pad=1; pad <=6; pad++ ) {
  
  } // for ( Int_t hType = 0; hType < 2; hType++ ) {


  if (cSpectra) cSpectra->SaveAs(Form("TriggerHists_Pt_%d.png",++fgImgIdx));
  if (cSpectraEff) cSpectraEff->SaveAs(Form("TriggerHists_PtEff_%d.png",++fgImgIdx));

  if (cMult) cMult->SaveAs(Form("TriggerHists_Mult_%d_B.png",++fgImgIdx));
  if (cMultEff) cMultEff->SaveAs(Form("TriggerHists_MultEff_%d_B.png",++fgImgIdx));

  return;
}

// ----------------------------------------------------------------------
// --                      Helper Functions                            --
// ----------------------------------------------------------------------

//#######################################################################
void SetHist(TH1F* hist, Int_t hType, Char_t* dataType, Char_t* histType, Float_t minY, Float_t maxY ) {

  if (hType == 0)
    hist->SetTitle(Form("Ratio %s - %s Distribution - Event Selection", dataType, histType));
  else
    hist->SetTitle(Form("Ratio %s - %s Distribution - Track Selection", dataType, histType));

  hist->SetAxisRange(minY, maxY, "Y");
  hist->GetYaxis()->SetTitle(Form("Ratio %s", dataType));
}

//#######################################################################
void DivideHist(TH1F* h0, TH1F* h1, TH1F* h2) {
  // Divides h1 / h2, does proper error handling and saves in h0

  for(Int_t iBin=1; iBin < h0->GetXaxis()->GetNbins(); iBin++) {
    Float_t hist1Content = h1->GetBinContent(iBin);
    Float_t hist1Error = h1->GetBinError(iBin);
    Float_t hist2Content = h2->GetBinContent(iBin);
    Float_t hist2Error = h2->GetBinError(iBin);
    
    if ( hist1Content == 0. && hist1Content == 0. ){
      h0->SetBinContent(iBin, 1.);
      h0->SetBinError(iBin, 0.);
    }

    if(hist1Content<=0.0) continue;
    if(hist2Content<=0.0) continue;

    Float_t ratio = hist1Content/hist2Content;
    Float_t relError1 = hist1Error/hist1Content;
    Float_t relError2 = hist2Error/hist2Content;
    Float_t error = ratio*TMath::Sqrt(relError1*relError1+relError2*relError2);

    h0->SetBinContent(iBin, ratio);
    h0->SetBinError(iBin, error);
  }

  return;
}

//#######################################################################
void FillIntegrated(TH1F* hist, TH1F* histInt ) {
  // see header file for class documentation

  Double_t oldBin = 0.;
  for ( Int_t idx = 0; idx < hist->GetNbinsX() ; ++idx ){
    histInt->SetBinContent(idx, hist->GetBinContent(idx) + oldBin);
    oldBin = hist->GetBinContent(idx) + oldBin;
  }
  histInt->SetEntries(hist->GetNbinsX());

  return;
}

//#######################################################################
void FillReduxHistograms( TH1F* hN, TH1F* hRedux ) {

  Float_t base       = hN->GetBinContent(1);
  Float_t baseRelErr = hN->GetBinError(1) / base;

  for ( Int_t idx = 2; idx <= 2*(fgkNTrigger+1) ; ++idx ){

    Int_t bin = idx-1;    
    if ( idx == fgkNTrigger+2 ) {
      base = hN->GetBinContent(fgkNTrigger+2);
      baseRelErr = hN->GetBinError(fgkNTrigger+2) / base;
      continue;
    }
    else if( idx > fgkNTrigger+2 )
      bin = idx-2;

    Float_t cont       = hN->GetBinContent(idx);

    if (cont <=0.) continue;
    Float_t contRelErr = hN->GetBinError(idx) / cont;
    
    Float_t ratio      = base / cont;
    Float_t error      = ratio * TMath::Sqrt(contRelErr*contRelErr + baseRelErr*baseRelErr);
    
    hRedux->SetBinContent(bin, ratio);
    hRedux->SetBinError(bin, error);
  }

  return;
}

//#######################################################################
void FillReduxHistogramsPt( TH1F* hN, TH1F* hRedux, TH1F* hReduxW ) {

  Float_t base       = hN->GetBinContent(1);
  Float_t baseRelErr = hN->GetBinError(1) / base;

  TH1F* hist = hRedux;
  Int_t trgIdx = 0;
 
  for ( Int_t idx = 2; idx <= 2*(fgkNTrigger+1) ; ++idx ){

    if ( idx == fgkNTrigger+2 ) {
      base = hN->GetBinContent(fgkNTrigger+2);
      baseRelErr = hN->GetBinError(fgkNTrigger+2) / base;
      hist = hReduxW;
      trgIdx = 0;
      continue;
    }

    if ( fgkTriggerPt[trgIdx] > 100. ) {
      ++trgIdx;
      continue;
    }

    Int_t   nbin        = hist->FindBin(fgkTriggerPt[trgIdx]);
    ++trgIdx;

    Float_t cont       = hN->GetBinContent(idx);

    if (cont <=0.) continue;
    Float_t contRelErr = hN->GetBinError(idx) / cont;
    
    Float_t ratio      = base / cont;
    Float_t error      = ratio * TMath::Sqrt(contRelErr*contRelErr + baseRelErr*baseRelErr);
    
    hist->SetBinContent(nbin, ratio);
    hist->SetBinError(nbin, error);
  }
    
  return;
}

//#######################################################################
void FillRatioHistograms( TH1F* hN, TH1F* hF, TH1F* hRatio ) {

  for ( Int_t idx = 2; idx <= 2*(fgkNTrigger+1) ; ++idx ){

    Int_t bin = idx-1;    
    if ( idx == fgkNTrigger+2 )
      continue;
    else if( idx > fgkNTrigger+2 )
      bin = idx-2;

    Float_t h1Cont = hF->GetBinContent(idx);
    Float_t h2Cont = hN->GetBinContent(idx);

    if ( h2Cont <= 0. || h1Cont <= 0. ) continue;

    Float_t h1RelErr = hF->GetBinError(idx) / h1Cont;
    Float_t h2RelErr = hN->GetBinError(idx) / h2Cont;

    Float_t ratio    = (h2Cont - h1Cont) / h2Cont ;
    Float_t error    = h1Cont / h2Cont * TMath::Sqrt(h1RelErr*h2RelErr + h1RelErr*h2RelErr);

    hRatio->SetBinContent(bin, ratio);
    hRatio->SetBinError(bin, error);
  }

  return;
}

//#######################################################################
void FillRatioHistogramsPt( TH1F* hN, TH1F* hF, TH1F* hRatio ) {

  Int_t trgIdx = 0;

  for ( Int_t idx = 2; idx < fgkNTrigger+2 ; ++idx ){

    Int_t   nbin = hRatio->FindBin(fgkTriggerPt[trgIdx]);
    ++trgIdx;

    Float_t h1Cont = hF->GetBinContent(idx);
    Float_t h2Cont = hN->GetBinContent(idx);

    if ( h2Cont <= 0. || h1Cont <= 0. ) continue;

    Float_t h1RelErr = hF->GetBinError(idx) / h1Cont;
    Float_t h2RelErr = hN->GetBinError(idx) / h2Cont;

    Float_t ratio    = (h2Cont - h1Cont) / h2Cont ;
    Float_t error    = h1Cont / h2Cont * TMath::Sqrt(h1RelErr*h2RelErr + h1RelErr*h2RelErr);

    hRatio->SetBinContent(nbin, ratio);
    hRatio->SetBinError(nbin, error);
  }

  return;
}

// ----------------------------------------------------------------------
// --                       Draw Functions                             --
// ----------------------------------------------------------------------

//#######################################################################
void DrawHistogramTrigger( TCanvas* canvas, Int_t idx, TH1* hist, Int_t iLayer, Bool_t bLogY) {
  // see header file for class documentation

  if ( hist == NULL )
    return;

  TVirtualPad* pad = canvas->cd(idx);
  
  if ( bLogY && hist->GetEntries() != 0 )
    pad->SetLogy();

  pad->SetGridy();
  pad->SetGridx();

  switch(iLayer) {
  case 0 : 
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(kBlack);
    break;;
  case 1 : 
    hist->SetLineColor(kBlack-6); 
    hist->SetMarkerStyle(24);
    hist->SetMarkerColor(kBlack-6);
    break;;
  case 2 : 
    hist->SetLineColor(kBlue); 
    hist->SetMarkerStyle(21);
    hist->SetMarkerColor(kBlue);
    break;;
  case 3 : 
    hist->SetLineColor(kBlue-6); 
    hist->SetMarkerStyle(25);
    hist->SetMarkerColor(kBlue-6);
    break;;
  case 4 : 
    hist->SetLineColor(kGreen); 
    hist->SetMarkerStyle(22);
    hist->SetMarkerColor(kGreen);
    break;;
  case 5 : 
    hist->SetLineColor(kGreen-6); 
    hist->SetMarkerStyle(26);
    hist->SetMarkerColor(kGreen-6);
    break;;
  }

  hist->SetStats(kFALSE);

  if ( !iLayer )
    hist->Draw("E");
  else 
    hist->Draw("E,SAME");
    
  return;
}

//#######################################################################
void DrawHistogram( TCanvas* canvas, Int_t idx, TH1* hist, Int_t iLayer, 
		    Bool_t bScale, Bool_t bLogY, Bool_t bLogX) {
  // Draw histogram and set pad properties

  if ( hist == NULL )
    return;

  TVirtualPad* pad = canvas->cd(idx);
  
  if ( bScale ) 
    hist->Scale( 1./hist->GetBinWidth(0) );
  
  if ( bLogY && hist->GetEntries() != 0 )
    pad->SetLogy();

  if ( bLogX && hist->GetEntries() != 0 )
    pad->SetLogx();
  
  pad->SetGridy();
  pad->SetGridx();

  if ( !strcmp(canvas->GetName(),"cEvCuts") || !strcmp(canvas->GetName(),"cEvCuts2") ) {
    switch(iLayer) {
    case 0 : hist->SetLineColor(kRed); break;;
    case 1 : hist->SetLineColor(kBlue); break;;
    case 2 : hist->SetLineColor(kBlack); break;;
    case 3 : hist->SetLineColor(kBlack); break;;
    }
  }
  else   if ( !strcmp(canvas->GetName(),"cTrCuts") || !strcmp(canvas->GetName(),"cTrCuts2") ) {
    switch(iLayer) {
    case 0 : hist->SetLineColor(kBlack); break;;
    case 1 : hist->SetLineColor(kBlue); break;;
    case 2 : hist->SetLineColor(kGreen); break;;
    case 3 : hist->SetLineColor(kRed); break;;
    case 4 : hist->SetLineColor(kBlue-6); break;;
    case 5 : hist->SetLineColor(kGreen-6); break;;
    case 6 : hist->SetLineColor(kRed-6); break;;
    }
    
  }
  else {
    switch(iLayer) {
    case 0 : hist->SetLineColor(kBlack); break;;
    case 1 : hist->SetLineColor(kBlue); break;;
    case 2 : hist->SetLineColor(kGreen); break;;
    case 3 : hist->SetLineColor(kRed); break;;
    case 4 : hist->SetLineColor(kOrange); break;;
    case 5 : hist->SetLineColor(kBlue-6); break;;
    case 6 : hist->SetLineColor(kGreen-6); break;;
    case 7 : hist->SetLineColor(kRed-6); break;;
    }
  }
  
  hist->SetStats(kFALSE);

  if ( !iLayer )
    hist->Draw("E");
  else 
    hist->Draw("E,SAME");
  
  return;
}
