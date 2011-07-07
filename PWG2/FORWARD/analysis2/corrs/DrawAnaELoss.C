/**
 * Script to draw the energy loss fits from the output file of 
 * AliFMDELossFitter(Task). 
 * 
 *
 * @ingroup pwg2_forward_analysis_scripts_corr
 */
#include <TFile.h>
#include <THStack.h>
#include <TList.h>
#include <TError.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLegend.h>
#include <TMath.h>

//____________________________________________________________________
// Global 
TList* fitter = 0;
TCanvas* canvas = 0;
const char* pdfName = "FitResults.pdf";

//____________________________________________________________________
/** 
 * Open a file.  The file is expected to contain the directory 
 * structure 
 *
 * @verbatim 
 *  file
 *   +- ForwardResults
 *       +- fmdEnergyFitter 
 *           +- chi2   (THStack)
 *           +- c      (THStack)
 *           +- delta  (THStack)
 *           +- xi     (THStack)
 *           +- sigma  (THStack)
 *           +- sigman (THStack)
 *           +- n      (THStack)
 *           +- a2     (THStack)
 *           +- ...    (THStack)
 *           +- an     (THStack)
 *           +- FMD1I (TList)
 *           |   +- FMD1I_edist (TH1)
 *           |   +- EDists      (TList)
 *           ...
 * @endverbatim
 * 
 * @param fname File to open  
 * 
 * @return Pointer to the list of objects 
 * 
 * @ingroup pwg2_forward_analysis_scripts_corr
 */
TList* OpenFile(const char* fname)
{
  TFile* file = TFile::Open(fname, "READ");
  if (!file) {
    Error("DrawFits", "Couldn't open %s", fname);
    return 0;
  }
    
  TList* forward = static_cast<TList*>(file->Get("ForwardResults"));
  // static_cast<TList*>(file->Get("PWG2forwardDnDeta/Forward"));
  if (!forward) { 
    Error("DrawFits", "Couldn't get forward list from %s", fname);
    return 0;
  }

  fitter = static_cast<TList*>(forward->FindObject("fmdEnergyFitter"));
  if (!fitter) { 
    Error("DrawFits", "Couldn't get fitter folder");
    return 0;
  }
  return fitter;
}
//____________________________________________________________________
/** 
 * Open file if not done already 
 * 
 * @param fname File to open
 * 
 * @return List of fits 
 * 
 * @ingroup pwg2_forward_analysis_scripts_corr
 */
TList* CheckFitter(const char* fname="AnalysisResults.root")
{
  if (!fitter) return OpenFile(fname);
  return fitter;
}

//____________________________________________________________________
/** 
 * Make canvas if not done already 
 * 
 * 
 * @return Canvas 
 * 
 * @ingroup pwg2_forward_analysis_scripts_corr
 */
TCanvas* CheckCanvas()
{
  if (canvas) return canvas;

  gStyle->SetOptFit(111111);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0);
  gStyle->SetTitleY(.9);
  gStyle->SetTitleW(.4);
  gStyle->SetTitleH(.1);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.2);
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(1111);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.15);
  gStyle->SetStatColor(0);
  gStyle->SetStatBorderSize(1);

  canvas = new TCanvas("c", "C", Int_t(800 / TMath::Sqrt(2)), 800);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetBottomMargin(0.15);
  return canvas;
}
  
//____________________________________________________________________
/** 
 * Draw summary 
 * 
 * @param fname 
 * 
 * @ingroup pwg2_forward_analysis_scripts_corr
 */
void DrawSummary(const char* fname="forward_eloss.root")
{
  if (!CheckFitter(fname)) {
    Error("DrawFits", "File not opened");
    return;
  }
  if (!CheckCanvas()) {
    Error("DrawFits", "No canvas");
    return;
  }
  canvas->Clear();

  THStack* chi2nu;
  THStack* c;
  THStack* delta;
  THStack* xi;
  THStack* sigma;
  THStack* sigman;
  THStack* n;
  TList stacks;
  stacks.Add(chi2nu = static_cast<THStack*>(fitter->FindObject("chi2")));
  stacks.Add(c      = static_cast<THStack*>(fitter->FindObject("c")));
  stacks.Add(delta  = static_cast<THStack*>(fitter->FindObject("delta")));
  stacks.Add(xi     = static_cast<THStack*>(fitter->FindObject("xi")));
  stacks.Add(sigma  = static_cast<THStack*>(fitter->FindObject("sigma")));
  stacks.Add(sigman = static_cast<THStack*>(fitter->FindObject("sigman")));
  stacks.Add(n      = static_cast<THStack*>(fitter->FindObject("n")));
  Int_t baseA = stacks.GetEntries()+1;
  Int_t i=2;
  while (true) { 
    TObject* o = fitter->FindObject(Form("a%d",i++));
    if (!o) break;
    Info("DrawFits", "Adding %s", o->GetName());
    stacks.Add(o);
  }
  // stacks.ls();
  Int_t nMax = stacks.GetEntries();
  for (Int_t i = nMax-1; i >= baseA; i--) { 
    THStack* stack   = static_cast<THStack*>(stacks.At(i));
    TIter    nextH(stack->GetHists());
    TH1*     hist    = 0;
    Bool_t   hasData = kFALSE;
    while ((hist = static_cast<TH1*>(nextH()))) 
      if (hist->Integral() > 0) hasData = kTRUE;
    if (!hasData) nMax--;
  }

  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.05);
  canvas->Divide(2, (nMax+1)/2, 0.1, 0, 0);

  TIter next(&stacks);
  THStack* stack = 0;
  i = 0;
  Int_t b = 1;
  while ((stack = static_cast<THStack*>(next()))) {
    if (i > nMax) break;
    TVirtualPad* p = canvas->cd(1+i/5 + 2*(i%5));
    p->SetLeftMargin(.15);
    p->SetFillColor(0);
    p->SetFillStyle(0);
    p->SetGridx();
    stack->Draw("nostack");
    stack->GetHistogram()->SetYTitle(stack->GetTitle());
    stack->GetHistogram()->SetXTitle("#eta");

    TAxis* yaxis = stack->GetHistogram()->GetYaxis();
    if (i == 0)                     yaxis->SetRangeUser(0,20); // Chi2
    if (i == 1)                     stack->SetMaximum(1);      // c
    if (i == 2)                     stack->SetMaximum(1);      // delta
    if (i == 3)                     stack->SetMaximum(0.1);   // xi
    if (i == 4 || i == 5)           stack->SetMaximum(0.5);    // sigma{,n}
    if (i == 7)                     stack->SetMaximum(0.5);    // a
    yaxis->SetTitleSize(0.15);
    yaxis->SetLabelSize(0.08);
    yaxis->SetTitleOffset(0.35);
    yaxis->SetNdivisions(5);

    TAxis* xaxis = stack->GetHistogram()->GetXaxis();
    xaxis->SetTitleSize(0.15);
    xaxis->SetLabelSize(0.08);
    xaxis->SetTitleOffset(0.35);
    xaxis->SetNdivisions(10);

    // Redraw 
    stack->Draw("nostack");
    i++;
    if (i >= 5) b = 2;
    p->cd();
  }
  canvas->Print(pdfName, "Title:Fit summary");
}

//____________________________________________________________________
/** 
 * Draw ring fits 
 * 
 * @param fname 
 * 
 * @ingroup pwg2_forward_analysis_scripts_corr
 */
void DrawRings(const char* fname="AnalysisResults.root")
{
  if (!CheckFitter(fname)) {
    Error("DrawFits", "File not opened");
    return;
  }
  if (!CheckCanvas()) {
    Error("DrawFits", "No canvas");
    return;
  }
  canvas->Clear();
  
  canvas->Clear();
  canvas->Divide(1, 5,0,0);

  const char* dets[] = { "FMD1I", "FMD2I", "FMD2O", "FMD3I", "FMD3O", 0 };
  for (Int_t i = 0; i < 5; i++) { 
    TVirtualPad* p = canvas->cd(i+1);
    p->SetGridx();
    p->SetFillColor(0);
    p->SetFillStyle(0);
    p->SetLogy();
    TList* d = static_cast<TList*>(fitter->FindObject(dets[i]));
    if (!d) { 
      Warning("DrawFits", "List %s not found", dets[i]);
      continue;
    }
    TH1* edist = static_cast<TH1*>(d->FindObject(Form("%s_edist", dets[i])));
    if (!edist) {
      Warning("DrawFits", "Histogram %s_edist not found", dets[i]);
      continue;
    }
    edist->Draw();
    TF1*   f = 0;
    TIter  nextF(edist->GetListOfFunctions());
    while ((f = static_cast<TF1*>(nextF()))) {
      Double_t chi2 = f->GetChisquare();
      Int_t    ndf  = f->GetNDF();
      Printf("%s %s:\n  Range: %f-%f\n" 
             "chi^2/ndf= %f / %d = %f", 
	     edist->GetName(), f->GetName(), 
             f->GetXmin(), f->GetXmax(), chi2, ndf, 
	     (ndf > 0) ? chi2/ndf : 0);
      for (Int_t j = 0; j < f->GetNpar(); j++) { 
	Printf("  %-20s : %9.4f +/- %9.4f", 
	       f->GetParName(j), f->GetParameter(j), f->GetParError(j));
      }
    }
    p->cd();
  }
  canvas->cd();
  canvas->Print(pdfName, "Title:Fit to rings");
}

//____________________________________________________________________
/** 
 * Draw fits in eta bins 
 * 
 * @param fname 
 * 
 * @ingroup pwg2_forward_analysis_scripts_corr
 */
void DrawEtaBins(const char* fname="AnalysisResults.root")
{
  if (!CheckFitter(fname)) {
    Error("DrawFits", "File not opened");
    return;
  }
  if (!CheckCanvas()) {
    Error("DrawFits", "No canvas");
    return;
  }
  canvas->Clear();
  canvas->Divide(2,2,0,0);

  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nQ = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nQ; q++) { 
      Char_t r = (q == 0 ? 'I' : 'O');
      
      TList* ring = 
	static_cast<TList*>(fitter->FindObject(Form("FMD%d%c",d,r)));
      if (!ring) { 
	Error("PrintFits", "Couldn't get ring FMD%d%c", d,r);
	continue; 
      }
      TList* edists = static_cast<TList*>(ring->FindObject("EDists"));
      if (!edists) { 
	Error("PrintFits", "Couldn't get EDists list for FMD%d%c", d,r);
	continue; 
      }
      TIter next(edists);
      TH1*  dist = 0;
      Int_t i    = 0;
      Int_t j    = 1;
      while ((dist = static_cast<TH1*>(next()))) { 
	if (i == 4) { 
	  i = 0;
	  j++;
	  canvas->Print(pdfName, Form("Title:FMD%d%c page %2d", d,r,j));
	}
	TVirtualPad* p = canvas->cd(++i);
	p->SetFillColor(kWhite);
	p->SetFillStyle(0);
	p->SetBorderSize(0);
	p->SetLogy();
	dist->SetMaximum(15);
	dist->Draw();
	
      }
      if (i != 0) {
	i++;
	for (; i <= 4; i++) {
	  TVirtualPad* p = canvas->cd(i);
	  p->Clear();
	  p->SetFillColor(kMagenta-3);
	  p->SetFillStyle(0);
	  p->SetBorderSize(0);
	}
	canvas->Print(pdfName, Form("FMD%d%c page %2d", d,r,j++));
      }
    }
  }
}   

//____________________________________________________________________
/** 
 * Draw energy loss fits to a multi-page PDF
 *
 * The input file is the result of running AliFMDELossFitter - 
 * either directly via AliFMDELossFitterTask or as part of a larger 
 * train (AliForwardMultiplicityTask or AliForwardMCMultiplicityTask). 
 * 
 * @verbatim 
 *  file
 *   +- ForwardResults 
 *       +- fmdEnergyFitter 
 *           +- chi2   (THStack)
 *           +- c      (THStack)
 *           +- delta  (THStack)
 *           +- xi     (THStack)
 *           +- sigma  (THStack)
 *           +- sigman (THStack)
 *           +- n      (THStack)
 *           +- a2     (THStack)
 *           +- ...    (THStack)
 *           +- an     (THStack)
 *           +- FMD1I (TList)
 *           |   +- FMD1I_edist (TH1)
 *           |   +- EDists      (TList)
 *           ...
 * @endverbatim
 *
 * @param fname 
 * 
 * @ingroup pwg2_forward_analysis_scripts_corr
 */
void
DrawAnaELoss(const char* fname="forward_eloss.root")
{
  if (!CheckCanvas()) {
    Error("DrawFits", "No canvas");
    return;
  }
  canvas->Print(Form("%s[", pdfName));
  DrawSummary(fname);
  DrawRings(fname);
  DrawEtaBins(fname);
  canvas->Print(Form("%s]", pdfName));
}
//
// EOF
//
