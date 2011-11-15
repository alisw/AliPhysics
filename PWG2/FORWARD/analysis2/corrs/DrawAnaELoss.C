/**
 * Script to draw the energy loss fits from the output file of 
 * AliFMDELossFitter(Task). 
 * 
 *
 * @ingroup pwg2_forward_scripts_corr
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
bool landscape = true;

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
 * @ingroup pwg2_forward_scripts_corr
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
 * @ingroup pwg2_forward_scripts_corr
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
 * @ingroup pwg2_forward_scripts_corr
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

  Int_t w = Int_t(800 / TMath::Sqrt(2));
  Int_t h = 800;
  if (landscape) { 
    Int_t tmp = h;
    h         = w;
    w         = tmp;
  }
  canvas = new TCanvas("c", "C", w, h);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(0);
  canvas->SetBottomMargin(0.15);
  return canvas;
}

//____________________________________________________________________
void CleanStack(THStack* stack)
{
  TIter next(stack->GetHists());
  TObject* o = 0;
  while ((o = next())) { 
    TString name(o->GetName());
    if (name.Contains("_t_")) 
      stack->RecursiveRemove(o);
  }
}
  

//____________________________________________________________________
THStack*
AddToStack(TList* stacks, TList* list, const char* name)
{
  TObject* o = list->FindObject(name);
  if (!o) { 
    Warning("AddToStack", "Object %s not found in %s", name, 
	    list->GetName());
    // list->ls();
    return 0;
  }
  THStack* toAdd = static_cast<THStack*>(o);
  CleanStack(toAdd);
  Info("AddToStack", "Adding %s to stacks", name);
  stacks->Add(toAdd);
  return toAdd;
}

  
//____________________________________________________________________
/** 
 * Draw summary 
 * 
 * @param fname 
 * 
 * @ingroup pwg2_forward_scripts_corr
 */
void DrawSummary(const char* fname="forward_eloss.root", 
		 bool onlySummary=true)
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

  TList    stacks;
  /* THStack* chi2nu  = */ AddToStack(&stacks, fitter, "chi2");
  /* THStack* c       = */ AddToStack(&stacks, fitter, "c");
  /* THStack* delta   = */ AddToStack(&stacks, fitter, "delta");
  /* THStack* xi      = */ AddToStack(&stacks, fitter, "xi");
  /* THStack* sigma   = */ AddToStack(&stacks, fitter, "sigma");
  /* THStack* sigman  = */ AddToStack(&stacks, fitter, "sigman");
  /* THStack* n       = */ AddToStack(&stacks, fitter, "n");
  Int_t    baseA   = stacks.GetEntries()+1;
  Int_t    i       = 2;
  while (true) { 
    if (!AddToStack(&stacks, fitter, Form("a%d",i++)))
      break;
  }
  // stacks.ls();
  Int_t nMax = stacks.GetEntries();
  for (i = nMax-1; i >= baseA; i--) { 
    THStack* stack   = static_cast<THStack*>(stacks.At(i));
    TIter    nextH(stack->GetHists());
    TH1*     hist    = 0;
    Bool_t   hasData = kFALSE;
    while ((hist = static_cast<TH1*>(nextH()))) 
      if (hist->Integral() > 0) hasData = kTRUE;
    if (!hasData) nMax--;
  }

  canvas->SetRightMargin(0.01);
  canvas->SetTopMargin(0.01);
  Int_t nL = (nMax+1) / 2;
  Int_t nX = 2;
  Int_t nY = nL;
  if (landscape) { 
    Int_t tmp = nY;
    nY        = nX;
    nX        = tmp;
  }

  canvas->Divide(nX, nY, 0.1, 0, 0);

  TIter next(&stacks);
  THStack* stack = 0;
  i = 0;
  // Int_t b = 1;
  while ((stack = static_cast<THStack*>(next()))) {
    if (i > nMax) break;
    Int_t ipad = 1+i/nL + 2 * (i % nL);
    Info("DrawSummary", "cd'ing to canvas %d for %s", ipad, 
	 stack->GetName());
    TVirtualPad* p = canvas->cd(ipad);
    p->SetLeftMargin(.6/nL);
    p->SetTopMargin(.01);
    p->SetRightMargin(.01);
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
    if (i == 0) p->SetLogy();
    yaxis->SetTitleSize(0.3/nL);
    yaxis->SetLabelSize(0.08);
    yaxis->SetTitleOffset(2.5/nL);
    yaxis->SetNdivisions(5);
    yaxis->SetTitleFont(42);
    yaxis->SetLabelFont(42);
    yaxis->SetDecimals();
    
    TAxis* xaxis = stack->GetHistogram()->GetXaxis();
    xaxis->SetTitleSize(0.3/nL);
    xaxis->SetLabelSize(0.08);
    xaxis->SetTitleOffset(2./nL);
    xaxis->SetNdivisions(10);
    xaxis->SetTitleFont(42);
    xaxis->SetLabelFont(42);
    xaxis->SetDecimals();

    // Redraw 
    stack->Draw("nostack");
    i++;
    // if (i >= 5) b = 2;
    p->cd();
  }
  canvas->SaveAs("fit_results.png");
  if (!onlySummary) 
    canvas->Print(pdfName, "Title:Fit summary");
}

//____________________________________________________________________
/** 
 * Draw ring fits 
 * 
 * @param fname 
 * 
 * @ingroup pwg2_forward_scripts_corr
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
 * @ingroup pwg2_forward_scripts_corr
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
      
      Info("DrawEtaBins", "Drawing for FMD%d%c", d, r);
      TIter next(edists);
      TH1*  dist = 0;
      Int_t i    = 0;
      Int_t j    = 1;
      while ((dist = static_cast<TH1*>(next()))) { 
	Info("DrawEtaBins", "FMD%d%c: %s", d, r, dist->GetName());
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
 * @ingroup pwg2_forward_scripts_corr
 */
void
DrawAnaELoss(const char* fname="forward_eloss.root", bool onlySummary=true)
{
  if (!CheckCanvas()) {
    Error("DrawFits", "No canvas");
    return;
  }
  if (!onlySummary) canvas->Print(Form("%s[", pdfName));
  DrawSummary(fname, onlySummary);
  if (onlySummary) return;
  DrawRings(fname);
  DrawEtaBins(fname);
  canvas->Print(Form("%s]", pdfName));
}
//
// EOF
//
