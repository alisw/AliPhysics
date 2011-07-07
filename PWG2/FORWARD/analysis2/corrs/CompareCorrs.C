/**
 * @file   CompareCorrs.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Fri Jan 28 23:01:59 2011
 * 
 * @brief  Utilities for comparing correction objects 
 * 
 * 
 */
#ifndef __CINT__
#include <TCanvas.h>
#include <TFile.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include "AliForwardCorrectionManager.h"
#endif


//======================================================================
struct Canvas 
{
  //____________________________________________________________________
  Canvas(const char* name, const char* title,
	 const char* n1,   const char* n2)
    : fName(name),
      fTitle(title),
      fN1(n1), 
      fN2(n2),
      fCanvas(0),
      fBody(0)
  {
    gStyle->SetPalette(1);
    gStyle->SetTitleX(.10);
    gStyle->SetTitleY(.99);
    gStyle->SetTitleW(.85);
    gStyle->SetTitleH(.085);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTitleBorderSize(0);
  }
  //____________________________________________________________________
  void Open() 
  {
    fCanvas = new TCanvas(fName, fTitle, 800, TMath::Sqrt(2)*800);
    fCanvas->SetFillColor(0);
  
    fCanvas->Print("comparison.pdf[", "pdf");
  }
  //____________________________________________________________________
  TPad* 
  Clear(UShort_t nPad, UShort_t d, Char_t r)
  {
    fCanvas->Clear();
    TPad* top = new TPad("top", "Top", 0, .95, 1, 1, 0, 0);
    top->Draw();
    top->cd();

    TLatex* l = new TLatex(.5, .5, Form("%s for FMD%d%c (%s / %s)", 
					fTitle, d, r, fN1, fN2));
    l->SetNDC();
    l->SetTextAlign(22);
    l->SetTextSize(0.3);
    l->Draw();
  
    fCanvas->cd();
    fBody = new TPad("body", "Body", 0, 0, 1, .95, 0, 0);
    fBody->SetTopMargin(0.05);
    fBody->SetRightMargin(0.05);
    fBody->Divide(2, (nPad+1)/2, 0.001, 0.001);
    fBody->Draw();
    
    return fBody;
  }  
  //____________________________________________________________________
  TVirtualPad* cd(Int_t i) 
  {
    if (!fBody) return 0;
    
    return fBody->cd(i);
  }
  //____________________________________________________________________
  void Print(UShort_t d, Char_t r, const char* extra="")
  {
    fCanvas->Print("comparison.pdf", 
		   Form("Title:FMD%d%c %s", d, r, extra));
  }
  //____________________________________________________________________
  void Close()
  {
    fCanvas->Print("comparison.pdf]", "pdf");
  }    
  //____________________________________________________________________
  const char* fName;
  const char* fTitle;
  const char* fN1;
  const char* fN2;
  TCanvas*    fCanvas;
  TPad*       fBody;
};

//======================================================================
void
GetObjects(UShort_t    what, 
	   const char* fn1, const char* fn2, 
	   TObject*&   o1,  TObject*&   o2)
{
  // --- Open files --------------------------------------------------
  TFile* file1 = TFile::Open(fn1, "READ");
  TFile* file2 = TFile::Open(fn2, "READ");

  if (!file1) { 
    Error("CompareSecMaps", "File %s cannot be opened", fn1);
    return;
  }

  if (!file2) { 
    Error("CompareSecMaps", "File %s cannot be opened", fn2);
    return;
  }
  
  // --- Find Objects ------------------------------------------------
  AliForwardCorrectionManager::ECorrection ewhat = what;
  // (AliForwardCorrectionManager::ECorrection)what;
  const char* objName = 
    AliForwardCorrectionManager::Instance().GetObjectName(ewhat);
  
  o1 = file1->Get(objName);
  o2 = file2->Get(objName);
  
  if (!o1) {
    Error("CompareSecMaps", "File %s does not contain an object named %s", 
	  fn1, objName);
    return;
  }
  if (!o2) {
    Error("CompareSecMaps", "File %s does not contain an object named %s", 
	  fn2, objName);
    return;
  }
};



