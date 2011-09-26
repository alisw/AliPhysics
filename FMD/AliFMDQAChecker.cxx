/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//__________________________________________________________________
//
// Yves?
// What 
// is 
// this 
// class 
// supposed 
// to
// do?
//__________________________________________________________________
//
// --- ROOT system ---
#include <TClass.h>
#include <TH1F.h> 
#include <TH2.h> 
#include <TH1I.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 
#include <iostream>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TLatex.h>

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliFMDQAChecker.h"
#include "AliRecoParam.h"

ClassImp(AliFMDQAChecker)
#if 0
; // This is for Emacs! - do not delete
#endif

//__________________________________________________________________
Double_t 
AliFMDQAChecker::CheckOne(AliQAv1::ALITASK_t          what,
			  AliRecoParam::EventSpecie_t specie, 
			  TH1*                        hist) const
{
  if(what == AliQAv1::kESD) return CheckESD(specie, hist);
  if(what == AliQAv1::kRAW) return CheckRaw(specie, hist);
  if(what == AliQAv1::kSIM) return CheckSim(specie, hist);
  if(what == AliQAv1::kREC) return CheckRec(specie, hist);
  return 0;
}
//__________________________________________________________________
Double_t 
AliFMDQAChecker::CheckESD(AliRecoParam::EventSpecie_t /* specie*/, 
			  TH1*                        hist) const
{
  return (hist->GetMean() > 0 ? 1 : 0);
}
//__________________________________________________________________
Double_t 
AliFMDQAChecker::CheckRaw(AliRecoParam::EventSpecie_t /* specie*/, 
			  TH1*                        hist) const
{
  return (hist->GetMean() > 0 ? 1 : 0);
}
//__________________________________________________________________
Double_t 
AliFMDQAChecker::CheckSim(AliRecoParam::EventSpecie_t /* specie*/, 
			  TH1*                        hist) const
{
  return (hist->GetMean() > 0 ? 1 : 0);
}
//__________________________________________________________________
Double_t 
AliFMDQAChecker::CheckRec(AliRecoParam::EventSpecie_t /* specie*/, 
			  TH1*                        hist) const
{
  return (hist->GetMean() > 0 ? 1 : 0);
}

//__________________________________________________________________
void AliFMDQAChecker::Check(Double_t*                   rv, 
			    AliQAv1::ALITASK_t          what, 
			    TObjArray**                 list, 
			    const AliDetectorRecoParam* /*t*/) 
{
  // 
  // Member function called to do the actual checking
  //
  // Parameters: 
  //    rv   Array of return values. 
  //    what What to check 
  //    list Array of arrays of histograms.  There's one arrat for
  //         each 'specie'
  //    t    Reconstruction parameters - not used. 
  //
  
  // Double_t* rv = new Double_t[AliRecoParam::kNSpecies] ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    // Int_t count   = 0;
    rv[specie]    = 0.; 

    if (!AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue ;
    
    if(!list[specie]) continue;
    
    TH1* hist  = 0;
    Int_t nHist = list[specie]->GetEntriesFast();
    for(Int_t i= 0; i< nHist; i++) {
      
      if (!(hist = static_cast<TH1*>(list[specie]->At(i)))) continue;
      
      rv[specie] += CheckOne(what, AliRecoParam::ConvertIndex(specie), hist);
    } // for (int i ...)
    // if (count != 0) rv[specie] /= count;
  }
  // return rv;
}

namespace {
  Int_t CheckForLog(TAxis*       axis,
		    TVirtualPad* pad, 
		    Int_t        xyz)
  {
    Int_t ret = 0;
    TString t(axis->GetTitle());
    if (!t.Contains("[log]", TString::kIgnoreCase)) return 0;
    t.ReplaceAll("[log]", "");
    switch (xyz) { 
    case 1: pad->SetLogx(); ret |= 0x1; break;
    case 2: pad->SetLogy(); ret |= 0x2; break;
    case 3: pad->SetLogz(); ret |= 0x4; break;
    }
    axis->SetTitle(t);
    return ret;
  }
  void RestoreLog(TAxis* axis, Bool_t log) 
  {
    if (!log) return;
    TString t(axis->GetTitle());
    t.Append("[log]");
    axis->SetTitle(t);
  }
}

namespace {
  void FindMinMax(TH1* h, Double_t& min, Double_t& max)
  {
    Double_t tmin = 1e9;
    Double_t tmax = 0;
    for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
      Double_t c = h->GetBinContent(i);
      if (c < 1e-8) continue;
      tmin = TMath::Min(tmin, c);
      tmax = TMath::Max(tmax, c);
    }
    min = tmin;
    max = tmax;
  }
}
//____________________________________________________________________________ 
void 
AliFMDQAChecker::MakeImage(TObjArray** list, 
			   AliQAv1::TASKINDEX_t task, 
			   AliQAv1::MODE_t mode) 
{
  // makes the QA image for sim and rec
  // 
  // Parameters: 
  //    task What to check 
  //    list Array of arrays of histograms.  There's one array for
  //         each 'specie'
  //    t    Reconstruction parameters - not used. 
  // 
  Int_t    nImages = 0 ;
  Double_t max     = 0;
  Double_t min     = 10000;

  // Loop over all species 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    AliRecoParam::EventSpecie_t spe = AliRecoParam::ConvertIndex(specie);
    if (!AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))
	->IsEventSpecieSet(spe)) 
      continue;
									
    // If nothing is defined for this specie, go on. 
    if(!list[specie] || list[specie]->GetEntriesFast() == 0) continue;

    // Loop over the histograms and figure out how many histograms we
    // have and the min/max 
    TH1* hist  = 0;
    Int_t nHist = list[specie]->GetEntriesFast();
    for(Int_t i= 0; i< nHist; i++) {
      hist = static_cast<TH1F*>(list[specie]->At(i));
      if (hist && hist->TestBit(AliQAv1::GetImageBit())) {
        nImages++; 
	TString name(hist->GetName());
	if (name.Contains("readouterrors", TString::kIgnoreCase)) continue;

	// Double_t hMax = hist->GetMaximum(); 
	// hist->GetBinContent(hist->GetMaximumBin());
	// Double_t hMin = hist->GetMinimum();
	// hist->GetBinContent(hist->GetMinimumBin());
	Double_t hMax, hMin;
	FindMinMax(hist, hMin, hMax);
	max = TMath::Max(max, hMax);
	min = TMath::Min(min, hMin);
      }
    }
    break ; 
  }
  min = TMath::Max(1e-6, min);
  max = TMath::Max(1.0,  max);

  // IF no images, go on. 
  if (nImages == 0) {
    AliDebug(AliQAv1::GetQADebugLevel(), 
	     Form("No histogram will be plotted for %s %s\n", GetName(), 
		  AliQAv1::GetTaskName(task).Data()));
    return;
  }

  AliDebug(AliQAv1::GetQADebugLevel(), 
	   Form("%d histograms will be plotted for %s %s\n", 
		nImages, GetName(), AliQAv1::GetTaskName(task).Data()));  
  gStyle->SetOptStat(0);
  
  // Again loop over species and draw a canvas 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (!AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))
	->IsEventSpecieSet(specie)) continue;

    // if Nothing here, go on
    if(!list[specie] || list[specie]->GetEntries() <= 0 || 
       nImages <= 0) continue;

    // Form the title 
    const Char_t * title = Form("QA_%s_%s_%s", GetName(), 
				AliQAv1::GetTaskName(task).Data(), 
				AliRecoParam::GetEventSpecieName(specie)); 
    if (!fImage[specie]) fImage[specie] = new TCanvas(title, title) ;
    fImage[specie]->Clear() ; 
    fImage[specie]->SetTitle(title) ; 
    fImage[specie]->cd() ; 

    // Put something in the canvas - even if empty 
    TPaveText someText(0.015, 0.015, 0.98, 0.98) ;
    someText.AddText(title) ;
    someText.SetFillColor(0);
    someText.SetFillStyle(0);
    someText.SetBorderSize(0);
    someText.SetTextColor(kRed+1);
    someText.Draw() ; 
    TString outName(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), 
			 AliQAv1::GetModeName(mode), 
			 AliQAChecker::Instance()->GetRunNumber(), 
			 AliQAv1::GetImageFileFormat()));
    fImage[specie]->Print(outName, "ps") ; 


    // Now set some parameters on the canvas 
    fImage[specie]->Clear(); 
    fImage[specie]->SetTopMargin(0.10);
    fImage[specie]->SetBottomMargin(0.15);
    fImage[specie]->SetLeftMargin(0.15);
    fImage[specie]->SetRightMargin(0.05);
    
    // Put title on top 
    const char* topT = Form("Mode: %s, Task: %s, Specie: %s, Run: %d",
			    AliQAv1::GetModeName(mode), 
			    AliQAv1::GetTaskName(task).Data(), 
			    AliRecoParam::GetEventSpecieName(specie),
			    AliQAChecker::Instance()->GetRunNumber());
    TLatex* topText = new TLatex(.5, .99, topT);
    topText->SetTextAlign(23);
    topText->SetTextSize(.038);
    topText->SetTextFont(42);
    topText->SetTextColor(kBlue+3);
    topText->SetNDC();
    topText->Draw();
				 
    // Divide canvas 
    Int_t nx = int(nImages + .5) / 2;
    Int_t ny = 2;
    fImage[specie]->Divide(nx, ny, 0, 0);
    
    
    // Loop over histograms 
    TH1*  hist  = 0;
    Int_t nHist = list[specie]->GetEntriesFast();
    Int_t j     = 0;
    for (Int_t i = 0; i < nHist; i++) { 
      hist = static_cast<TH1*>(list[specie]->At(i));
      if (!hist || !hist->TestBit(AliQAv1::GetImageBit())) continue;

      // Go to sub-pad 
      TVirtualPad* pad = fImage[specie]->cd(++j);
      pad->SetRightMargin(0.01);

      // Check for log scale 
      Int_t logOpts = 0;
      logOpts |= CheckForLog(hist->GetXaxis(), pad, 1);
      logOpts |= CheckForLog(hist->GetYaxis(), pad, 2);
      logOpts |= CheckForLog(hist->GetZaxis(), pad, 3);

      // Figure out special cases 
      TString opt;
      TString name(hist->GetName());
      if (name.Contains("readouterrors", TString::kIgnoreCase)) {
	pad->SetRightMargin(0.15);
	pad->SetBottomMargin(0.10);
	pad->SetTopMargin(0.02);
	opt="COLZ";
      }
      else {
	pad->SetGridx();
	pad->SetGridy();
	hist->SetMinimum(min);
	hist->SetMaximum(max);

      }
      // Draw (As a copy)
      hist->DrawCopy(opt);

      // Special cases 
      if (name.Contains("readouterrors", TString::kIgnoreCase)) {
	for (Int_t kk = 1; kk <= 3; kk++) {
	  TH1* proj = static_cast<TH2*>(hist)->ProjectionY("",kk,kk);
	Double_t m = proj->GetMean(); 
	TLatex* l = new TLatex(kk, 30, Form("Mean: %f", m));
	l->SetTextAngle(90);
	l->SetTextColor(m > 10 ? kRed+1 : m > .7 ? kOrange+2 :kGreen+2);
	l->Draw();
      }
      }
      else {
	gStyle->SetOptTitle(0);
	TPad* insert = new TPad("insert", "Zoom", 
				.4,.4, .99, .95, 0, 0, 0);
	insert->SetTopMargin(0.01);
	insert->SetRightMargin(0.01);
	insert->SetFillColor(0);
	insert->SetBorderSize(1);
	insert->SetBorderMode(0);
	insert->Draw();
	insert->cd();
	if (logOpts & 0x1) insert->SetLogx();
	if (logOpts & 0x2) insert->SetLogy();
	if (logOpts & 0x4) insert->SetLogz();
	hist->GetXaxis()->SetRange(1, hist->GetNbinsX()/8);
	TH1* copy = hist->DrawCopy(opt);
	copy->GetXaxis()->SetNdivisions(408, false);
	// Restore full range 
	hist->GetXaxis()->SetRange(0, 0);
	gStyle->SetOptTitle(1);
      }
      pad->cd();
      // Possibly restore the log options 
      RestoreLog(hist->GetXaxis(), logOpts & 0x1);
      RestoreLog(hist->GetYaxis(), logOpts & 0x2);
      RestoreLog(hist->GetZaxis(), logOpts & 0x4);
    }
    // Print to a post-script file 
    fImage[specie]->Print(outName, "ps");
  }
}

//__________________________________________________________________
//
// EOF
//
