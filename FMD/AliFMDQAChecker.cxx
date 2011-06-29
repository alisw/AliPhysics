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
#include <TH1I.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 
#include <iostream>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TStyle.h>

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
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    AliRecoParam::EventSpecie_t spe = AliRecoParam::ConvertIndex(specie);
    if (!AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))
	->IsEventSpecieSet(spe)) 
      continue;
    // if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue ;
    
    if(!list[specie] || list[specie]->GetEntriesFast() == 0) continue;

    TH1* hist  = 0;
    Int_t nHist = list[specie]->GetEntriesFast();
    for(Int_t i= 0; i< nHist; i++) {
      hist = static_cast<TH1F*>(list[specie]->At(i));
      if (hist && hist->TestBit(AliQAv1::GetImageBit())) {
        nImages++; 
	TString name(hist->GetName());
	if (name.Contains("readouterrors", TString::kIgnoreCase)) continue;
	Double_t hMax = hist->GetMaximum(); 
	// hist->GetBinContent(hist->GetMaximumBin());
	Double_t hMin = hist->GetMinimum();
	// hist->GetBinContent(hist->GetMinimumBin());
	max = TMath::Max(max, hMax);
	min = TMath::Min(min, hMin);
	AliInfo(Form("Histogram %30s min=%f, max=%f (min=%f,max=%f)",
		     name.Data(), hMin, hMax, min, max));
      }
    }
    break ; 
  }
  min = TMath::Max(0.1, min);
  max = TMath::Max(1.0, max);

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
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (!AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))
	->IsEventSpecieSet(specie)) continue;
    
    if(!list[specie] || list[specie]->GetEntries() <= 0) continue;

    const Char_t * title = Form("QA_%s_%s_%s", GetName(), 
				AliQAv1::GetTaskName(task).Data(), 
				AliRecoParam::GetEventSpecieName(specie)); 
    if (!fImage[specie]) 
      fImage[specie] = new TCanvas(title, title) ;
    fImage[specie]->Clear() ; 
    fImage[specie]->SetTitle(title) ; 
    fImage[specie]->cd() ; 

    TPaveText someText(0.015, 0.015, 0.98, 0.98) ;
    someText.AddText(title) ;
    someText.Draw() ; 
    TString outName(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), 
			 AliQAv1::GetModeName(mode), 
			 AliQAChecker::Instance()->GetRunNumber(), 
			 AliQAv1::GetImageFileFormat()));
    fImage[specie]->Print(outName, "ps") ; 
    fImage[specie]->Clear(); 

    Int_t nx = int(nImages + .5) / 2;
    Int_t ny = 2;
    fImage[specie]->Divide(nx, ny, 0, 0);
    
    
    TH1*  hist  = 0;
    Int_t nHist = list[specie]->GetEntriesFast();
    // list[specie]->Print();
    Int_t j     = 0;
    for (Int_t i = 0; i < nHist; i++) { 
      hist = static_cast<TH1*>(list[specie]->At(i));
      if (!hist || !hist->TestBit(AliQAv1::GetImageBit())) continue;
      
      TVirtualPad* pad = fImage[specie]->cd(++j);
      pad->SetRightMargin(0.001);
      Int_t logOpts = 0;
      logOpts |= CheckForLog(hist->GetXaxis(), pad, 1);
      logOpts |= CheckForLog(hist->GetYaxis(), pad, 2);
      logOpts |= CheckForLog(hist->GetZaxis(), pad, 3);
      if (hist->GetEntries() > 0)
	AliInfo(Form("Drawing %30s (min=%f, max=%f) in pad %d",
		     hist->GetName(), min, max, j));

      TString opt;
      TString name(hist->GetName());
      if (name.Contains("readouterrors", TString::kIgnoreCase)) {
	pad->SetRightMargin(0.15);
	pad->SetBottomMargin(0.10);
	pad->SetTopMargin(0.02);
	opt="COLZ";
      }
      else {
	hist->SetMinimum(min);
	hist->SetMaximum(max);

      }
      hist->DrawCopy(opt);

      if (name.Contains("readouterrors", TString::kIgnoreCase)) 
	continue;

      TPad* insert = new TPad("insert", "Zoom", 
			      .4,.4, .99, .95, 0, 0, 0);
      insert->SetTopMargin(0);
      insert->SetRightMargin(0.001);
      insert->SetFillColor(0);
      insert->SetBorderSize(1);
      insert->SetBorderMode(0);
      insert->Draw();
      insert->cd();
      if (logOpts & 0x1) insert->SetLogx();
      if (logOpts & 0x2) insert->SetLogy();
      if (logOpts & 0x4) insert->SetLogz();
      hist->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetXmin(), 
				     hist->GetXaxis()->GetXmax()/8);
      hist->DrawCopy(opt);
      pad->cd();
    }

    fImage[specie]->Print(outName, "ps");
  }
}

//__________________________________________________________________
//
// EOF
//
