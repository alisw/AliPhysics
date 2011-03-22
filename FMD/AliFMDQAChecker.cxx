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
    
    TH1F* hist  = 0;
    Int_t nHist = list[specie]->GetEntriesFast();
    for(Int_t i= 0; i< nHist; i++) {
      
      if (!(hist = static_cast<TH1F*>(list[specie]->At(i)))) continue;
      
      if(what == AliQAv1::kESD) 
	rv[specie] += (hist->GetMean() > 0 ? 1 : 0);
      if(what == AliQAv1::kRAW) 
	rv[specie] += (hist->GetMean() > 0 ? 1 : 0);
      if(what == AliQAv1::kSIM) 
	rv[specie] += (hist->GetMean() > 0 ? 1 : 0);
      if(what == AliQAv1::kREC) 
	rv[specie] += (hist->GetMean() > 0 ? 1 : 0);
    } // for (int i ...)
    // if (count != 0) rv[specie] /= count;
  }
  // return rv;
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
    if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue ;
    
    if(!list[specie]) continue;

    TH1F* hist  = 0;
    Int_t nHist = list[specie]->GetEntriesFast();
    for(Int_t i= 0; i< nHist; i++) {
      hist = static_cast<TH1F*>(list[specie]->At(i));
      if (hist && hist->TestBit(AliQAv1::GetImageBit())) {
        nImages++; 
	max = TMath::Max(max, hist->GetMaximum());
	min = TMath::Min(min, hist->GetMinimum());
      }
    }
    break ; 
  }
  min = TMath::Max(0.1, min);
  max = TMath::Min(1.0, max);

  if (nImages == 0) {
    AliDebug(AliQAv1::GetQADebugLevel(), 
	     Form("No histogram will be plotted for %s %s\n", GetName(), 
		  AliQAv1::GetTaskName(task).Data()));
    return;
  }

  AliDebug(AliQAv1::GetQADebugLevel(), 
	   Form("%d histograms will be plotted for %s %s\n", 
		nImages, GetName(), AliQAv1::GetTaskName(task).Data()));  

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (!AliQAv1::Instance()->IsEventSpecieSet(specie)) continue ;
    
    if(!list[specie]) continue;

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
    fImage[specie]->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), 
				AliQAv1::GetModeName(mode), 
				AliQAChecker::Instance()->GetRunNumber(), 
				AliQAv1::GetImageFileFormat()), "ps") ; 
    fImage[specie]->Clear(); 

    Int_t nx = int(nImages + .5) / 2;
    Int_t ny = 2;
    fImage[specie]->Divide(nx, ny, 0, 0);
    
    
    TH1F* hist  = 0;
    Int_t nHist = list[specie]->GetEntriesFast();
    Int_t j     = 0;
    for (Int_t i = 0; i < nHist; i++) { 
      hist = static_cast<TH1F*>(list[specie]->At(i));
      if (!(hist && hist->TestBit(AliQAv1::GetImageBit()))) continue;
      
      TVirtualPad* pad = fImage[specie]->cd(++j);
      pad->SetLogy();
      hist->SetMinimum(min);
      hist->SetMaximum(max);
      hist->DrawCopy();
    }

    fImage[specie]->Print(Form("%s%s%d.%s", AliQAv1::GetImageFileName(), 
			       AliQAv1::GetModeName(mode), 
			       AliQAChecker::Instance()->GetRunNumber(), 
			       AliQAv1::GetImageFileFormat()), "ps"); 
  }
}

//__________________________________________________________________
//
// EOF
//
