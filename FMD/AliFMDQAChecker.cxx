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
Double_t* AliFMDQAChecker::Check(AliQAv1::ALITASK_t what, TObjArray ** list, AliDetectorRecoParam* /*t*/) {
  
  Double_t* rv = new Double_t[AliRecoParam::kNSpecies] ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    rv[specie] = 0.0 ; 
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue ;
    
    if(!list[specie]) continue;
    
    TIter next1(list[specie]);
    TH1F*  hist = 0;
    
    for(Int_t i= 0; i<list[specie]->GetEntriesFast(); i++) {
      
      hist = (TH1F*)list[specie]->At(i);
      if(!hist) continue;
      
      if(what == AliQAv1::kESD) {
	if(hist->GetMean() > 0)
	  rv[specie] = 1.;
	else rv[specie] = 0.;
      }
      if(what == AliQAv1::kRAW) {
	if(hist->GetMean() > 0)
	  rv[specie] = 1.;
	else rv[specie] = 0.;
      }
      if(what == AliQAv1::kSIM) {
	if(hist->GetMean() > 0)
	  rv[specie] = 1.;
	else rv[specie] = 0.;
      }
      if(what == AliQAv1::kREC) {
	if(hist->GetMean() > 0)
	  rv[specie] = 1.;
	else rv[specie] = 0.;
      }
      
    }
  }
  
  return rv;
  
}



//__________________________________________________________________
//
// EOF
//
