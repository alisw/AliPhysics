/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Checks the quality assurance.                                         //
//  By comparing with reference data                                      //
//  S.Radomski Uni-Heidelberg October 2007                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include <TClass.h>
#include <TH1F.h> 
#include <TH1I.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliTRDQAChecker.h"

ClassImp(AliTRDQAChecker)

//__________________________________________________________________

void AliTRDQAChecker::Check(Double_t * test, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam* /*param*/) 
{

  // Super-basic check on the QA histograms on the input list: 

  for(Int_t i=0; i<AliRecoParam::kNSpecies; i++) test[i] = 0.5; 

  //Int_t count[AliRecoParam::kNSpecies] = { 0 }; 

  if (index != AliQAv1::kREC) return;

  const Double_t lowAmp = 30;
  const Double_t highAmp = 50;

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    
    TH1D *hist = (TH1D*) list[specie]->At(12);
    if (!hist) continue;
    
    Double_t value = hist->Integral(hist->FindBin(lowAmp), hist->FindBin(highAmp));
    if (hist->GetSum())
      test[specie] = value / hist->GetSum();

  }
}  

//____________________________________________________________________________
