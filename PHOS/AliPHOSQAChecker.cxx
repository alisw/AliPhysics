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

/*
  Checks the quality assurance. 
  By comparing with reference data
  Y. Schutz CERN July 2007
*/

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
#include "AliPHOSQAChecker.h"

ClassImp(AliPHOSQAChecker)

//__________________________________________________________________

AliPHOSQAChecker & AliPHOSQAChecker::operator = (const AliPHOSQAChecker &)
{
  Fatal("operator =", "not implemented");
  return *this;
}

//____________________________________________________________________________
void AliPHOSQAChecker::Check(Double_t * test, AliQAv1::ALITASK_t task, TObjArray ** list, const AliDetectorRecoParam * /* recoParam */) 
{
  // Performs a basic checking
  // Compares all the histograms in the list

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    test[specie] = 1.0;
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue ; 
    // checking for empty histograms
    // if (CheckEntries(list[specie]) == 0)  {
    //   AliWarning("histograms are empty");
    //   test[specie] = 0.4;//-> Corresponds to kWARNING see AliQACheckerBase::Run
    // }
  
    // checking raw data
    if(task == AliQAv1::kRAW){
      if(AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCalib    ||
         AliRecoParam::ConvertIndex(specie) == AliRecoParam::kHighMult ||
         AliRecoParam::ConvertIndex(specie) == AliRecoParam::kLowMult  ||
	 AliRecoParam::ConvertIndex(specie) == AliRecoParam::kDefault) {
	// list[specie]->Print();
	TH1F *hHighNtot = (TH1F*)list[specie]->At(13);
	if (hHighNtot!=0) {
	  if (hHighNtot->GetMean() < 1000) test[specie]=1;
	}
	else test[specie]=0.1;
      }
    }

    //default check response. It will be changed when reasonable checks will be considered
    else test[specie] = 0.7 ; // /-> Corresponds to kINFO see AliQACheckerBase::Run 
  } // species loop
}
