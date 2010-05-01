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

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  Checks the quality assurance.                                  //
//  By analysis of the histograms & comparing with reference data  //
//  S.Arcelli                                                      //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "TH1.h"
#include "TObjArray.h"

#include "AliLog.h"
//#include "AliQAv1.h"
//#include "AliQAChecker.h"

#include "AliTOFQAChecker.h"

ClassImp(AliTOFQAChecker)

//____________________________________________________________________________
void AliTOFQAChecker::Check(Double_t * test, AliQAv1::ALITASK_t /*index*/,
				  TObjArray ** list,
				  const AliDetectorRecoParam * /*recoParam*/) 
{
  // Super-basic check on the QA histograms on the input list: 
  // look whether they are empty!

  Int_t count[AliRecoParam::kNSpecies] = { 0 }; 

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (! AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(specie)) ) 
      continue ;
    test[specie] = 1.0 ; 
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue ; 
    if (list[specie]->GetEntries() == 0){  
      test[specie] = 0.0 ; // nothing to check
    }
    else {
      TIter next(list[specie]) ; 
      TH1 * hdata ;
      count[specie] = 0 ; 
      while ( (hdata = static_cast<TH1 *>(next())) ) {
        if (hdata && hdata->InheritsFrom("TH1")) { 
          Double_t rv = 0.;
          if(hdata->GetEntries()>0)rv=1; 
          AliDebug(AliQAv1::GetQADebugLevel(), Form("%s -> %f", hdata->GetName(), rv)) ; 
          count[specie]++ ; 
          test[specie] += rv ; 
        }
        else{
          AliError("Data type cannot be processed") ;
        }
      }
      if (count[specie] != 0) { 
        if (test[specie]==0) {
          AliWarning("Histograms are there, but they are all empty: setting flag to kWARNING");
          test[specie] = 0.5;  //upper limit value to set kWARNING flag for a task
        }
        else {
        test[specie] /= count[specie] ;
        }
        AliDebug(AliQAv1::GetQADebugLevel(), Form("Test Result = %f", test[specie])) ; 
      }
    }
  }
}  



