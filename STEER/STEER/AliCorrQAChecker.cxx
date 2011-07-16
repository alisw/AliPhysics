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


/* $Id: $ */

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
#include <TNtupleD.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliCorrQAChecker.h"

ClassImp(AliCorrQAChecker)

//__________________________________________________________________
Double_t * AliCorrQAChecker::CheckN(AliQAv1::ALITASK_t index, TNtupleD ** nData, AliDetectorRecoParam * /*recoParam*/) 
{
 // check the QA of correlated data stored in a ntuple
  
  Double_t * test = new Double_t[AliRecoParam::kNSpecies] ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    test[specie] = 0. ; 
    
  if ( index != AliQAv1::kRAW ) {
    AliWarning("Checker not implemented") ; 
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
      test[specie] = 1. ; 
    return test ; 
  }
//	if (!fRefSubDir) {
//		test = 1. ; // no reference data
//	} else {
    if ( ! nData ) {
      AliError(Form("nRawCorr not found in %s", fDataSubDir->GetName())) ; 
    } else {
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
        if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
          continue ; 
        TObjArray * bList = nData[specie]->GetListOfBranches() ; 
        for (Int_t b = 0 ; b < bList->GetEntries() ; b++) {
          AliInfo(Form("Ntuple %s parameter name %d : %s", nData[specie]->GetName(), b, bList->At(b)->GetName())) ;  
        }
      }
    }
 // }
  return test ; 
}
//__________________________________________________________________
void   AliCorrQAChecker::Run(AliQAv1::ALITASK_t tsk, TNtupleD ** nt, AliDetectorRecoParam * recoParam) 
{
    // special run for TNtupleD
	AliDebug(AliQAv1::GetQADebugLevel(), Form("Processing %s", AliQAv1::GetAliTaskName(tsk))) ; 
  
	Double_t * rv = NULL ;
  rv = CheckN(tsk, nt, recoParam) ;
	SetQA(tsk, rv) ; 	
	
  AliDebug(AliQAv1::GetQADebugLevel(), Form("Test result of %s", AliQAv1::GetAliTaskName(tsk))) ;
	
  if (rv) 
    delete [] rv ; 
  Finish() ; 
} 

