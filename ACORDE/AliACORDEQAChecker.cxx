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
//...
//  Checks the quality assurance for ACORDE. 
//  Default implementation
//  Authors:
//	Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch> (FCFM-BUAP) 
//	Luciano Diaz Gonzalez <luciano.diaz@nucleares.unam.mx> (ICN-UNAM)
//	Arturo Fernandez <afernan@mail.cern.ch> (FCFM-BUAP)
//...

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
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliACORDEQAChecker.h"

ClassImp(AliACORDEQAChecker)

//__________________________________________________________________

Double_t AliACORDEQAChecker::Check(AliQA::ALITASK_t /*index*/, TObjArray * list)
{


// Super-basic check on the QA histograms on the input list: 
  // look whether they are empty!
  Double_t test = 0.0  ;
  Int_t count = 0 ; 
  
  if (list->GetEntries() == 0){  
    test = 1. ; // nothing to check
  }
  else {
    TIter next(list) ; 
    TH1 * hdata ;
    count = 0 ; 
    while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
      if (hdata) { 
	Double_t rv = 0.;
	if(hdata->GetEntries()>0)rv=1; 
	AliInfo(Form("%s -> %f", hdata->GetName(), rv)) ; 
	count++ ; 
	test += rv ; 
      }
      else{
	AliError("Data type cannot be processed") ;
      }
      
    }
    if (count != 0) { 
      if (test==0) {
	AliWarning("Histograms are there, but they are all empty: setting flag to kWARNING");
	test = 0.5;  //upper limit value to set kWARNING flag for a task
      }
      else {
	test /= count ;
      }
    }
  }

  AliInfo(Form("Test Result = %f", test)) ; 
  return test ; 



}

