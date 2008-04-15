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
//  Checks the quality assurance. 
//  By comparing with reference data
//  Skeleton for T0
//...

// --- ROOT system ---
#include <TClass.h>
#include <TH1F.h> 
#include <TH1I.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 
#include <TMath.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliT0QAChecker.h"

ClassImp(AliT0QAChecker)

//__________________________________________________________________
AliT0QAChecker& AliT0QAChecker::operator = (const AliT0QAChecker& qac )
{
  // Equal operator.
  this->~AliT0QAChecker();
  new(this) AliT0QAChecker(qac);
  return *this;
}
//__________________________________________________________________

const Double_t AliT0QAChecker::Check(TObjArray * list)
{

  // Super-basic check on the QA histograms on the input list:
  // look whether they are empty!

  Double_t test = 0.0  ;
  Int_t count = 0 ;
  Double_t nent[100];
  memset(nent,0,100*sizeof(Double_t));
  Double_t w[100];
  memset(w,1.,100*sizeof(Double_t));


  if (list->GetEntries() == 0){
    test = 1. ; // nothing to check
  }
  else {
    
    TIter next(list) ;
    TH1 * hdata ;
    count = 0 ;
    while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
      if (hdata) {
	nent[count] = hdata->GetEntries();
       	AliDebug(1,Form("count %i %s -> %f",count, hdata->GetName(),nent[count])) ;

        Double_t rv = 0.;
        if(hdata->GetEntries()>0) rv = 1;
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
	AliDebug(10,Form(" MaxElement %f ", TMath::MaxElement(count,nent)));	
	if(TMath::MaxElement(count,nent) > 1000) {
	Double_t mean = TMath::Mean(count,nent,w);
	AliDebug(10,Form(" Mean %f ", mean));	
	for (Int_t i=0; i<count; i++) 
	  {
	    Double_t diff = TMath::Abs(nent[i]-mean);
	    if (diff > 0.1*mean )
	      AliInfo(Form("Problem in Number of entried in hist %i  is %f\n", i, nent[i])) ; 
	  }
	}
      }
    }
  }
  AliInfo(Form("Test Result = %f", test)) ;
  return test ;
}
