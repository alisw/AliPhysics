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


/*
  Checks the quality assurance. Under construction. 
  By comparing with reference data

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
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliVZEROQAChecker.h"
//#include "AliCDBEntry.h"
//#include "AliCDBManager.h"

ClassImp(AliVZEROQAChecker)

//__________________________________________________________________
const Double_t AliVZEROQAChecker::Check(AliQA::ALITASK_t index, TObjArray * list) 
{

// Main check function: Depending on the TASK, different checks will be applied
// Check for empty histograms 

//   AliDebug(1,Form("AliVZEROChecker"));
//   AliCDBEntry *QARefRec = AliCDBManager::Instance()->Get("VZERO/QARef/RAW");
//   if( !QARefRec){
//     AliInfo("QA reference data NOT retrieved for QA check...");
//     return 1.;
//   }
     
  if ( index == AliQA::kRAW ) 
  {
       return CheckEntries(list);
  }
  
  AliWarning(Form("Checker for task %d not implemented for the moment",index));
  return 0.0;    

}
//_________________________________________________________________
Double_t AliVZEROQAChecker::CheckEntries(TObjArray * list) const
{
  
  //  check on the QA histograms on the input list: 
   

  Double_t test = 0.0  ;
  Int_t   count = 0 ; 

  if (list->GetEntries() == 0){  
      test = 1.0; 
      AliInfo(Form("There are no entries to be checked..."));
  }
  else {
      TIter next(list) ; 
      TH1 * hdata ;
      count = 0 ; 
      while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
        if (hdata) { 
	  Double_t rv = 0.0;
//          Printf("Histogram %s     has entries: %f ",hdata->GetName(),hdata->GetEntries());
	  if(hdata->GetEntries()>0)rv=1.0; 
	  count++ ; 
	  test += rv ; 
        }
        else{
	  AliError(Form("Data type cannot be processed"));
        }      
      }
      if (count != 0) { 
        if (test==0.0) {
	    AliInfo(Form("Histograms are booked for THIS specific Task, but they are all empty"));
        }
        else test = 1.0;
      }
  }

  return test ; 
}  
//______________________________________________________________________________
void AliVZEROQAChecker::SetQA(AliQA::ALITASK_t index, const Double_t value) const
{
// sets the QA word according to return value of the Check

  AliQA * qa = AliQA::Instance(index);
  
  qa->UnSet(AliQA::kFATAL);
  qa->UnSet(AliQA::kWARNING);
  qa->UnSet(AliQA::kERROR);
  qa->UnSet(AliQA::kINFO);
  
  if (value == 1.0)      {qa->Set(AliQA::kINFO);}
  else if (value == 0.0) {qa->Set(AliQA::kFATAL);}
  else if (value > 0.5)  {qa->Set(AliQA::kWARNING);}
  else                   {qa->Set(AliQA::kERROR);}
  
}
