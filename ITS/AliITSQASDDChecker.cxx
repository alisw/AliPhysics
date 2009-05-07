
/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

// *****************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  P. Cerello Apr 2008
//  INFN Torino

// --- ROOT system ---
#include "TH1.h"
#include <TCanvas.h>

// --- AliRoot header files ---
#include "AliITSQASDDChecker.h"
#include "AliLog.h"
#include "AliCDBEntry.h"
#include "AliQAManager.h"
#include "AliQACheckerBase.h"
#include "TSystem.h"


ClassImp(AliITSQASDDChecker)
//__________________________________________________________________
AliITSQASDDChecker& AliITSQASDDChecker::operator = (const AliITSQASDDChecker& qac ) 
{
  // Equal operator.
  this->~AliITSQASDDChecker();
  new(this) AliITSQASDDChecker(qac);
  return *this;
}

//__________________________________________________________________
Double_t AliITSQASDDChecker::Check(AliQAv1::ALITASK_t index, TObjArray * list) 
{  
  AliDebug(1,Form("AliITSQASDDChecker called with offset: %d\n", fSubDetOffset));
  char * detOCDBDir = Form("ITS/%s/%s", AliQAv1::GetRefOCDBDirName(), AliQAv1::GetRefDataDirName()) ; 
  AliCDBEntry *QARefObj = AliQAManager::QAManager()->Get(detOCDBDir);
  if( !QARefObj){
    AliError("Calibration object retrieval failed! SDD will not be processed");
    return 1.;
  }

  Double_t test = 0.0;
  Int_t offset = 0;

  if(index==AliQAv1::kRAW){  //analizing RAWS
    TH1F *ModPattern = (TH1F*)QARefObj->GetObject();
    if (list->GetEntries() == 0){
      test = 1. ; // nothing to check
    }
    else {
      TIter next(list) ;
      TH1 * hdata ;
      for(offset =0;offset < fSubDetOffset; offset++){
        hdata = dynamic_cast<TH1 *>(next());  
      }

      while ( (hdata = dynamic_cast<TH1 *>(next())) && offset >= fSubDetOffset){
	if (hdata) {
	  if(offset == fSubDetOffset){ //ModPattern check
	    if ( hdata->Integral() == 0 ) {
	      AliWarning(Form("Spectrum %s is empty", hdata->GetName())) ; 
	      return 0.5 ;
	    }
	    test = hdata->Chi2Test(ModPattern,"UU,p");
	  }  // ModPattern check
	}
	else{
	  AliError("Data type cannot be processed") ;
	}
	offset++;
      }  //SDD histo

      while ( (hdata = dynamic_cast<TH1 *>(next()))) {
	offset++;
      }
    } //else entries !=0
  AliInfo(Form("Test Result for RAWS = %f", test)) ;
  } // if(index==0)

  
  if( index==AliQAv1::kREC){ //analizing RECP
    //printf("analizing recp, offset %d \n",fSubDetOffset);
    if (list->GetEntries() == 0){
      test = 1. ; // nothing to check
    }
    else {
      TIter next(list) ;
      TH1 * hdata ;
      for(offset =0;offset < fSubDetOffset; offset++){
	hdata = dynamic_cast<TH1 *>(next());    // magari TIter++ ?? 
	//printf("Skipping histo %s, offset %d \n",hdata->GetName(),fSubDetOffset);
      }
        
      while ( (hdata = dynamic_cast<TH1 *>(next())) && offset >= fSubDetOffset ){
	if (hdata) { // offset=9 ModPatternRP
	  //printf("Treating histo %s, offset %d \n",hdata->GetName(),fSubDetOffset);
	  if( offset == 9 && hdata->GetEntries()>0)test = 0.1;	  
	}
	else{
	  AliError("Data type cannot be processed") ;
	}
	offset++;
      }
    } // GetEntries loop
  AliInfo(Form("Test Result for RECP = %f", test)) ; 
  } // if(index==2) loop

  return test;	
}
 
//__________________________________________________________________
void AliITSQASDDChecker::SetTaskOffset(Int_t TaskOffset)
{
  fSubDetOffset = TaskOffset;
}
