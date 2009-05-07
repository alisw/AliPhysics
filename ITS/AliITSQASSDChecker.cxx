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
#include "TString.h"
#include "Riostream.h"

// --- AliRoot header files ---
#include "AliITSQASSDChecker.h"
#include "AliLog.h"

ClassImp(AliITSQASSDChecker)
//__________________________________________________________________
AliITSQASSDChecker& AliITSQASSDChecker::operator = (const AliITSQASSDChecker& qac ) 
{
  // Equal operator.
  this->~AliITSQASSDChecker();
  new(this) AliITSQASSDChecker(qac);
  return *this;
}

//__________________________________________________________________
Double_t AliITSQASSDChecker::Check(AliQAv1::ALITASK_t /*index*/, TObjArray * list) {  
  AliDebug(1,Form("AliITSQASSDChecker called with offset: %d\n", fSubDetOffset));
  
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
	TString histname = hdata->GetName();
	if(!histname.Contains("fHistSSD")) continue;
        Double_t rv = 0.;
        if(hdata->GetEntries()>0) rv = 1;
        //AliInfo(Form("%s -> %f", hdata->GetName(), rv)) ;
	//cout<<hdata->GetName()<<" - "<<rv<<endl;
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
  
  //AliInfo(Form("Test Result = %f", test)) ;
  //cout<<"Test result: "<<test<<endl;

  return test ;

  //return 0.;
  
}

//__________________________________________________________________
void AliITSQASSDChecker::SetTaskOffset(Int_t TaskOffset)
{
  fSubDetOffset = TaskOffset;
}
