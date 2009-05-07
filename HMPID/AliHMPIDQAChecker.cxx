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

//...
//  Checks the quality assurance. 
//  By comparing with reference data
//  Skeleton for HMPID
//...

// --- ROOT system ---
#include <TClass.h>
#include <TH1F.h> 
#include <TH1I.h> 
#include <TF1.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliHMPIDQAChecker.h"
#include "AliCDBEntry.h"
#include "AliQAManager.h"

ClassImp(AliHMPIDQAChecker)

//____________________________________________________________________________
Double_t * AliHMPIDQAChecker::Check(AliQAv1::ALITASK_t /*index*/)
{
  Double_t * rv = new Double_t[AliRecoParam::kNSpecies] ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    rv[specie] = 0.0 ; 
  return rv ;  
}

//_________________________________________________________________
Double_t * AliHMPIDQAChecker::Check(AliQAv1::ALITASK_t index, TObjArray ** list) 
{
//
// Main check function: Depending on the TASK, different checks are applied
// At the moment:       check for empty histograms and checks for RecPoints

  Double_t * check = new Double_t[AliRecoParam::kNSpecies] ; 
  
  AliInfo(Form("Fix needed ....."));
  char * detOCDBDir = Form("HMPID/%s/%s", AliQAv1::GetRefOCDBDirName(), AliQAv1::GetRefDataDirName()) ; 
  AliCDBEntry *QARefRec = AliQAManager::QAManager()->Get(detOCDBDir);
  if( !QARefRec){
    AliInfo("QA reference data NOT retrieved for Reconstruction check. No HMPIDChecker  ...exiting");
    return check;
  }

// checking for empy histograms
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    check[specie] = 1.0;
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue ; 
    if(CheckEntries(list[specie]) == 0)  {
      AliWarning("histograms are empty");
      check[specie] = 0.4;//-> Corresponds to kWARNING see AliQACheckerBase::Run
    }
  
    // checking rec points
    if(index == AliQAv1::kREC) check[specie] = CheckRecPoints(list[specie],(TObjArray *)QARefRec->GetObject());

    //default check response. It will be changed when reasonable checks will be considered
    else check[specie] = 0.7 ; // /-> Corresponds to kINFO see AliQACheckerBase::Run 
  }
  return check;

}
//_________________________________________________________________
Double_t AliHMPIDQAChecker::CheckEntries(TObjArray * list) const
{
  //
  //  check on the QA histograms on the input list: 
  // 

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
        //Printf("hitogram %s     has entries: %f ",hdata->GetName(),hdata->GetEntries());
	if(hdata->GetEntries()>0)rv=1; 
	count++ ; 
	test += rv ; 
      }
      else{
	AliError("Data type cannot be processed") ;
      }
      
    }
    if (count != 0) { 
      if (test==0) {
	AliWarning("Histograms are booked for THIS specific Task, but they are all empty: setting flag to kWARNING");
	test = 0.;  //upper limit value to set kWARNING flag for a task
      }
      else test = 1 ;
    }
  }

  return test ; 
}  
//_________________________________________________________________

Double_t AliHMPIDQAChecker::CheckRecPoints(TObjArray *listrec, TObjArray *listref) const
{
  //
  //  check on the HMPID RecPoints by using expo fit and Kolmogorov Test:
  //

   Float_t checkresponse = 0;

   Float_t counter = 0 ;
   TIter next(listrec) ;
   TH1* histo;
   while ( (histo = dynamic_cast<TH1 *>(next())) ) {
   if( histo->GetEntries() < 3 ) counter++;
   else {
    TString h = histo->GetTitle();
    if(h.Contains("Zoom")){
    histo->Fit("expo","Q0","",1,50);
    if(histo->GetFunction("expo")->GetParameter(1) !=0 ) if(TMath::Abs((-1./(histo->GetFunction("expo"))->GetParameter(1)) - 35 ) > 5) counter++;
    }
    if(h.Contains("size  MIP"))   if(TMath::Abs(histo->GetMean()-5) > 2) counter++;
    if(h.Contains("size  Phots")) if(TMath::Abs(histo->GetMean()-2) > 2) counter++;
    if(h.Contains("distribution")) if(histo->KolmogorovTest((TH1F *)listref->At(0))<0.8) counter++;
    AliDebug(1,Form(" Kolm. test : %f ",histo->KolmogorovTest((TH1F *)listref->At(0))));  
   }
  }
 Float_t response = counter/(7.+7.+42.+42.); // 7.+7.+42 +42 = N checked histograms (-> To be replaced by listrec->GetEntries())
 
 if(response < 0.1) checkresponse = 0.7;      // <10% of the check histograms show a failing check -> Corresponds to kINFO see AliQACheckerBase::Run
 else if(response < 0.5) checkresponse = 0.4; //  50%  of the check histograms show a failing check -> Corresponds to kWARNING see AliQACheckerBase::Run
 else checkresponse = 0.001;                  // > 50% of the check histograms show a failing check -> Corresponds to kERROR see AliQACheckerBase::Run
 return checkresponse;
}
//________________________________________________________________

