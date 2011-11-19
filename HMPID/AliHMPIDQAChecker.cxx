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
 //_________________________________________________________________
AliHMPIDQAChecker::AliHMPIDQAChecker() : 
AliQACheckerBase("HMPID","HMPID Quality Assurance Data Checker"), 
fNoReference(kTRUE),
fQARefRec(NULL)
{
    //ctor, fetches the reference data from OCDB 
  char * detOCDBDir = Form("HMPID/%s/%s", AliQAv1::GetRefOCDBDirName(), AliQAv1::GetRefDataDirName()) ; 
  AliCDBEntry * QARefRec = AliQAManager::QAManager()->Get(detOCDBDir);
  if(QARefRec) {
    fQARefRec = dynamic_cast<TObjArray*> (QARefRec->GetObject()) ; 
    if (fQARefRec)
      if (fQARefRec->GetEntries()) 
        fNoReference = kFALSE ;            
    if (fNoReference) 
      AliInfo("QA reference data NOT retrieved for Reconstruction check. No HMPIDChecker!");
  }
}

//_________________________________________________________________
AliHMPIDQAChecker::AliHMPIDQAChecker(const AliHMPIDQAChecker& qac) : 
AliQACheckerBase(qac.GetName(), qac.GetTitle()), 
fNoReference(qac.fNoReference), 
fQARefRec(NULL)
{
  fNoReference = qac.fNoReference ; 
  if (qac.fQARefRec) {
    fQARefRec = new TObjArray(qac.fQARefRec->GetEntries()) ; 
    for (Int_t index=0; index < qac.fQARefRec->GetEntries(); index++) 
      fQARefRec->Add(qac.fQARefRec->At(index)) ; 
  }
}

//_________________________________________________________________
AliHMPIDQAChecker::~AliHMPIDQAChecker() 
{
  if(fQARefRec) { fQARefRec->Delete() ;   delete fQARefRec ; }
}
//_________________________________________________________________
void AliHMPIDQAChecker::Check(Double_t *  check, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * /*recoParam*/) 
{
//
// Main check function: Depending on the TASK, different checks are applied
// At the moment:       check for empty histograms and checks for RecPoints

  if(fNoReference)  

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    check[specie] = 1.0;
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue ; 
    // checking for empy histograms
    if(CheckEntries(list[specie]) == 0)  {
      AliWarning("histograms are empty");
      check[specie] = 0.4;//-> Corresponds to kWARNING see AliQACheckerBase::Run
    }
  
    //check sim
    if(index == AliQAv1::kSIM) check[specie] = CheckSim(list[specie], fQARefRec);

    // checking rec points
    if(index == AliQAv1::kREC) check[specie] = CheckRec(list[specie], fQARefRec);
   
    //default check response. It will be changed when reasonable checks will be considered
    else check[specie] = 0.7 ; // /-> Corresponds to kINFO see AliQACheckerBase::Run 
  } // species loop
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
Double_t AliHMPIDQAChecker::CheckSim(TObjArray *listsim, TObjArray *listref) const
{
  //
  //  check on the HMPID RecPoints by using expo fit and Kolmogorov Test:
  //

   Float_t checkresponse = 0;

   Float_t counter = 0 ;
   TIter next(listsim) ;
   TH1* histo;
   while ( (histo = dynamic_cast<TH1 *>(next())) ) {
     //PH The histogram should have at least 10 bins with at least 5 entries
     Int_t nbinsabove = 0;
     for (Int_t ibin=histo->FindBin(1); ibin<=histo->FindBin(50); ibin++) { 
       if (histo->GetBinContent(ibin)>5) nbinsabove++;
     }

   if( nbinsabove < 10 ) counter++;
   else {
    TString h = histo->GetTitle();
    if(h.Contains("Zoom")){
    histo->Fit("expo","LQ0","",5,50);
    if(histo->GetFunction("expo")->GetParameter(1) !=0 ) if(TMath::Abs((-1./(histo->GetFunction("expo"))->GetParameter(1)) - 35 ) > 5) counter++;
    }
    if(h.Contains("size  MIP"))   if(TMath::Abs(histo->GetMean()-5) > 2) counter++;
    if(h.Contains("size  Phots")) if(TMath::Abs(histo->GetMean()-2) > 2) counter++;
    if(h.Contains("distribution")) if(histo->KolmogorovTest((TH1F *)listref->At(0))<0.8) counter++;
    AliDebug(AliQAv1::GetQADebugLevel(),Form(" Kolm. test : %f ",histo->KolmogorovTest((TH1F *)listref->At(0))));  
   }
  }
 Float_t response = counter/(7.+7.+42.+42.); // 7.+7.+42 +42 = N checked histograms (-> To be replaced by listsim->GetEntries())
 
 if(response < 0.1) checkresponse = 0.7;      // <10% of the check histograms show a failing check -> Corresponds to kINFO see AliQACheckerBase::Run
 else if(response < 0.5) checkresponse = 0.4; //  50%  of the check histograms show a failing check -> Corresponds to kWARNING see AliQACheckerBase::Run
 else checkresponse = 0.001;                  // > 50% of the check histograms show a failing check -> Corresponds to kERROR see AliQACheckerBase::Run
 return checkresponse;
}

//___________________________________________________________________________________________________
Double_t AliHMPIDQAChecker::CheckRec(TObjArray *listrec, TObjArray *listref) const
{
  //
  //  check on the HMPID RecPoints by using expo fit and Kolmogorov Test:
  //

   Float_t checkresponse = 0;

   Float_t counter = 0 ;
   TIter next(listrec) ;
   TH1* histo;
   while ( (histo = dynamic_cast<TH1 *>(next())) ) {
     //PH The histogram should have at least 10 bins with at least 5 entries
     Int_t nbinsabove = 0;
     for (Int_t ibin=histo->FindBin(1); ibin<=histo->FindBin(50); ibin++) { 
       if (histo->GetBinContent(ibin)>5) nbinsabove++;
     }

   if( nbinsabove < 10 ) counter++;
   else {
    TString h = histo->GetTitle();
    if(h.Contains("Zoom")){
    histo->Fit("expo","LQ0","",5,50);
    if(histo->GetFunction("expo")->GetParameter(1) !=0 ) if(TMath::Abs((-1./(histo->GetFunction("expo"))->GetParameter(1)) - 35 ) > 5) counter++;
    }
    if(h.Contains("size  MIP"))   if(TMath::Abs(histo->GetMean()-5) > 2) counter++;
    if(h.Contains("size  Phots")) if(TMath::Abs(histo->GetMean()-2) > 2) counter++;
    if(h.Contains("distribution")) if(histo->KolmogorovTest((TH1F *)listref->At(0))<0.8) counter++;
    AliDebug(AliQAv1::GetQADebugLevel(),Form(" Kolm. test : %f ",histo->KolmogorovTest((TH1F *)listref->At(0))));  
   }
  }
 Float_t response = counter/(7.+7.+42.+42.); // 7.+7.+42 +42 = N checked histograms (-> To be replaced by listrec->GetEntries())
 
 if(response < 0.1) checkresponse = 0.7;      // <10% of the check histograms show a failing check -> Corresponds to kINFO see AliQACheckerBase::Run
 else if(response < 0.5) checkresponse = 0.4; //  50%  of the check histograms show a failing check -> Corresponds to kWARNING see AliQACheckerBase::Run
 else checkresponse = 0.001;                  // > 50% of the check histograms show a failing check -> Corresponds to kERROR see AliQACheckerBase::Run
 return checkresponse;
}

