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
  Base class for detectors quality assurance checkers 
  Compares Data made by QualAssDataMakers with reference data
  Y. Schutz CERN August 2007
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
#include "AliQualAss.h"
#include "AliQualAssChecker.h"
#include "AliQualAssCheckerBase.h"
#include "AliQualAssDataMaker.h"

ClassImp(AliQualAssCheckerBase)

           
//____________________________________________________________________________ 
AliQualAssCheckerBase::AliQualAssCheckerBase(const char * name, const char * title) : 
  TNamed(name, title), 
  fDataSubDir(0x0),
  fRefSubDir(0x0) 
{
  // ctor
}

//____________________________________________________________________________ 
AliQualAssCheckerBase::AliQualAssCheckerBase(const AliQualAssCheckerBase& qadm) :
  TNamed(qadm.GetName(), qadm.GetTitle()),
  fDataSubDir(qadm.fDataSubDir), 
  fRefSubDir(qadm.fRefSubDir)
{
  //copy ctor
    
}

//____________________________________________________________________________
AliQualAssCheckerBase& AliQualAssCheckerBase::operator = (const AliQualAssCheckerBase& qadm )
{
  // Equal operator.
  this->~AliQualAssCheckerBase();
  new(this) AliQualAssCheckerBase(qadm);
  return *this;
}

//____________________________________________________________________________
const Double_t AliQualAssCheckerBase::Check() 
{
  // Performs a basic checking
  // Compares all the histograms stored in the directory

   Double_t test = 0.0  ;
   Int_t count = 0 ; 

   if (!fDataSubDir)  
     test = 1. ; // nothing to check
   else 
     if (!fRefSubDir)
       test = -1 ; // no reference data
     else {
       TList * keyList = fDataSubDir->GetListOfKeys() ; 
       TIter next(keyList) ; 
       TKey * key ;
       count = 0 ; 
       while ( (key = static_cast<TKey *>(next())) ) {
	 TObject * odata = fRefSubDir->Get(key->GetName()) ; 
	 if ( odata->IsA()->InheritsFrom("TH1") ) {
	   TH1 * hdata = static_cast<TH1*>(odata) ; 
	   TH1 * href  = static_cast<TH1*>(fRefSubDir->Get(key->GetName())) ;
	   if (!href) 
	     test = -1 ; // no reference data ; 
	   else {
	     Double_t rv =  DiffK(hdata, href) ;
	     AliInfo(Form("%s ->Test = %f", hdata->GetName(), rv)) ; 
	     test += rv ; 
	     count++ ; 
	   }
	 }
	 else
	   AliError(Form("%s Is a Classname that cannot be processed", key->GetClassName())) ;
       }

     }
   if (count != 0) 
     test /= count ;
   
   return test ; 
}  

//____________________________________________________________________________ 
const Double_t AliQualAssCheckerBase::DiffC(const TH1 * href, const TH1 * hin) const
{
  // compares two histograms using the Chi2 test
  if ( hin->Integral() == 0 ) {
    AliWarning(Form("Spectrum %s in %s is empty", hin->GetName(), AliQualAss::GetDataName())) ; 
    return 0. ;
  }
    
  return hin->Chi2Test(href) ;  
}

//____________________________________________________________________________ 
const Double_t AliQualAssCheckerBase::DiffK(const TH1 * href, const TH1 * hin) const
{
  // compares two histograms using the Kolmogorov test
  if ( hin->Integral() == 0 ) {
    AliWarning(Form("Spectrum %s in %s is empty", hin->GetName(), AliQualAss::GetDataName())) ; 
    return 0. ;
  }
    
  return hin->KolmogorovTest(href) ;  
}

//____________________________________________________________________________ 
void AliQualAssCheckerBase::Init(const AliQualAss::DETECTORINDEX det)
{
  AliQualAss::Instance(det) ; 
}
 
//____________________________________________________________________________
void AliQualAssCheckerBase::Run(AliQualAss::ALITASK index) 
{ 
  AliInfo(Form("Processing %s", AliQualAss::GetAliTaskName(index))) ; 

  AliQualAss * qa = AliQualAss::Instance(index) ; 

  Double_t rv = Check() ;   
  if ( rv <= 0.) 
    qa->Set(AliQualAss::kFATAL) ; 
  else if ( rv > 0 && rv <= 0.2 )
    qa->Set(AliQualAss::kERROR) ; 
  else if ( rv > 0.2 && rv <= 0.5 )
    qa->Set(AliQualAss::kWARNING) ;
  else if ( rv > 0.5 && rv < 1 ) 
    qa->Set(AliQualAss::kINFO) ; 
  AliInfo(Form("Test result of %s", AliQualAss::GetAliTaskName(index))) ;
  Finish() ; 
}

//____________________________________________________________________________
void AliQualAssCheckerBase::Finish() const 
{
  // wrap up and save QA in proper file
    
  AliQualAss * qa = AliQualAss::Instance() ; 
  qa->Show() ;
  AliQualAssChecker::GetQAResultFile()->cd() ; 
  qa->Write(qa->GetName(), kWriteDelete) ;   
  AliQualAssChecker::GetQAResultFile()->Close() ;  
}
