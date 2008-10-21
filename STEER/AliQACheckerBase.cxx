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

//
//  Base class for detectors quality assurance checkers 
//  Compares Data made by QADataMakers with reference data
//  Y. Schutz CERN August 2007
//

// --- ROOT system ---
#include <TClass.h>
#include <TH1F.h> 
#include <TH1I.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 
#include <TList.h>
#include <TNtupleD.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliQACheckerBase.h"
#include "AliQADataMaker.h"

ClassImp(AliQACheckerBase)

           
//____________________________________________________________________________ 
AliQACheckerBase::AliQACheckerBase(const char * name, const char * title) : 
  TNamed(name, title), 
  fDataSubDir(0x0),
  fRefSubDir(0x0), 
  fRefOCDBSubDir(0x0)
{
  // ctor
}

//____________________________________________________________________________ 
AliQACheckerBase::AliQACheckerBase(const AliQACheckerBase& qac) :
  TNamed(qac.GetName(), qac.GetTitle()),
  fDataSubDir(qac.fDataSubDir), 
  fRefSubDir(qac.fRefSubDir), 
  fRefOCDBSubDir(qac.fRefOCDBSubDir)
{
  //copy ctor
    
}

//____________________________________________________________________________
AliQACheckerBase& AliQACheckerBase::operator = (const AliQACheckerBase& qadm )
{
  // Equal operator.
  this->~AliQACheckerBase();
  new(this) AliQACheckerBase(qadm);
  return *this;
}

//____________________________________________________________________________
Double_t AliQACheckerBase::Check(AliQA::ALITASK_t /*index*/) 
{
  // Performs a basic checking
  // Compares all the histograms stored in the directory
  // With reference histograms either in a file of in OCDB  

	Double_t test = 0.0  ;
	Int_t count = 0 ; 

	if (!fDataSubDir)  
		test = 1. ; // nothing to check
	else 
		if (!fRefSubDir && !fRefOCDBSubDir)
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
					TH1 * href = NULL ; 
					if (fRefSubDir) 
						href  = static_cast<TH1*>(fRefSubDir->Get(key->GetName())) ;
					else if (fRefOCDBSubDir) {
						href  = static_cast<TH1*>(fRefOCDBSubDir->FindObject(key->GetName())) ;
					}
					if (!href) 
						test = -1 ; // no reference data ; 
					else {
						Double_t rv =  DiffK(hdata, href) ;
						AliInfo(Form("%s ->Test = %f", hdata->GetName(), rv)) ; 
						test += rv ; 
						count++ ; 
					}
				} else
					AliError(Form("%s Is a Classname that cannot be processed", key->GetClassName())) ;
			}
		} 
	
	if (count != 0) 
		test /= count ;
   
	return test ;
}  

//____________________________________________________________________________
Double_t AliQACheckerBase::Check(AliQA::ALITASK_t /*index*/, TObjArray * list) 
{
  // Performs a basic checking
  // Compares all the histograms in the list

	Double_t test = 0.0  ;
	Int_t count = 0 ; 

	if (list->GetEntries() == 0)  
		test = 1. ; // nothing to check
	else {
		if (!fRefSubDir)
			test = -1 ; // no reference data
		else {
			TIter next(list) ; 
			TH1 * hdata ;
			count = 0 ; 
			while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
				if ( hdata) { 
					TH1 * href = NULL ; 
					if (fRefSubDir) 
						href  = static_cast<TH1*>(fRefSubDir->Get(hdata->GetName())) ;
					else if (fRefOCDBSubDir)
						href  = static_cast<TH1*>(fRefOCDBSubDir->FindObject(hdata->GetName())) ;
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
					AliError("Data type cannot be processed") ;
			}
		}
	}
	if (count != 0) 
		test /= count ;
	return test ;
}  

//____________________________________________________________________________ 
Double_t AliQACheckerBase::DiffC(const TH1 * href, const TH1 * hin) const
{
  // compares two histograms using the Chi2 test
  if ( hin->Integral() == 0 ) {
    AliWarning(Form("Spectrum %s is empty", hin->GetName())) ; 
    return 0. ;
  }
    
  return hin->Chi2Test(href) ;  
}

//____________________________________________________________________________ 
Double_t AliQACheckerBase::DiffK(const TH1 * href, const TH1 * hin) const
{
  // compares two histograms using the Kolmogorov test
  if ( hin->Integral() == 0 ) {
    AliWarning(Form("Spectrum %s is empty", hin->GetName())) ; 
    return 0. ;
  }
    
  return hin->KolmogorovTest(href) ;  
}

//____________________________________________________________________________ 
void AliQACheckerBase::Init(const AliQA::DETECTORINDEX_t det)
{
  AliQA::Instance(det) ; 
}
 

//____________________________________________________________________________
void AliQACheckerBase::Run(AliQA::ALITASK_t index, TObject * obj) 
{ 
	AliDebug(1, Form("Processing %s", AliQA::GetAliTaskName(index))) ; 
  
	Double_t rv = -1 ;	
  if ( !obj ) {
    rv = Check(index) ;
  } else { 
    TString className(obj->IsA()->GetName()) ; 
    if (className.Contains("TObjArray")) {
      rv = Check(index, static_cast<TObjArray *>(obj)) ;
    } else if (className.Contains("TNtupleD")) { 
      rv = Check(index, static_cast<TNtupleD *>(obj)) ;
    }  else {
      AliError(Form("%s class not implemented", className.Data())) ; 
    }
  }
  
	SetQA(index, rv) ; 	
	
  AliDebug(1, Form("Test result of %s", AliQA::GetAliTaskName(index))) ;
	
  Finish() ; 
}

//____________________________________________________________________________
void AliQACheckerBase::Finish() const 
{
	// wrap up and save QA in proper file
	AliQA * qa = AliQA::Instance() ; 
	qa->Show() ;
	AliQA::GetQAResultFile()->cd() ; 
	qa->Write(qa->GetName(), kWriteDelete) ;   
	AliQA::GetQAResultFile()->Close() ; 
}

//____________________________________________________________________________
void AliQACheckerBase::SetQA(AliQA::ALITASK_t index, const Double_t value) const
{
	// sets the QA according the return value of the Check

	AliQA * qa = AliQA::Instance(index) ; 

	if ( value <= 0.) 
		qa->Set(AliQA::kFATAL) ; 
	else if ( value > 0 && value <= 0.0002 )
		qa->Set(AliQA::kERROR) ; 
	else if ( value > 0.0002 && value <= 0.5 )
		qa->Set(AliQA::kWARNING) ;
	else if ( value > 0.5 && value < 1 ) 
		qa->Set(AliQA::kINFO) ; 		
}
