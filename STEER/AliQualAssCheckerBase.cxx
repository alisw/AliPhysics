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
  fData(0x0), 
  fDetectorDir(0x0),
  fRef(0x0) 
{
  // ctor
  Init() ; 
}

//____________________________________________________________________________ 
AliQualAssCheckerBase::AliQualAssCheckerBase(const AliQualAssCheckerBase& qadm) :
  TNamed(qadm.GetName(), qadm.GetTitle()),
  fData(qadm.fData),
  fDetectorDir(qadm.fDetectorDir), 
  fRef(qadm.fRef)
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
const Double_t AliQualAssCheckerBase::DiffC(const TH1 * href, const TH1 * hin) const
{
  // compares two histograms using the Chi2 test
  if ( hin->Integral() == 0 ) {
    AliWarning(Form("Spectrum %s in %s is empty", hin->GetName(), fData->GetName())) ; 
    return 0. ;
  }
    
  return hin->Chi2Test(href) ;  
}

//____________________________________________________________________________ 
const Double_t AliQualAssCheckerBase::DiffK(const TH1 * href, const TH1 * hin) const
{
  // compares two histograms using the Kolmogorov test
  if ( hin->Integral() == 0 ) {
    AliWarning(Form("Spectrum %s in %s is empty", hin->GetName(), fData->GetName())) ; 
    return 0. ;
  }
    
  return hin->KolmogorovTest(href) ;  
}

//____________________________________________________________________________ 
void AliQualAssCheckerBase::Init()
{
  //open files and search for the appropriate detector directory
  
  fRef = AliQualAssChecker::GetRefFile() ;
  if (!fRef)
    AliFatal(Form("Reference file %s not found", AliQualAssChecker::GetRefFileName())) ; 
	
  fDetectorDir = fRef->GetDirectory(AliQualAssDataMaker::GetDetectorDirName()) ; 
  if (!fDetectorDir)
    AliFatal(Form("Directory %s not found in reference file %s not found", AliQualAssDataMaker::GetDetectorDirName(), AliQualAssChecker::GetRefFileName())) ; 

  fData = AliQualAssChecker::GetDataFile() ;
  if (!fData)
    AliFatal(Form("Reference file %s not found", AliQualAss::GetOutputName())) ; 

  if (! fData->GetDirectory(AliQualAssDataMaker::GetDetectorDirName())) ; 
    AliFatal(Form("Directory %s not found in reference file %s not found", AliQualAssDataMaker::GetDetectorDirName(), AliQualAss::GetOutputName())) ; 
  
  AliQualAss::Instance(AliQualAss::kSIM) ; 
}
 
//____________________________________________________________________________
void AliQualAssCheckerBase::Exec(const Option_t * opt) 
{ 
  AliInfo(Form("Processing %s", opt)) ; 
// loop over detectors
  AliQualAss * qa = AliQualAss::Instance(AliQualAss::kPHOS) ; 
  Double_t rv = Check(opt) ;   
 if ( rv <= 0.) 
    qa->Set(AliQualAss::kFATAL) ; 
  else if ( rv > 0 && rv <= 0.2 )
    qa->Set(AliQualAss::kERROR) ; 
  else if ( rv > 0.2 && rv <= 0.5 )
    qa->Set(AliQualAss::kWARNING) ;
  else if ( rv > 0.5 && rv < 1 ) 
    qa->Set(AliQualAss::kINFO) ; 
  AliInfo(Form("Test result of %s in PHOS", opt)) ;
  Finish() ; 
// endloop     
}

//____________________________________________________________________________
void AliQualAssCheckerBase::Finish() const 
{
  // wrap up and save QA in proper file
    
  AliQualAss * qa = AliQualAss::Instance() ; 
  qa->Show() ;
  AliQualAssChecker::GetOutFile()->cd() ; 
  qa->Write(qa->GetName(), kWriteDelete) ;   
  AliQualAssChecker::GetOutFile()->Close() ;  
}
