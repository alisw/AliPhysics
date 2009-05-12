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
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliQACheckerBase.h"
#include "AliQADataMaker.h"

ClassImp(AliQACheckerBase)

           
//____________________________________________________________________________ 
AliQACheckerBase::AliQACheckerBase(const char * name, const char * title) : 
  TNamed(name, title), 
  fDataSubDir(0x0),
  fRefSubDir(0x0), 
  fRefOCDBSubDir(0x0), 
  fLowTestValue(0x0),
  fUpTestValue(0x0)
{
  // ctor
  fLowTestValue = new Float_t[AliQAv1::kNBIT] ; 
  fUpTestValue  = new Float_t[AliQAv1::kNBIT] ; 
  fLowTestValue[AliQAv1::kINFO]    =  0.5   ; 
  fUpTestValue[AliQAv1::kINFO]     = 1.0 ; 
  fLowTestValue[AliQAv1::kWARNING] =  0.002 ; 
  fUpTestValue[AliQAv1::kWARNING]  = 0.5 ; 
  fLowTestValue[AliQAv1::kERROR]   =  0.0   ; 
  fUpTestValue[AliQAv1::kERROR]    = 0.002 ; 
  fLowTestValue[AliQAv1::kFATAL]   = -1.0   ; 
  fUpTestValue[AliQAv1::kFATAL]    = 0.0 ; 
  
  AliDebug(AliQAv1::GetQADebugLevel(), "Default setting is:") ;
  if ( AliDebugLevel()  == AliQAv1::GetQADebugLevel() ) {
    printf( "                      INFO    -> %1.5f <  value <  %1.5f \n", fLowTestValue[AliQAv1::kINFO], fUpTestValue[AliQAv1::kINFO]) ; 
    printf( "                      WARNING -> %1.5f <  value <= %1.5f \n", fLowTestValue[AliQAv1::kWARNING], fUpTestValue[AliQAv1::kWARNING]) ; 
    printf( "                      ERROR   -> %1.5f <  value <= %1.5f \n", fLowTestValue[AliQAv1::kERROR], fUpTestValue[AliQAv1::kERROR]) ; 
    printf( "                      FATAL   -> %1.5f <= value <  %1.5f \n", fLowTestValue[AliQAv1::kFATAL], fUpTestValue[AliQAv1::kFATAL]) ; 
  }
}

//____________________________________________________________________________ 
AliQACheckerBase::AliQACheckerBase(const AliQACheckerBase& qac) :
  TNamed(qac.GetName(), qac.GetTitle()),
  fDataSubDir(qac.fDataSubDir), 
  fRefSubDir(qac.fRefSubDir), 
  fRefOCDBSubDir(qac.fRefOCDBSubDir), 
  fLowTestValue(qac.fLowTestValue),
  fUpTestValue(qac.fLowTestValue)
{
  //copy ctor
  for (Int_t index = 0 ; index < AliQAv1::kNBIT ; index++) {
    fLowTestValue[index]  = qac.fLowTestValue[index] ; 
    fUpTestValue[index] = qac.fUpTestValue[index] ; 
  }
    
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
AliQACheckerBase::~AliQACheckerBase()
{
  delete [] fLowTestValue ; 
  delete [] fUpTestValue ; 
}

//____________________________________________________________________________
Double_t * AliQACheckerBase::Check(AliQAv1::ALITASK_t /*index*/) 
{
  // Performs a basic checking
  // Compares all the histograms stored in the directory
  // With reference histograms either in a file of in OCDB  

	Double_t * test = new Double_t[AliRecoParam::kNSpecies] ;
	Int_t count[AliRecoParam::kNSpecies]   = { 0 }; 

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    test[specie] = 1.0 ; 
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue ; 
    if (!fDataSubDir) {
      test[specie] = 0. ; // nothing to check
    } else if (!fRefSubDir && !fRefOCDBSubDir) {
        test[specie] = -1 ; // no reference data
    } else {
      TList * keyList = fDataSubDir->GetListOfKeys() ; 
      TIter next(keyList) ; 
      TKey * key ;
      count[specie] = 0 ; 
      while ( (key = static_cast<TKey *>(next())) ) {
        TObject * odata = fRefSubDir->Get(key->GetName()) ; 
        if ( odata->IsA()->InheritsFrom("TH1") ) {
          TH1 * hdata = static_cast<TH1*>(odata) ;
          TH1 * href = NULL ; 
          if (fRefSubDir) 
            href  = static_cast<TH1*>(fRefSubDir->Get(key->GetName())) ;
          else if (fRefOCDBSubDir[specie]) {  
            href  = static_cast<TH1*>(fRefOCDBSubDir[specie]->FindObject(key->GetName())) ;
          }
          if (!href) 
            test[specie] = -1 ; // no reference data ; 
          else {
            Double_t rv =  DiffK(hdata, href) ;
            AliDebug(AliQAv1::GetQADebugLevel(), Form("%s ->Test = %f", hdata->GetName(), rv)) ; 
            test[specie] += rv ; 
            count[specie]++ ; 
          }
        } else
          AliError(Form("%s Is a Classname that cannot be processed", key->GetClassName())) ;
      }
      if (count[specie] != 0) 
        test[specie] /= count[specie] ;
    }
  }
 	return test ;
}  

//____________________________________________________________________________
Double_t * AliQACheckerBase::Check(AliQAv1::ALITASK_t /*index*/, TObjArray ** list) 
{
  // Performs a basic checking
  // Compares all the histograms in the list

	Double_t * test = new Double_t[AliRecoParam::kNSpecies] ;
	Int_t count[AliRecoParam::kNSpecies]   = { 0 }; 

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    test[specie] = 1.0 ; 
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue ; 
    if (list[specie]->GetEntries() == 0)  
      test[specie] = 0. ; // nothing to check
    else {
      if (!fRefSubDir && !fRefOCDBSubDir)
        test[specie] = -1 ; // no reference data
      else {
        TIter next(list[specie]) ; 
        TH1 * hdata ;
        count[specie] = 0 ; 
        while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
          if ( hdata) { 
            if ( hdata->TestBit(AliQAv1::GetExpertBit()) )  // does not perform the test for expert data
              continue ; 
            TH1 * href = NULL ; 
            if (fRefSubDir) 
              href  = static_cast<TH1*>(fRefSubDir->Get(hdata->GetName())) ;
            else if (fRefOCDBSubDir[specie])
              href  = static_cast<TH1*>(fRefOCDBSubDir[specie]->FindObject(hdata->GetName())) ;
            if (!href) 
              test[specie] = -1 ; // no reference data ; 
            else {
              Double_t rv =  DiffK(hdata, href) ;
              AliDebug(AliQAv1::GetQADebugLevel(), Form("%s ->Test = %f", hdata->GetName(), rv)) ; 
              test[specie] += rv ; 
              count[specie]++ ; 
            }
          } 
          else
            AliError("Data type cannot be processed") ;
        }
        if (count[specie] != 0) 
          test[specie] /= count[specie] ;
      }
    }
  }
	return test ;
}  


//____________________________________________________________________________ 
Double_t AliQACheckerBase::DiffC(const TH1 * href, const TH1 * hin) const
{
  // compares two histograms using the Chi2 test
  if ( hin->Integral() == 0 ) {
    AliDebug(AliQAv1::GetQADebugLevel(), Form("Spectrum %s is empty", hin->GetName())) ; 
    return 0. ;
  }
    
  return hin->Chi2Test(href) ;  
}

//____________________________________________________________________________ 
Double_t AliQACheckerBase::DiffK(const TH1 * href, const TH1 * hin) const
{
  // compares two histograms using the Kolmogorov test
  if ( hin->Integral() == 0 ) {
    AliDebug(AliQAv1::GetQADebugLevel(), Form("Spectrum %s is empty", hin->GetName())) ; 
    return 0. ;
  }
    
  return hin->KolmogorovTest(href) ;  
}

//____________________________________________________________________________
void AliQACheckerBase::Run(AliQAv1::ALITASK_t index, TObjArray ** list) 
{ 
	AliDebug(AliQAv1::GetQADebugLevel(), Form("Processing %s", AliQAv1::GetAliTaskName(index))) ; 
  
	Double_t * rv = NULL ;
  if ( !list) 
    rv = Check(index) ;
  else 
    rv = Check(index, list) ;
	SetQA(index, rv) ; 	
	
  AliDebug(AliQAv1::GetQADebugLevel(), Form("Test result of %s", AliQAv1::GetAliTaskName(index))) ;
	
  if (rv) 
    delete [] rv ; 
  Finish() ; 
}

//____________________________________________________________________________
void AliQACheckerBase::Finish() const 
{
	// wrap up and save QA in proper file
	AliQAv1 * qa = AliQAv1::Instance() ; 
	qa->Show() ;
	AliQAv1::GetQAResultFile()->cd() ; 
	qa->Write(qa->GetName(), kWriteDelete) ;   
	AliQAv1::GetQAResultFile()->Close() ; 
}

//____________________________________________________________________________
void AliQACheckerBase::SetHiLo(Float_t * hiValue, Float_t * lowValue) 
{
  AliDebug(AliQAv1::GetQADebugLevel(), "Previous setting was:") ;
  if ( AliDebugLevel() == AliQAv1::GetQADebugLevel() ) {
    printf( "                      INFO    -> %1.5f <  value <  %1.5f \n", fLowTestValue[AliQAv1::kINFO], fUpTestValue[AliQAv1::kINFO]) ; 
    printf( "                      WARNING -> %1.5f <  value <= %1.5f \n", fLowTestValue[AliQAv1::kWARNING], fUpTestValue[AliQAv1::kWARNING]) ; 
    printf( "                      ERROR   -> %1.5f <  value <= %1.5f \n", fLowTestValue[AliQAv1::kERROR], fUpTestValue[AliQAv1::kERROR]) ; 
    printf( "                      FATAL   -> %1.5f <= value <  %1.5f \n", fLowTestValue[AliQAv1::kFATAL], fUpTestValue[AliQAv1::kFATAL]) ; 
  }
  
  for (Int_t index = 0 ; index < AliQAv1::kNBIT ; index++) {
    fLowTestValue[index]  = lowValue[index] ; 
    fUpTestValue[index] = hiValue[index] ; 
  }
  AliDebug(AliQAv1::GetQADebugLevel(), "Current setting is:") ;
  if ( AliDebugLevel()  == AliQAv1::GetQADebugLevel() ) {
    printf( "                      INFO    -> %1.5f <  value <  %1.5f \n", fLowTestValue[AliQAv1::kINFO], fUpTestValue[AliQAv1::kINFO]) ; 
    printf( "                      WARNING -> %1.5f <  value <= %1.5f \n", fLowTestValue[AliQAv1::kWARNING], fUpTestValue[AliQAv1::kWARNING]) ; 
    printf( "                      ERROR   -> %1.5f <  value <= %1.5f \n", fLowTestValue[AliQAv1::kERROR], fUpTestValue[AliQAv1::kERROR]) ; 
    printf( "                      FATAL   -> %1.5f <= value <  %1.5f \n", fLowTestValue[AliQAv1::kFATAL], fUpTestValue[AliQAv1::kFATAL]) ; 
  }
}

//____________________________________________________________________________
void AliQACheckerBase::SetQA(AliQAv1::ALITASK_t index, Double_t * value) const
{
	// sets the QA according the return value of the Check

  AliQAv1 * qa = AliQAv1::Instance(index) ;
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (  value == NULL ) { // No checker is implemented, set all QA to Fatal
      qa->Set(AliQAv1::kFATAL, specie) ; 
    } else {
      if ( value[specie] >= fLowTestValue[AliQAv1::kFATAL] && value[specie] < fUpTestValue[AliQAv1::kFATAL] ) 
        qa->Set(AliQAv1::kFATAL, specie) ; 
      else if ( value[specie] > fLowTestValue[AliQAv1::kERROR] && value[specie] <= fUpTestValue[AliQAv1::kERROR]  )
        qa->Set(AliQAv1::kERROR, specie) ; 
      else if ( value[specie] > fLowTestValue[AliQAv1::kWARNING] && value[specie] <= fUpTestValue[AliQAv1::kWARNING]  )
        qa->Set(AliQAv1::kWARNING, specie) ;
      else if ( value[specie] > fLowTestValue[AliQAv1::kINFO] && value[specie] <= fUpTestValue[AliQAv1::kINFO] ) 
        qa->Set(AliQAv1::kINFO, specie) ; 	
    }
  }
}
