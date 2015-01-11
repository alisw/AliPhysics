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
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliADQAChecker.h"
#include "AliADQADataMakerRec.h"

ClassImp(AliADQAChecker)

//__________________________________________________________________
AliADQAChecker::AliADQAChecker() : AliQACheckerBase("AD","AD Quality Assurance Data Checker"),
  fLowEventCut(1000),
  fORvsANDCut(0.2),
  fBGvsBBCut(0.2)
{
  // Default constructor
  // Nothing else here
}

//__________________________________________________________________
void AliADQAChecker::Check(Double_t * check, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * /*recoParam*/) 
{
  // Main check function: Depending on the TASK, different checks will be applied
  // Check for missing channels and check on the trigger type for raw data
  // Check for missing disk or rings for esd (to be redone)

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    check[specie] = 1.0;
    // no check on cosmic or calibration events
    if (AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCosmic || AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCalib)
      continue;
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue;
    if (index == AliQAv1::kRAW) {
      check[specie] =  CheckRaws(list[specie]);
    } else if (index == AliQAv1::kESD) {
      // Check for one disk missing (FATAL) or one ring missing (ERROR) in ESDs (to be redone)
      check[specie] =  CheckEsds(list[specie]);
    } 
  }
}

//_________________________________________________________________
Double_t AliADQAChecker::CheckRaws(TObjArray * list) const
{

  //  Check on the QA histograms on the raw-data input list:
  //  Two things are checked: the presence of data in all channels and
  //  the ratio between different trigger types

  Double_t test = 1.0;
  if (list->GetEntries() == 0){  
    AliWarning("There are no histograms to be checked");
  } 
  return test ; 
}  

//_________________________________________________________________
Double_t AliADQAChecker::CheckEsds(TObjArray * list) const
{
  
//  check the ESDs for missing disk or ring
//  printf(" Number of entries in ESD list = %d\n", list->GetEntries()); 
//  list->Print();

  Double_t test     = 1.0;     // initialisation to OK

  return test ; 
} 

//______________________________________________________________________________
void AliADQAChecker::Init(const AliQAv1::DETECTORINDEX_t det) 
{
  // intialises QA and QA checker settings
  AliQAv1::Instance(det) ; 
  Float_t * hiValue = new Float_t[AliQAv1::kNBIT] ; 
  Float_t * lowValue = new Float_t[AliQAv1::kNBIT] ;
  lowValue[AliQAv1::kINFO]      = 0.5   ; 
  hiValue[AliQAv1::kINFO]       = 1.0 ; 
  lowValue[AliQAv1::kWARNING]   = 0.2 ; 
  hiValue[AliQAv1::kWARNING]    = 0.5 ; 
  lowValue[AliQAv1::kERROR]     = 0.0   ; 
  hiValue[AliQAv1::kERROR]      = 0.2 ; 
  lowValue[AliQAv1::kFATAL]     = -1.0   ; 
  hiValue[AliQAv1::kFATAL]      = 0.0 ; 
  SetHiLo(hiValue, lowValue) ; 
  delete [] hiValue;
  delete [] lowValue;
}

//______________________________________________________________________________
void AliADQAChecker::SetQA(AliQAv1::ALITASK_t index, Double_t * value) const
{
// sets the QA word according to return value of the Check
  AliQAv1 * qa = AliQAv1::Instance(index);
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    qa->UnSet(AliQAv1::kFATAL, specie);
    qa->UnSet(AliQAv1::kWARNING, specie);
    qa->UnSet(AliQAv1::kERROR, specie);
    qa->UnSet(AliQAv1::kINFO, specie);
    if ( ! value ) { // No checker is implemented, set all QA to Fatal
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
  
