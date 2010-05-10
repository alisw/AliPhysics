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
#include "AliVZEROQAChecker.h"
#include "AliVZEROQADataMakerRec.h"
//#include "AliCDBEntry.h"
//#include "AliCDBManager.h"

ClassImp(AliVZEROQAChecker)

//__________________________________________________________________
void AliVZEROQAChecker::Check(Double_t * check, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * /*recoParam*/) 
{

// Main check function: Depending on the TASK, different checks will be applied
// Check for empty histograms 

//   AliDebug(AliQAv1::GetQADebugLevel(),Form("AliVZEROChecker"));
//   AliCDBEntry *QARefRec = AliCDBManager::Instance()->Get("VZERO/QARef/RAW");
//   if( !QARefRec){
//     AliDebug(AliQAv1::GetQADebugLevel(), "QA reference data NOT retrieved for QA check...");
//     return 1.;
//   }

//   Check that histos are filled - (FATAL) set if empty
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    check[specie]    = 1.0 ; 
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue ; 
    if(CheckEntries(list[specie]) == 0.0){
        check[specie] =  CheckEntries(list[specie]);
    } else {
      //   Check for one disk missing (FATAL) or one ring missing (ERROR) in ESDs     
      if(index == AliQAv1::kESD) 
          check[specie] =  CheckEsds(list[specie]);
    } 
  }
}

//_________________________________________________________________
Double_t AliVZEROQAChecker::CheckEntries(TObjArray * list) const
{

  //  check on the QA histograms on the input list: list
//  list->Print();

  Double_t test = 0.0 ;
  Int_t   count = 0 ; 

  if (list->GetEntries() == 0){  
	test = 1.0; 
    AliDebug(AliQAv1::GetQADebugLevel(), Form("There are NO ENTRIES to be checked..."));
  } else {
	TIter next(list) ; 
	TH1 * hdata ;
	count = 0 ; 
	while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
		if (hdata) { 	   
			Double_t rv = 0.0;
			if(hdata->GetEntries()>0) rv=1.0;
//	   AliDebug(AliQAv1::GetQADebugLevel(), Form("%s -> %f", hdata->GetName(), rv)); 
			count++ ;        // number of histos
			test += rv ;     // number of histos filled
        }else{
			AliError(Form("Data type cannot be processed"));
        }      
	}
	if (count != 0) { 
		if (test==0.0) {
			AliWarning(Form("Histograms are BOOKED for this specific task, but they are all EMPTY"));
		} else {
			test /= count; 
		}
	}
  }
  return test ; 
}  

//_________________________________________________________________
Double_t AliVZEROQAChecker::CheckEsds(TObjArray * list) const
{
  
//  check the ESDs for missing disk or ring

//  printf(" Number of entries in ESD list = %d\n", list->GetEntries()); 
//  list->Print();

  Double_t test     = 1.0;     // initialisation to OK
  Int_t    histonb =   0; 
  Double_t multV0A = 0.0;
  Double_t multV0C = 0.0;
  Double_t v0ABBRing[4], v0CBBRing[4];
  Double_t v0ABGRing[4], v0CBGRing[4];
  for (Int_t i=0; i<4; i++) { 
       v0ABBRing[i]= v0CBBRing[i]= 0.0;
       v0ABGRing[i]= v0CBGRing[i]= 0.0;
  }  
  TIter next(list) ; 
  TH1 * hdata ;
  
  while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
	  if (hdata) {
		  switch (histonb) {
		  case AliVZEROQADataMakerRec::kCellMultiV0A:
			  multV0A  = hdata->GetMean();
			  break;
		  case AliVZEROQADataMakerRec::kCellMultiV0C:
			  multV0C  = hdata->GetMean();
			  break;
		  case AliVZEROQADataMakerRec::kBBFlag:
	          for (Int_t i=0; i<8; i++) {         
				  if(i<4) v0CBBRing[i]  = hdata->Integral((i*8)+1, (i*8) +8);
				  else v0ABBRing[i-4]  = hdata->Integral((i*8)+1, (i*8) +8);
			  }	      
			  break;
		  case AliVZEROQADataMakerRec::kBGFlag:
	          for (Int_t i=0; i<8; i++) {         
				  if(i<4) v0CBGRing[i]  = hdata->Integral((i*8)+1, (i*8) +8);
				  else v0ABGRing[i-4]  = hdata->Integral((i*8)+1, (i*8) +8);
			  }	      
			  break;
		  }
	  }
	  histonb++;
  }
  
  if(multV0A == 0.0 || multV0C == 0.0) {
     AliWarning(Form("One of the two disks is missing !") );
     test = 0.0; // bit FATAL set
  }
  if( v0ABBRing[0]+v0ABGRing[0] == 0.0 || 
      v0ABBRing[1]+v0ABGRing[1] == 0.0 || 
      v0ABBRing[2]+v0ABGRing[2] == 0.0 || 
      v0ABBRing[3]+v0ABGRing[3] == 0.0 || 
      v0CBBRing[0]+v0CBGRing[0] == 0.0 || 
      v0CBBRing[1]+v0CBGRing[1] == 0.0 || 
      v0CBBRing[2]+v0CBGRing[2] == 0.0 || 
      v0CBBRing[3]+v0CBGRing[3] == 0.0  ){    
      AliWarning(Form("One ring is missing !") );
      test = 0.1;   // bit ERROR set
  }

  return test ; 
} 

//______________________________________________________________________________
void AliVZEROQAChecker::Init(const AliQAv1::DETECTORINDEX_t det) 
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
}

//______________________________________________________________________________
void AliVZEROQAChecker::SetQA(AliQAv1::ALITASK_t index, Double_t * value) const
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
  
