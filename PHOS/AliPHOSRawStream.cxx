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

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to PHOS digits in raw data.
///
/// It loops over all PHOS digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE.
/// Several getters provide information about the current digit.
/// usage: 
/// root > AliRawReaderFile rawReader ; 
/// root > AliPHOSRawStream input(&rawReader) ; 
/// root > while (input.Next()) ..... 
///////////////////////////////////////////////////////////////////////////////

#include "TClonesArray.h"
#include "TClass.h"



#include "AliPHOSRawStream.h"
#include "AliRawReader.h"
#include "AliPHOSConTableDB.h"
#include "AliPHOSDigit.h"

#ifdef ALI_DATE
#include "event.h"
#else
#include "AliPHOSevent.h"
#endif

ClassImp(AliPHOSRawStream)

//_____________________________________________________________________________
AliPHOSRawStream::AliPHOSRawStream(AliRawReader* rawReader) : TObject()
{
  fRawReader = rawReader ;
  fctdb = 0 ;
}
//_____________________________________________________________________________
Bool_t AliPHOSRawStream::ReadDigits(TClonesArray * digits){

  Bool_t isOK = kFALSE ;
  if(!fctdb){
    Error("ReadDigits","Connection table not set") ;
    return kFALSE ;
  }

  if(!digits){
    Error("ReadDigits","Output array not created") ;
    return kFALSE ;
  }

  if(!(digits->GetClass()->InheritsFrom("AliPHOSDigit"))){
    Error("ReadDigits","Digits contanier made for %s, not AliPHOSDigits",digits->GetClass()->GetName()) ;
    return kFALSE ;
  }

  digits->Clear() ;

  //Check, if current event - PHYSICS event
  if(!((fRawReader->GetType() & EVENT_TYPE_MASK)==PHYSICS_EVENT)){
    return kFALSE ;
  }

  //Scan subevents until subevent with digits
  while(fRawReader->ReadNextData(fData)){    
    switch (fRawReader->GetEquipmentType()){
    case kPattUnitMarker:
      if(fRawReader->GetEquipmentId() == kPattUnitEquipId){ //Read PHOS trigger
	Int_t * patt = (Int_t *)fData;
	if(fRawReader->GetEquipmentSize() >= (Int_t)sizeof(Int_t))
	  fTrig = patt[0];
	else
	  fTrig = 0 ;
      }
      break;
    case kPhosAdcMarker:
      if(fRawReader->GetEquipmentId() == kPhosAdcEquipId){
	Int_t ndigits = fRawReader->GetEquipmentSize()/sizeof(Int_t);      
	digits->Expand(ndigits) ;
	for(Int_t i=0; i<ndigits; i++){
	  Int_t * amps = (Int_t *)fData ;
	  Int_t absID = fctdb->Raw2AbsId(i) ;
	  Int_t time = 0;
	  if(absID>0) //Only real digits are filled, last ADC numbers (if any) are scipped
	    new((*digits)[i]) AliPHOSDigit( -1, absID, amps[i], time) ;
	}
	digits->Sort() ;
	digits->Compress() ;  
	digits->Expand(digits->GetEntriesFast()) ;

	//Set indexes in list of digits and make true digitization of the energy
	for (Int_t id = 0 ; id < digits->GetEntriesFast() ; id++) { 
	  AliPHOSDigit * digit = dynamic_cast<AliPHOSDigit*>( digits->At(id) ) ; 
	  digit->SetIndexInList(id) ;     
	}

	isOK = kTRUE ;
      }
      break;
    case kTdcMarker:
      if(fRawReader->GetEquipmentId() == kTdcEquipId){
	// Probably sometime we will need to handle these data 
      }
      break;
    case kChargeAdcMarker:
      if(fRawReader->GetEquipmentId() == kChargeAdcEquipId){
	//Probably sometime we will need to handle these data 
      }
      break;
    case kScalerMarker:
      if(fRawReader->GetEquipmentId() == kScalerEquipId){
	//Probably sometime we will need to handle these data 
      }
      break;
    default:
      break;
    }
  }
  return isOK ;
}

