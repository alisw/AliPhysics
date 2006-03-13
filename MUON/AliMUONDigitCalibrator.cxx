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

// $Id$

#include "AliMUONDigitCalibrator.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliLog.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONConstants.h"
#include "AliMUONData.h"
#include "AliMUONDigit.h"
#include "AliMUONVCalibParam.h"
#include "AliMpDEManager.h"
#include "AliMpPad.h"
#include "AliMpPlaneType.h"
#include "AliMpStationType.h"
#include "AliMpVSegmentation.h"
#include "Riostream.h"
#include "TClonesArray.h"

ClassImp(AliMUONDigitCalibrator)

//_____________________________________________________________________________
AliMUONDigitCalibrator::AliMUONDigitCalibrator(AliMUONData* muonData,
                                              AliMUONCalibrationData* calib)
: TTask("AliMUONDigitCalibrator","Subtract pedestal from digit charge"),
  fData(muonData),
  fCalibrationData(calib)
{
    //
    // ctor. This class need the muonData to get access to the digit,
    // and the calibrationData to get access to calibration parameters.
    //
}

//______________________________________________________________________________
AliMUONDigitCalibrator::AliMUONDigitCalibrator(const AliMUONDigitCalibrator& right) 
  : TTask(right) 
{  
/// Protected copy constructor (not implemented)

  AliFatal("Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMUONDigitCalibrator::~AliMUONDigitCalibrator()
{
  //
  // empty dtor.
  //
}

//______________________________________________________________________________
AliMUONDigitCalibrator& 
AliMUONDigitCalibrator::operator=(const AliMUONDigitCalibrator& right)
{
/// Protected assignement operator (not implemented)

  // check assignement to self
  if (this == &right) return *this;

  AliFatal("Assignement operator not provided.");
    
  return *this;  
}    

//_____________________________________________________________________________
void
AliMUONDigitCalibrator::Exec(Option_t*)
{
  //
  // Main method.
  // We loop on tracking chambers (i.e. we do nothing for trigger)
  // and for each digit in that chamber, we calibrate it :
  // a) if the corresponding channel is known to be bad, we set the signal to 0
  //    (so that digit can be suppressed later on)
  // b) we then apply pedestal and gain corrections.
  
  for ( Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ++ch )
  {
    TClonesArray* digitArray = fData->Digits(ch);
    Int_t nDigits = digitArray->GetEntriesFast();
    for ( Int_t d = 0; d < nDigits; ++d )
    {
      AliMUONDigit* digit = 
        static_cast<AliMUONDigit*>(digitArray->UncheckedAt(d));
 
      // Very first check is whether this channel is known to be bad,
      // in which case we set the signal to zero.
      AliMUONVCalibParam* dead = static_cast<AliMUONVCalibParam*>
        (fCalibrationData->DeadChannel(digit->DetElemId(),digit->ManuId()));
      if ( dead && dead->ValueAsInt(digit->ManuChannel()) )
      {
        AliDebug(10,Form("Removing dead channel detElemId %d manuId %d "
                        "manuChannel %d",digit->DetElemId(),digit->ManuId(),
                        digit->ManuChannel()));
        digit->SetSignal(0);
        continue;
      }
          
      // If the channel is good, go on with the calibration itself.
      
      AliMUONVCalibParam* pedestal = static_cast<AliMUONVCalibParam*>
        (fCalibrationData->Pedestal(digit->DetElemId(),digit->ManuId()));
      
      AliMUONVCalibParam* gain = static_cast<AliMUONVCalibParam*>
        (fCalibrationData->Gain(digit->DetElemId(),digit->ManuId()));
      
      if (!pedestal)
      {
        AliFatal(Form("Got a null ped object for DE,manu=%d,%d",
                      digit->DetElemId(),digit->ManuId()));
        
      }
      if (!gain)
      {
        AliFatal(Form("Got a null gain object for DE,manu=%d,%d",
                      digit->DetElemId(),digit->ManuId()));        
      }
      
      Int_t manuChannel = digit->ManuChannel();
      Int_t adc = digit->Signal();
      Float_t padc = adc-pedestal->ValueAsFloat(manuChannel,0);
      if ( padc < 3.0*pedestal->ValueAsFloat(manuChannel,1) ) 
      {
        padc = 0.0;
      }
      Float_t charge = padc*gain->ValueAsFloat(manuChannel,0);
      Int_t signal = TMath::Nint(charge);
      digit->SetSignal(signal);
      Int_t saturation = gain->ValueAsInt(manuChannel,1);
      if ( signal >= saturation )
      {
        digit->Saturated(kTRUE);
      }
    }
  }
}
