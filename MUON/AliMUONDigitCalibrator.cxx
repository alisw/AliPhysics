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
#include "AliMUONCalibParam.h"
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
}

//_____________________________________________________________________________
AliMUONDigitCalibrator::~AliMUONDigitCalibrator()
{
}

//_____________________________________________________________________________
void
AliMUONDigitCalibrator::Exec(Option_t*)
{
  for ( Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ++ch )
  {
    TClonesArray* digitArray = fData->Digits(ch);
    Int_t nDigits = digitArray->GetEntriesFast();
    for ( Int_t d = 0; d < nDigits; ++d )
    {
      AliMUONDigit* digit = 
        static_cast<AliMUONDigit*>(digitArray->UncheckedAt(d));
 
      AliMUONCalibParam* pedestal = static_cast<AliMUONCalibParam*>
        (fCalibrationData->Pedestal(digit->DetElemId(),
                                    digit->ManuId(),digit->ManuChannel()));
      
      AliMUONCalibParam* gain = static_cast<AliMUONCalibParam*>
        (fCalibrationData->Gain(digit->DetElemId(),
                                    digit->ManuId(),digit->ManuChannel()));
      if (!pedestal)
      {
        AliFatal(Form("Got a null ped object for DE,manu,channel=%d,%d,%d",
                      digit->DetElemId(),digit->ManuId(),digit->ManuChannel()));
        
      }
      if (!gain)
      {
        AliFatal(Form("Got a null gain object for DE,manu,channel=%d,%d,%d",
                      digit->DetElemId(),digit->ManuId(),digit->ManuChannel()));
        
      }
      
      Int_t adc = digit->Signal();
      Float_t padc = adc-pedestal->Mean();
      if ( padc < 3.0*pedestal->Sigma() ) 
      {
        padc = 0.0;
      }
      Float_t charge = padc*gain->Mean();
      Int_t signal = TMath::Nint(charge);
      digit->SetSignal(signal);
    }
  }
}
