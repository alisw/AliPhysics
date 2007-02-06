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

#include "AliLog.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONConstants.h"
#include "AliMUONData.h"
#include "AliMUONDigit.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONPadStatusMapMaker.h"
#include "AliMUONV2DStore.h"
#include "AliMUONVCalibParam.h"
#include "TClonesArray.h"

/// \class AliMUONDigitCalibrator
/// Class used to calibrate digits (either real or simulated ones).
///
/// The calibration consists of subtracting the pedestal
/// and multiplying by a gain, so that
/// Signal = (ADC-pedestal)*gain
///
/// Please note also that for the moment, if a digit lies on a dead channel
/// we remove this digit from the list of digits.
/// FIXME: this has to be revisited. By using the AliMUONDigit::fFlags we
/// should in principle flag a digit as bad w/o removing it, but this 
/// then requires some changes in the cluster finder to deal with this extra
/// information correctly (e.g. to set a quality for the cluster if it contains
/// bad digits).
///
/// \author Laurent Aphecetche


/// \cond CLASSIMP
ClassImp(AliMUONDigitCalibrator)
/// \endcond

//_____________________________________________________________________________
AliMUONDigitCalibrator::AliMUONDigitCalibrator(AliMUONData* muonData,
                                               AliMUONCalibrationData* calib)
: TTask("AliMUONDigitCalibrator","Raw digit calibration"),
  fData(muonData),
  fCalibrationData(calib),
  fStatusMap(0x0)
{
    /// ctor. This class needs the muonData to get access to the digit,
    /// and the calibrationData to get access to calibration parameters.
    
    if (!calib) throw;
    
    AliMUONPadStatusMaker maker(*calib);
    
    // this is here that we decide on our "goodness" policy, i.e.
    // what do we call an invalid pad (a pad maybe bad because it's HV
    // was too low, or its pedestals too high, etc..)
    //
    maker.SetHVSt12Limits(1300,1600);
    maker.SetHVSt345Limits(1500,2000);
    maker.SetPedMeanLimits(50,200);
    maker.SetPedSigmaLimits(0.1,3);
    
    // From this set of limits, compute the status of all tracker pads.
    AliMUONV2DStore* status = maker.MakeStatus();
    
    AliMUONPadStatusMapMaker mapMaker;
    
    Int_t mask(0x8000000); 
      //FIXME: fake one (consider dead only if ped mean too high or hv switch off)
    
    fStatusMap = mapMaker.MakePadStatusMap(*status,mask);
    
    delete status;
}

//_____________________________________________________________________________
AliMUONDigitCalibrator::~AliMUONDigitCalibrator()
{
  /// dtor.
  delete fStatusMap;
}

//_____________________________________________________________________________
void
AliMUONDigitCalibrator::Exec(Option_t*)
{
  /// Main method.
  /// We loop on tracking chambers (i.e. we do nothing for trigger)
  /// and for each digit in that chamber, we calibrate it :
  /// a) we set its status map and if status is bad, set the signal to zero
  /// b) we then apply pedestal and gain corrections.
  
  for ( Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ++ch )
  {
    TClonesArray* digitArray = fData->Digits(ch);
    Int_t nDigits = digitArray->GetEntriesFast();
    for ( Int_t d = 0; d < nDigits; ++d )
    {
      AliMUONDigit* digit = 
        static_cast<AliMUONDigit*>(digitArray->UncheckedAt(d));
 
      AliMUONVCalibParam* deadmap = static_cast<AliMUONVCalibParam*>
        (fStatusMap->Get(digit->DetElemId(),digit->ManuId()));
      Int_t statusMap = deadmap->ValueAsInt(digit->ManuChannel());
      digit->SetStatusMap(statusMap);
      if ( ( statusMap & AliMUONPadStatusMapMaker::SelfDeadMask() ) != 0 ) // pad itself is bad (not testing its neighbours at this stage)
      {
        digit->SetSignal(0);
        AliWarning(Form("Channel detElemId %d manuId %d "
                        "manuChannel %d is bad %x",digit->DetElemId(),digit->ManuId(),
                        digit->ManuChannel(),digit->StatusMap()));
        continue;
      }
          
      // If the channel is good, go on with the calibration itself.
      
      AliMUONVCalibParam* pedestal = static_cast<AliMUONVCalibParam*>
        (fCalibrationData->Pedestals(digit->DetElemId(),digit->ManuId()));
      
      AliMUONVCalibParam* gain = static_cast<AliMUONVCalibParam*>
        (fCalibrationData->Gains(digit->DetElemId(),digit->ManuId()));
      
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
      Float_t adc = digit->Signal();
      Float_t padc = adc-pedestal->ValueAsFloat(manuChannel,0);
      if ( padc < 3.0*pedestal->ValueAsFloat(manuChannel,1) ) 
      {
        padc = 0.0;
      }
      Float_t charge = padc*gain->ValueAsFloat(manuChannel,0);
      digit->SetSignal(charge);
      Int_t saturation = gain->ValueAsInt(manuChannel,1);
      if ( charge >= saturation )
      {
        digit->Saturated(kTRUE);
      }
    }
  }
}
