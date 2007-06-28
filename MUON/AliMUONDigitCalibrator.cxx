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
#include "AliMpConstants.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONLogger.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONPadStatusMapMaker.h"
#include "AliMUONVStore.h"
#include "AliMUONVCalibParam.h"

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
AliMUONDigitCalibrator::AliMUONDigitCalibrator(const AliMUONCalibrationData& calib,
                                               Bool_t createAndUseStatusMap)
: TObject(),
fCalibrationData(calib),
fStatusMap(0x0),
fLogger(new AliMUONLogger(1000))
{
  /// ctor
  if (createAndUseStatusMap) 
  {
    AliMUONPadStatusMaker maker(fCalibrationData);
    
    // this is here that we decide on our "goodness" policy, i.e.
    // what do we call an invalid pad (a pad maybe bad because its HV
    // was too low, or its pedestals too high, etc..)
    // FIXME: find a way not to hard-code the goodness policy (i.e. the limits)
    // here...
    maker.SetHVSt12Limits(1300,1600);
    maker.SetHVSt345Limits(1500,2000);
    maker.SetPedMeanLimits(50,200);
    maker.SetPedSigmaLimits(0.1,3);
    
    // From this set of limits, compute the status of all tracker pads.
    AliMUONVStore* status = maker.MakeStatus();      
    // we do not check that status is != 0x0, as this is supposed to be
    // the responsability of the padStatusMaker.
    
    AliMUONPadStatusMapMaker mapMaker(fCalibrationData);
    
    Int_t mask(0x8080); 
    //FIXME: kind of fake one for the moment, we consider dead only 
    // if ped and/or hv value missing.
    //WARNING : getting this mask wrong is a very effective way of getting
    //no digits at all out of this class ;-)
    
    fStatusMap = mapMaker.MakePadStatusMap(*status,mask);
    
    delete status;
    }
  else
  {
    // make a fake (empty) status map
    fStatusMap = AliMUONPadStatusMapMaker::MakeEmptyPadStatusMap();
  }
}

//_____________________________________________________________________________
AliMUONDigitCalibrator::~AliMUONDigitCalibrator()
{
  /// dtor.
  delete fStatusMap;
  
  AliInfo("Summary of messages:");
  fLogger->Print();

  delete fLogger;
}

//_____________________________________________________________________________
void
AliMUONDigitCalibrator::Calibrate(AliMUONVDigitStore& digitStore)
{
  /// Calibrate the digits contained in digitStore  
  TIter next(digitStore.CreateTrackerIterator());
  AliMUONVDigit* digit;
  
  while ( ( digit = static_cast<AliMUONVDigit*>(next() ) ) )
  {
    CalibrateDigit(*digit);
  }
}

//_____________________________________________________________________________
void
AliMUONDigitCalibrator::CalibrateDigit(AliMUONVDigit& digit)
{
  /// Calibrate one digit
  
  if ( digit.IsCalibrated() ) 
  {
    fLogger->Log("ERROR : trying to calibrate a digit twice");
    return;
  }
  
  AliMUONVCalibParam* deadmap = static_cast<AliMUONVCalibParam*>
  (fStatusMap->FindObject(digit.DetElemId(),digit.ManuId()));
  Int_t statusMap = deadmap->ValueAsInt(digit.ManuChannel());
  digit.SetStatusMap(statusMap);
  digit.Calibrated(kTRUE);
  if ( ( statusMap & AliMUONPadStatusMapMaker::SelfDeadMask() ) != 0 ) 
  {
    // pad itself is bad (not testing its neighbours at this stage)
    digit.SetCharge(0);
    fLogger->Log(Form("%s:%d:Channel detElemId %d manuId %d "
                    "manuChannel %d is bad %x",__FILE__,__LINE__,
                    digit.DetElemId(),digit.ManuId(),
                    digit.ManuChannel(),digit.StatusMap()));
  }
  else
  {
    // If the channel is good, go on with the calibration itself.

    AliMUONVCalibParam* pedestal = static_cast<AliMUONVCalibParam*>
    (fCalibrationData.Pedestals(digit.DetElemId(),digit.ManuId()));
    
    AliMUONVCalibParam* gain = static_cast<AliMUONVCalibParam*>
      (fCalibrationData.Gains(digit.DetElemId(),digit.ManuId()));
    
    if (!pedestal)
    {
      AliFatal(Form("Got a null ped object for DE,manu=%d,%d",
                    digit.DetElemId(),digit.ManuId()));
      
    }
    if (!gain)
    {
      AliFatal(Form("Got a null gain object for DE,manu=%d,%d",
                    digit.DetElemId(),digit.ManuId()));        
    }
    
    Int_t manuChannel = digit.ManuChannel();
    Float_t adc = digit.ADC();
    Float_t padc = adc-pedestal->ValueAsFloat(manuChannel,0);
    Float_t charge(0);
    if ( padc > 3.0*pedestal->ValueAsFloat(manuChannel,1) ) 
    {
      Float_t a0 = gain->ValueAsFloat(manuChannel,0);
      Float_t a1 = gain->ValueAsFloat(manuChannel,1);
      Int_t thres = gain->ValueAsInt(manuChannel,2);
      if ( padc < thres ) 
      {
        charge = a0*padc;
      }
      else
      {
        charge = a0*thres + a0*padc + a1*padc*padc;
      }
    }
    digit.SetCharge(charge);
    Int_t saturation = gain->ValueAsInt(manuChannel,4);
    if ( charge >= saturation )
    {
      digit.Saturated(kTRUE);
    }
  }
}
