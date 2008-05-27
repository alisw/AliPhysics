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
#include "AliMUONLogger.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONPadStatusMapMaker.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVStore.h"
#include "AliMpBusPatch.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDetElement.h"

//-----------------------------------------------------------------------------
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
//-----------------------------------------------------------------------------


/// \cond CLASSIMP
ClassImp(AliMUONDigitCalibrator)
/// \endcond

const Int_t AliMUONDigitCalibrator::fgkNoGain(0);
const Int_t AliMUONDigitCalibrator::fgkGainConstantCapa(1);
const Int_t AliMUONDigitCalibrator::fgkGain(2);

//_____________________________________________________________________________
AliMUONDigitCalibrator::AliMUONDigitCalibrator(const AliMUONCalibrationData& calib,
                                               const char* calibMode)
: TObject(),
fLogger(new AliMUONLogger(20000)),
fStatusMaker(0x0),
fStatusMapMaker(0x0),
fPedestals(0x0),
fGains(0x0),
fApplyGains(0),
fCapacitances(0x0)
{
  /// ctor
  
  TString cMode(calibMode);
  cMode.ToUpper();
  
  if ( cMode == "NOGAIN" ) 
  {
    fApplyGains = fgkNoGain;
    AliInfo("Will NOT apply gain correction");
  }
  else if ( cMode == "GAINCONSTANTCAPA" ) 
  {
    fApplyGains = fgkGainConstantCapa;
    AliInfo("Will apply gain correction, but with constant capacitance");
  }
  else if ( cMode == "GAIN" ) 
  {
    fApplyGains = fgkGain;
    AliInfo("Will apply gain correction, with measured capacitances");
  }
  else
  {
    AliError(Form("Invalid calib mode = %s. Will use NOGAIN instead",calibMode));
    fApplyGains = fgkNoGain;
  }
       
  fStatusMaker = new AliMUONPadStatusMaker(calib);
  
  // this is here that we decide on our "goodness" policy, i.e.
  // what do we call an invalid pad (a pad maybe bad because its HV
  // was too low, or its pedestals too high, etc..)
  // FIXME: find a way not to hard-code the goodness policy (i.e. the limits)
  // here...
  fStatusMaker->SetHVSt12Limits(1300,1600);
  fStatusMaker->SetHVSt345Limits(1500,2000);
  fStatusMaker->SetPedMeanLimits(50,200);
  fStatusMaker->SetPedSigmaLimits(0.1,3);
  
  Int_t mask(0x8080); 
  //FIXME: kind of fake one for the moment, we consider dead only 
  // if ped and/or hv value missing.
  //WARNING : getting this mask wrong is a very effective way of getting
  //no digits at all out of this class ;-)
  
  Bool_t deferredInitialization = kTRUE;
  
  fStatusMapMaker = new AliMUONPadStatusMapMaker(*fStatusMaker,mask,deferredInitialization);
  
  fPedestals = calib.Pedestals();

  fGains = calib.Gains(); // we get gains whatever the calibMode is, in order
  // to get the saturation value...

  if ( fApplyGains == fgkGain ) 
  {
    fCapacitances = calib.Capacitances();
  }
}

//_____________________________________________________________________________
AliMUONDigitCalibrator::~AliMUONDigitCalibrator()
{
  /// dtor.
  delete fStatusMaker;
  delete fStatusMapMaker;
  
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
  Int_t detElemId(-1);
  Double_t nsigmas(3.0);
  
  AliDebug(1,Form("# of digits = %d",digitStore.GetSize()));
  
  while ( ( digit = static_cast<AliMUONVDigit*>(next() ) ) )
  {
    if ( digit->DetElemId() != detElemId ) 
    {
      // Find out occupancy of that DE
      detElemId = digit->DetElemId();
      AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
      Double_t nchannels = de->NofChannels();
      Double_t occ = digitStore.GetSize(detElemId)/nchannels;
      if ( occ > 0.05 ) 
      {
        nsigmas = 10.0; // enlarge (a lot) sigma cut if occupancy is high
        // (which probably means zero suppression was not exactly OK).
        fLogger->Log(Form("Will use %5.0f*sigma cut for DE %04d "
                          "due to high occupancy",nsigmas,detElemId));
      }
      else
      {
        nsigmas = 3.0;
      }
    }

     CalibrateDigit(*digit,nsigmas);
  }
}

//_____________________________________________________________________________
void
AliMUONDigitCalibrator::CalibrateDigit(AliMUONVDigit& digit, Double_t nsigmas)
{
  /// Calibrate one digit
  
  if ( digit.IsCalibrated() ) 
  {
    fLogger->Log("ERROR : trying to calibrate a digit twice");
    return;
  }
  
  Int_t statusMap = fStatusMapMaker->StatusMap(digit.DetElemId(),
                                               digit.ManuId(),
                                               digit.ManuChannel());

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
    (fPedestals->FindObject(digit.DetElemId(),digit.ManuId()));

    if (!pedestal)
    {
      // no pedestal -> no charge
      digit.SetCharge(0);
      
      fLogger->Log(Form("Got a null pedestal object for DE,manu=%d,%d",
                        digit.DetElemId(),digit.ManuId()));        
      return;
    }
    
    
    AliMUONVCalibParam* gain = static_cast<AliMUONVCalibParam*>
        (fGains->FindObject(digit.DetElemId(),digit.ManuId()));

    if (!gain)
    {
      if ( fApplyGains != fgkNoGain )
      {
        // no gain -> no charge
        digit.SetCharge(0);

        fLogger->Log(Form("Got a null gain object for DE,manu=%d,%d",
                          digit.DetElemId(),digit.ManuId())); 
        return;
      }
    }

    Int_t manuChannel = digit.ManuChannel();
    Float_t adc = digit.ADC();
    Float_t padc = adc-pedestal->ValueAsFloat(manuChannel,0);
    Float_t charge(0);
    Float_t capa(1.0);
    
    if ( fApplyGains == fgkGainConstantCapa ) 
    {
      capa = 0.2; // pF
    }
    else if ( fApplyGains == fgkGain ) 
    {
      AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(digit.DetElemId());
      
      Int_t serialNumber = de->GetManuSerialFromId(digit.ManuId());

      AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fCapacitances->FindObject(serialNumber));
      
      if ( param )
      {
        capa = param->ValueAsFloat(digit.ManuChannel());
      }
      else
      {
        fLogger->Log(Form("No capa found for serialNumber=%d",serialNumber));
        capa = 0.0;
      }
    }
    
    if ( padc > nsigmas*pedestal->ValueAsFloat(manuChannel,1) ) 
    {
      if ( fApplyGains != fgkNoGain ) 
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
          charge = a0*thres + a0*(padc-thres) + a1*(padc-thres)*(padc-thres);
        }
      }
      else
      {
        charge = padc;
      }
    }
    
    charge *= capa;
    digit.SetCharge(charge);
    
    Int_t saturation(3000);
    
    if ( gain )
    {
      saturation = gain->ValueAsInt(manuChannel,4);
    }
    
    if ( padc >= saturation )
    {
      digit.Saturated(kTRUE);
    }
  }
}
