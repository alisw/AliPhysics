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
#include "AliMUONRecoParam.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVStore.h"
#include "AliMpBusPatch.h"
#include "AliMpConstants.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDetElement.h"
#include "AliMpManuStore.h"

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
                                               const AliMUONRecoParam* recoParams,
                                               const char* calibMode)
: TObject(),
fLogger(new AliMUONLogger(20000)),
fStatusMaker(0x0),
fStatusMapMaker(0x0),
fPedestals(0x0),
fGains(0x0),
fApplyGains(0),
fCapacitances(0x0),
fNumberOfBadPads(0),
fNumberOfPads(0),
fChargeSigmaCut(0)
{
  /// ctor
  
  Ctor(calibMode,calib,recoParams);
}

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
fCapacitances(0x0),
fNumberOfBadPads(0),
fNumberOfPads(0),
fChargeSigmaCut(0)
{
  /// ctor
  
  Ctor(calibMode,calib,0x0);
}

//_____________________________________________________________________________
void
AliMUONDigitCalibrator::Ctor(const char* calibMode,
                             const AliMUONCalibrationData& calib,
                             const AliMUONRecoParam* recoParams)
{
  /// designated ctor
  
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
       
  // Load mapping manu store
  if ( ! AliMpCDB::LoadManuStore() ) {
    AliFatal("Could not access manu store from OCDB !");
  }

  fStatusMaker = new AliMUONPadStatusMaker(calib);
  
  // Set default values, as loose as reasonable

  fChargeSigmaCut = 3.0;
  
	Int_t mask(0x8080); // reject pads where ped *or* hv are missing
	
	if ( recoParams )
	{
    // if we have reco params, we use limits and cuts from there :
    
		fStatusMaker->SetHVSt12Limits(recoParams->HVSt12LowLimit(),recoParams->HVSt12HighLimit());
		fStatusMaker->SetHVSt345Limits(recoParams->HVSt345LowLimit(),recoParams->HVSt345HighLimit());
		fStatusMaker->SetPedMeanLimits(recoParams->PedMeanLowLimit(),recoParams->PedMeanHighLimit());
		fStatusMaker->SetPedSigmaLimits(recoParams->PedSigmaLowLimit(),recoParams->PedSigmaHighLimit());
		fStatusMaker->SetGainA1Limits(recoParams->GainA1LowLimit(),recoParams->GainA1HighLimit());
		fStatusMaker->SetGainA2Limits(recoParams->GainA2LowLimit(),recoParams->GainA2HighLimit());
		fStatusMaker->SetGainThresLimits(recoParams->GainThresLowLimit(),recoParams->GainThresHighLimit());
		
    mask = recoParams->PadGoodnessMask();
		//WARNING : getting this mask wrong is a very effective way of getting
		//no digits at all out of this class ;-)
    
    fChargeSigmaCut = recoParams->ChargeSigmaCut();
	}
  
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

	AliInfo(Form("We have seen %g pads, and rejected %g (%7.2f %%)",
							 fNumberOfPads,fNumberOfBadPads,
							 ( fNumberOfPads > 0 ) ? fNumberOfBadPads*100.0/fNumberOfPads : 0 ));
	
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
  Double_t nsigmas = fChargeSigmaCut;
  
  AliDebug(1,Form("# of digits = %d",digitStore.GetSize()));
  
  while ( ( digit = static_cast<AliMUONVDigit*>(next() ) ) )
  {
    if ( digit->IsCalibrated() ) 
    {
      fLogger->Log("ERROR : trying to calibrate a digit twice");
      return;
    }
    
    digit->Calibrated(kTRUE);
    
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
        nsigmas = fChargeSigmaCut;
      }
    }

    Float_t charge(0.0);
    Int_t statusMap;
    Bool_t isSaturated(kFALSE);
    
    Bool_t ok = IsValidDigit(digit->DetElemId(),digit->ManuId(),digit->ManuChannel(),&statusMap);

    digit->SetStatusMap(statusMap);
    
    if (ok)
    {
      ++fNumberOfPads;
      charge = CalibrateDigit(digit->DetElemId(),digit->ManuId(),digit->ManuChannel(),
                              digit->ADC(),nsigmas,&isSaturated);
    }
    else
    {
      ++fNumberOfBadPads;
    }
    
    digit->SetCharge(charge);
    digit->Saturated(isSaturated);
  }
}

//_____________________________________________________________________________
Float_t
AliMUONDigitCalibrator::CalibrateDigit(Int_t detElemId, Int_t manuId, Int_t manuChannel,
                                       Float_t adc, Float_t nsigmas, 
                                       Bool_t* isSaturated) const

{
  /// Calibrate one digit
  
  
  AliMUONVCalibParam* pedestal = static_cast<AliMUONVCalibParam*>
  (fPedestals->FindObject(detElemId,manuId));
  
  if (!pedestal)
  {
    // no pedestal -> no charge    
    fLogger->Log(Form("Got a null pedestal object for DE,manu=%d,%d",detElemId,manuId));        
    return 0.0;
  }
  
  
  AliMUONVCalibParam* gain = static_cast<AliMUONVCalibParam*>
  (fGains->FindObject(detElemId,manuId));
  
  if (!gain)
  {
    if ( fApplyGains != fgkNoGain )
    {
      // no gain -> no charge
      fLogger->Log(Form("Got a null gain object for DE,manu=%d,%d",
                        detElemId,manuId)); 
      return 0.0;
    }
  }
  
  Float_t padc = adc-pedestal->ValueAsFloat(manuChannel,0);
  Float_t charge(0);
  Float_t capa(1.0);
  
  if ( fApplyGains == fgkGainConstantCapa ) 
  {
    capa = 0.2; // pF
  }
  else if ( fApplyGains == fgkGain ) 
  {
    

    Int_t serialNumber 
      = AliMpManuStore::Instance()->GetManuSerial(detElemId, manuId);
    
    AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fCapacitances->FindObject(serialNumber));
    
    if ( param )
    {
      capa = param->ValueAsFloat(manuChannel);
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
  
  if ( isSaturated ) 
  {
    Int_t saturation(3000);
  
    if ( gain && ( fApplyGains != fgkNoGain ) )
    {
      saturation = gain->ValueAsInt(manuChannel,4);
    }
  
    if ( padc >= saturation )
    {
      *isSaturated = kTRUE;
    }
    else
    {
      *isSaturated = kFALSE;
    }
  }
  
  return charge;
}

//_____________________________________________________________________________
Bool_t 
AliMUONDigitCalibrator::IsValidDigit(Int_t detElemId, Int_t manuId, Int_t manuChannel,
                                     Int_t* statusMap) const

{
  /// Check if a given pad is ok or not.
  
  // First a protection against bad input parameters
  AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
  if (!de) return kFALSE; // not existing DE
  if (!de->IsExistingChannel(manuId,manuChannel))
  {
    // non-existing (might happen when we get parity errors in read-out
    // that spoils the manuId
    return kFALSE;
  }
  if (!de->IsConnectedChannel(manuId,manuChannel))
  {
    // existing (in read-out), but not connected channel
    return kFALSE;
  }
  
  // ok, now we have a valid channel number, so let's see if that pad
  // behaves or not ;-)
  
  Int_t sm = fStatusMapMaker->StatusMap(detElemId,manuId,manuChannel);
  
  if (statusMap) *statusMap = sm;
  
  if ( ( sm & AliMUONPadStatusMapMaker::SelfDeadMask() ) != 0 ) 
  {
    // pad itself is bad (not testing its neighbours at this stage)
    return kFALSE;
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
Int_t 
AliMUONDigitCalibrator::PadStatus(Int_t detElemId, Int_t manuId, Int_t manuChannel) const
{
  /// Return the status of the given pad
  return fStatusMaker->PadStatus(detElemId,manuId,manuChannel);
}

//_____________________________________________________________________________
Int_t 
AliMUONDigitCalibrator::StatusMap(Int_t detElemId, Int_t manuId, Int_t manuChannel) const
{
  /// Return the status map of the given pad
  return fStatusMapMaker->StatusMap(detElemId,manuId,manuChannel);
  
}

