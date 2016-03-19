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
#include "AliLog.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONConstants.h"
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
/// and multiplying by a constant gain, so that
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

//_____________________________________________________________________________
AliMUONDigitCalibrator::AliMUONDigitCalibrator(Int_t runNumber)
: TObject(),
fLogger(new AliMUONLogger(20000)),
fStatusMaker(0x0),
fStatusMapMaker(0x0),
fPedestals(0x0),
fNumberOfBadPads(0),
fNumberOfPads(0),
fChargeSigmaCut(0),
fMask(0)
{
  /// ctor
  
  AliMUONRecoParam* recoParam(0x0);
  
  AliCDBEntry* e = AliCDBManager::Instance()->Get("MUON/Calib/RecoParam",runNumber);
  if (e)
  {
    TObject* o = e->GetObject();
    if ( o->IsA() == TObjArray::Class() )
    {
      TObjArray* a = static_cast<TObjArray*>(o);
//      a->SetOwner(kTRUE); // FIXME: this should be done but somehow makes the reco crash at the end at cleaning stage... investigate why ?
      TIter next(a);
      AliMUONRecoParam* p;
      while ( ( p = static_cast<AliMUONRecoParam*>(next()) ))
      {
        if ( p->IsDefault()) recoParam = p;
      }
    }
    else
    {
      recoParam = static_cast<AliMUONRecoParam*>(o);
    }
  }
  if (!recoParam)
  {
    AliError("Cannot get the recoParam. Failing");
    return;
  }
  
  // OK. Now get all we need and work...
  
  AliMUONCalibrationData calib(runNumber);
  
  Ctor(calib,recoParam,kFALSE);
}

//_____________________________________________________________________________
AliMUONDigitCalibrator::AliMUONDigitCalibrator(const AliMUONCalibrationData& calib,
                                               const AliMUONRecoParam* recoParams)
: TObject(),
fLogger(new AliMUONLogger(20000)),
fStatusMaker(0x0),
fStatusMapMaker(0x0),
fPedestals(0x0),
fNumberOfBadPads(0),
fNumberOfPads(0),
fChargeSigmaCut(0),
fMask(0)
{
  /// ctor
  
  Ctor(calib,recoParams);
}

//_____________________________________________________________________________
AliMUONDigitCalibrator::AliMUONDigitCalibrator(const AliMUONCalibrationData& calib, int /*b*/)
: TObject(),
fLogger(new AliMUONLogger(20000)),
fStatusMaker(0x0),
fStatusMapMaker(0x0),
fPedestals(0x0),
fNumberOfBadPads(0),
fNumberOfPads(0),
fChargeSigmaCut(0),
fMask(0)
{
  /// ctor
  
  Ctor(calib,0x0);
}

//_____________________________________________________________________________
void
AliMUONDigitCalibrator::Ctor(const AliMUONCalibrationData& calib,
                             const AliMUONRecoParam* recoParams,
                             Bool_t deferredInitialization)
{
  /// designated ctor
  
  // Load mapping manu store
  if ( ! AliMpCDB::LoadManuStore() ) {
    AliFatal("Could not access manu store from OCDB !");
  }
  
  fStatusMaker = new AliMUONPadStatusMaker(calib);
  
  // Set default values, as loose as reasonable
  
  fChargeSigmaCut = 3.0;
  
	fMask = 0x8080; // reject pads where ped *or* hv are missing
	
	if ( recoParams )
	{
    // if we have reco params, we use limits and cuts from there :
    
    fStatusMaker->SetLimits(*recoParams);
    
    fMask = recoParams->PadGoodnessMask();
		//WARNING : getting this mask wrong is a very effective way of getting
		//no digits at all out of this class ;-)
    
    fChargeSigmaCut = recoParams->ChargeSigmaCut();
	}
  else
  {
    fLogger->Log("No RecoParam available");
    fLogger->Log(Form("SigmaCut=%e",fChargeSigmaCut));
  }
  
  fStatusMapMaker = new AliMUONPadStatusMapMaker(*fStatusMaker,fMask,deferredInitialization);
  
  fPedestals = calib.Pedestals();
}

//_____________________________________________________________________________
AliMUONDigitCalibrator::~AliMUONDigitCalibrator()
{
  /// dtor.
  
  if ( fNumberOfPads > 0 ) 
  {
    if ( fStatusMaker ) 
    {
      fStatusMaker->Report(fMask);
    }
    
    AliInfo("Summary of messages:");

    fLogger->Print();
    
    AliInfo(Form("We have seen %g pads, and rejected %g (%7.2f %%)",
                 fNumberOfPads,fNumberOfBadPads,
                 ( fNumberOfPads > 0 ) ? fNumberOfBadPads*100.0/fNumberOfPads : 0 ));
	}

  delete fStatusMaker;
  delete fStatusMapMaker;
  delete fLogger;
}

//_____________________________________________________________________________
void
AliMUONDigitCalibrator::Calibrate(AliMUONVDigitStore& digitStore)
{
  /// Calibrate the digits contained in digitStore  
  TIter next(digitStore.CreateTrackerIterator());
  AliMUONVDigit* digit;
  
  fStatusMapMaker->RefreshRejectProbabilities(); // this will do something only for simulations
  // (and only for those simulations where the reject list contain probabilities which are
  // different from zero or one)
  
  AliDebug(1,Form("# of digits = %d",digitStore.GetSize()));
  
  while ( ( digit = static_cast<AliMUONVDigit*>(next() ) ) )
  {
    if ( digit->IsCalibrated() ) 
    {
      fLogger->Log("ERROR : trying to calibrate a digit twice");
      return;
    }
    
    digit->Calibrated(kTRUE);
    digit->ChargeInFC(kTRUE);
    
    Float_t charge(0.0);
    Int_t statusMap;
    Bool_t isSaturated(kFALSE);
    
    ++fNumberOfPads;
    
    Bool_t ok = IsValidDigit(digit->DetElemId(),digit->ManuId(),digit->ManuChannel(),&statusMap);
    
    digit->SetStatusMap(statusMap);

    if (ok)
    {
      charge = CalibrateDigit(digit->DetElemId(),digit->ManuId(),digit->ManuChannel(),
                              digit->ADC(),fChargeSigmaCut,&isSaturated);
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
  /// Return the digit charge, in fC
  
  if ( nsigmas < 0 ) 
  {
    nsigmas = fChargeSigmaCut;
  }

  fLogger->Log(Form("ChargeSigmaCut used = %e",nsigmas));

  AliMUONVCalibParam* pedestal = static_cast<AliMUONVCalibParam*>
  (fPedestals->FindObject(detElemId,manuId));
  
  if (!pedestal)
  {
    // no pedestal -> no charge    
    fLogger->Log(Form("Got a null pedestal object for DE,manu=%d,%d",detElemId,manuId));        
    return 0.0;
  }
  
  Float_t padc = adc-pedestal->ValueAsFloat(manuChannel,0);
  
	// Gain (mV/fC) = 1/(a0*capa) with a0~1.25 and capa~0.2 
  Float_t charge(0);
  Float_t capa(AliMUONConstants::DefaultCapa()); // capa = 0.2 and a0 = 1.25
  Float_t a0(AliMUONConstants::DefaultA0());  // is equivalent to gain = 4 mV/fC
  Float_t adc2mv(AliMUONConstants::DefaultADC2MV()); // 1 ADC channel = 0.61 mV
  //
  // Note that the ChargeMax (for one pad) is roughly 4096 * 0.61 mV/channel / 4 mV/fC = 625 fC
  
  if ( padc > nsigmas*pedestal->ValueAsFloat(manuChannel,1) ) 
  {
    charge = a0*padc*capa*adc2mv;
  }
  
  if ( isSaturated ) 
  {
    const Int_t saturation(3000);
    
    if ( padc >= saturation )
    {
      *isSaturated = kTRUE;
    }
    else
    {
      *isSaturated = kFALSE;
    }
  }
  
  return ( charge > 0.0 ? charge : 0.0 );
}

//_____________________________________________________________________________
Bool_t 
AliMUONDigitCalibrator::IsValidDigit(Int_t detElemId, Int_t manuId, Int_t manuChannel,
                                     Int_t* statusMap) const

{
  /// Check if a given pad is ok or not.
  
  // initialize the statusmap to dead by default
  if (statusMap) *statusMap = AliMUONPadStatusMapMaker::SelfDeadMask();

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
  
  Int_t sm = StatusMap(detElemId,manuId,manuChannel);
  
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

