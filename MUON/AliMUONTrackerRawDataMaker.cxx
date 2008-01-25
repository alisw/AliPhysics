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

#include "AliMUONTrackerRawDataMaker.h"

#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONDigitMaker.h"
#include "AliMUONDigitStoreV2R.h"
#include "AliMUONTrackerData.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMpDDLStore.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include <Riostream.h>

///\class AliMUONTrackerRawDataMaker
///
/// Creator of AliMUONVTrackerData from AliRawReader
/// 
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONTrackerRawDataMaker)
///\endcond

Int_t AliMUONTrackerRawDataMaker::fgkCounter(0);

//_____________________________________________________________________________
AliMUONTrackerRawDataMaker::AliMUONTrackerRawDataMaker(AliRawReader* reader,
                                                       const char* cdbpath)
: AliMUONVTrackerDataMaker(),
  fRawReader(reader),
  fAccumulatedData(0x0),
  fOneEventData(new AliMUON2DMap(true)),
  fIsOwner(kTRUE),
  fSource("unspecified"),
  fIsRunning(kFALSE),
  fDigitMaker(0x0),
  fDigitCalibrator(0x0),
  fCalibrationData(0x0),
  fDigitStore(0x0), 
  fCDBPath(cdbpath),
  fNumberOfEvents(0)
{
  /// Ctor
  reader->NextEvent(); // to be sure to get run number available
  
  Int_t runNumber = reader->GetRunNumber();
  
  ++fgkCounter;
  
  Bool_t calibrate = ( fCDBPath.Length() > 0 );
  
  TString name;
  
  if (!runNumber)
  {
    name = Form("%s(%d)",(calibrate ? "CAL" : "RAW"),fgkCounter);
  }
  else
  {
    name = Form("%s%d",(calibrate ? "CAL" : "RAW"),runNumber);
  }
  
  fAccumulatedData = new AliMUONTrackerData(name.Data(),"charge values",1);
  fAccumulatedData->SetDimensionName(0,(calibrate ? "Calibrated charge" : "Raw charge"));
  
  reader->RewindEvents();

  fDigitMaker = new AliMUONDigitMaker;
  fDigitMaker->SetMakeTriggerDigits(kFALSE);
  fDigitStore = new AliMUONDigitStoreV2R;

  if ( calibrate ) 
  {
    fCalibrationData = new AliMUONCalibrationData(runNumber);
    
    // force the reading of calibration NOW
    // FIXME: not really elegant and error prone (as we have the list of calib data twice, 
    // once here and once in the digitcalibrator class, hence the change of them getting
    // out of sync)
    // But with the current CDBManager implementation, I don't know how to solve
    // this better (e.g. to avoid clearing cache messages and so on).
    
    AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
    
    if ( storage->GetURI() != fCDBPath.Data() ) 
    {
      AliCDBManager::Instance()->SetDefaultStorage(fCDBPath.Data());
    }
    
    fCalibrationData->Pedestals();
    fCalibrationData->Gains();
    fCalibrationData->Neighbours();
    fCalibrationData->HV();
    
    if ( storage->GetURI() != fCDBPath.Data() ) 
    {
      AliCDBManager::Instance()->SetDefaultStorage(storage);
    }
    
    fDigitCalibrator = new AliMUONDigitCalibrator(*fCalibrationData);
  }
}

//_____________________________________________________________________________
AliMUONTrackerRawDataMaker::~AliMUONTrackerRawDataMaker()
{
  /// dtor
  delete fOneEventData;
  if ( fIsOwner ) delete fAccumulatedData;
  delete fRawReader;
  delete fDigitStore;
  delete fCalibrationData;
  delete fDigitMaker;
  delete fDigitCalibrator;
}

//_____________________________________________________________________________
Bool_t 
AliMUONTrackerRawDataMaker::NextEvent()
{
  /// Read next event
 
  static Int_t nphysics(0);
  static Int_t ngood(0);

  if ( !IsRunning() ) return kTRUE;
  
  Bool_t ok = fRawReader->NextEvent();

  if (!ok) 
  {
    fDigitMaker->Print();
    return kFALSE;
  }
  
  Int_t eventType = fRawReader->GetType();

  ++fNumberOfEvents;
  
  if (eventType != AliRawEventHeaderBase::kPhysicsEvent ) 
  {
    return kTRUE; // for the moment
  }

  ++nphysics;

  Int_t rv = fDigitMaker->Raw2Digits(fRawReader,fDigitStore);
  
  if ( ( rv & AliMUONDigitMaker::kTrackerBAD ) != 0 ) return kTRUE;

  if ( fDigitCalibrator ) 
  {
    fDigitCalibrator->Calibrate(*fDigitStore);
  }
  
  Bool_t dok = ConvertDigits();
  
  if ( dok )
    {
      ++ngood;
      fAccumulatedData->Add(*fOneEventData);
    }

  AliDebug(1,Form("n %10d nphysics %10d ngood %10d",fNumberOfEvents,nphysics,ngood));

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t 
AliMUONTrackerRawDataMaker::ConvertDigits()
{
  /// Convert digitstore into fOneEventData
  
  TIter next(fDigitStore->CreateIterator());
  AliMUONVDigit* digit;

  fOneEventData->Clear();
  
  while ( ( digit = static_cast<AliMUONVDigit*>(next())) )
  {
    Double_t value = ( digit->IsCalibrated() ? digit->Charge() : digit->ADC() );

    if ( value > 0 ) 
    {
      Int_t detElemId = digit->DetElemId();
      Int_t manuId = digit->ManuId();
    
      AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fOneEventData->FindObject(detElemId,manuId));
      if (!param)
      {
        param = new AliMUONCalibParamND(1,64,detElemId,manuId,
                                      AliMUONVCalibParam::InvalidFloatValue());
        fOneEventData->Add(param);
      }
    
      param->SetValueAsDouble(digit->ManuChannel(),0,value);
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
void
AliMUONTrackerRawDataMaker::Print(Option_t*) const
{
  /// Printout
  
  cout << "Source=" << Source() << " Running=" << ( IsRunning() ? "YES" : "NO")
  << endl;
  
}

//_____________________________________________________________________________
void 
AliMUONTrackerRawDataMaker::Rewind()
{
  /// Rewind events
  fRawReader->RewindEvents();
  fNumberOfEvents=0;
}
