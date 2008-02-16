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

#include "AliCodeTimer.h"
#include "AliLog.h"
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
#include "AliMUONRawStreamTracker.h"

///\class AliMUONTrackerRawDataMaker
///
/// Creator of raw AliMUONVTrackerData from AliRawReader
/// 
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONTrackerRawDataMaker)
///\endcond

Int_t AliMUONTrackerRawDataMaker::fgkCounter(0);

//_____________________________________________________________________________
AliMUONTrackerRawDataMaker::AliMUONTrackerRawDataMaker(AliRawReader* reader, Bool_t histogram)
: AliMUONVTrackerDataMaker(),
  fRawReader(reader),
  fAccumulatedData(0x0),
  fOneEventData(new AliMUON2DMap(true)),
  fIsOwner(kTRUE),
  fSource("unspecified"),
  fIsRunning(kFALSE),
  fNumberOfEvents(0)
{
  /// Ctor
  reader->NextEvent(); // to be sure to get run number available
  
  Int_t runNumber = reader->GetRunNumber();
  
  TString name;
  
  if (!runNumber)
  {
    ++fgkCounter;    
    name = Form("%sRAW(%d)",(histogram?"H":""),fgkCounter);
  }
  else
  {
    name = Form("%sRAW%d",(histogram?"H":""),runNumber);
  }
  
  fAccumulatedData = new AliMUONTrackerData(name.Data(),"charge values",1);
  fAccumulatedData->SetDimensionName(0,"Raw charge");
  if ( histogram ) 
  {
    fAccumulatedData->MakeHistogramForDimension(0,kTRUE);
  }
  
  reader->RewindEvents();
}

//_____________________________________________________________________________
AliMUONTrackerRawDataMaker::~AliMUONTrackerRawDataMaker()
{
  /// dtor
  delete fOneEventData;
  if ( fIsOwner ) delete fAccumulatedData;
}

//_____________________________________________________________________________
Bool_t 
AliMUONTrackerRawDataMaker::NextEvent()
{
  /// Read next event
 
  AliCodeTimerAuto("");
  
  static Int_t nphysics(0);
  static Int_t ngood(0);

  fOneEventData->Clear();
  
  if ( !IsRunning() ) return kTRUE;
  
  Bool_t ok = fRawReader->NextEvent();

  if (!ok) 
  {
    return kFALSE;
  }
  
  Int_t eventType = fRawReader->GetType();

  ++fNumberOfEvents;
  
  if (eventType != AliRawEventHeaderBase::kPhysicsEvent ) 
  {
    return kTRUE; // for the moment
  }

  ++nphysics;

  AliMUONVRawStreamTracker* stream = new AliMUONRawStreamTracker(fRawReader);
    
  stream->First();
    
  Int_t buspatchId;
  UShort_t manuId;
  UChar_t manuChannel;
	UShort_t adc;
  
  while ( stream->Next(buspatchId,manuId,manuChannel,adc) )
  {    
    Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(buspatchId);
    
    AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fOneEventData->FindObject(detElemId,manuId));
    if (!param)
    {
      param = new AliMUONCalibParamND(1,64,detElemId,manuId,
                                      AliMUONVCalibParam::InvalidFloatValue());
      fOneEventData->Add(param);
    }
    
    param->SetValueAsDouble(manuChannel,0,adc);    
  }    
  
  if ( !stream->IsErrorMessage() )
  {
    ++ngood;
    fAccumulatedData->Add(*fOneEventData);
  }

  AliDebug(1,Form("n %10d nphysics %10d ngood %10d",fNumberOfEvents,nphysics,ngood));

  delete stream;
  
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
