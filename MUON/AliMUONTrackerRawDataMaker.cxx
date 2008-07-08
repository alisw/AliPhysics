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

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliLog.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONRawStreamTracker.h"
#include "AliMUONRawStreamTrackerHP.h"
#include "AliMUONTrackerData.h"
#include "AliMpDDLStore.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawReader.h"
#include <Riostream.h>

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
AliMUONTrackerRawDataMaker::AliMUONTrackerRawDataMaker(TRootIOCtor*)
: AliMUONVTrackerDataMaker(),
fRawReader(0x0),
fIsOwnerOfRawReader(kFALSE),
fAccumulatedData(0x0),
fOneEventData(0x0),
fSource(""),
fIsRunning(kFALSE),
fNumberOfEvents(0),
fRunNumber(0),
fIsEventByEvent(kFALSE),
fUseHPDecoder(kTRUE)
{
  /// Ctor
  ++fgkCounter;
}

//_____________________________________________________________________________
AliMUONTrackerRawDataMaker::AliMUONTrackerRawDataMaker(AliRawReader* reader, 
                                                       Bool_t histogram,
                                                       Bool_t useHPdecoder)
: AliMUONVTrackerDataMaker(),
fRawReader(reader),
fIsOwnerOfRawReader(kTRUE),
fAccumulatedData(0x0),
fOneEventData(new AliMUON2DMap(true)),
fSource("unspecified"),
fIsRunning(kFALSE),
fNumberOfEvents(0),
fRunNumber(0),
fIsEventByEvent(kFALSE),
fUseHPDecoder(useHPdecoder)
{
  /// Ctor
  
  if (fRawReader)
  {
    fRawReader->NextEvent(); // to be sure to get run number available
    fRunNumber = fRawReader->GetRunNumber();
    fRawReader->RewindEvents();
  }
    
  Ctor(histogram);
}

//_____________________________________________________________________________
AliMUONTrackerRawDataMaker::AliMUONTrackerRawDataMaker(Int_t runNumber,
                                                       AliRawReader* reader, 
                                                       Bool_t histogram,
                                                       Bool_t useHPdecoder)
: AliMUONVTrackerDataMaker(),
fRawReader(reader),
fIsOwnerOfRawReader(kTRUE),
fAccumulatedData(0x0),
fOneEventData(new AliMUON2DMap(true)),
fSource("unspecified"),
fIsRunning(kFALSE),
fNumberOfEvents(0),
fRunNumber(runNumber),
fIsEventByEvent(kFALSE),
fUseHPDecoder(useHPdecoder)
{
  /// Ctor
    
  Ctor(histogram);
}

//_____________________________________________________________________________
void
AliMUONTrackerRawDataMaker::Ctor(Bool_t histogram)
{
  /// Designated ctor
  
  TString name;
  
  if (!fRunNumber)
  {
    ++fgkCounter;    
    name = Form("%sRAW_%d",(histogram?"H":""),fgkCounter);
  }
  else
  {
    name = Form("%sRAW%d",(histogram?"H":""),fRunNumber);
  }
  
  fAccumulatedData = new AliMUONTrackerData(name.Data(),"charge values",1);
  fAccumulatedData->SetDimensionName(0,"Raw charge");
  if ( histogram ) 
  {
    fAccumulatedData->MakeHistogramForDimension(0,kTRUE);
  }
}

//_____________________________________________________________________________
AliMUONTrackerRawDataMaker::~AliMUONTrackerRawDataMaker()
{
  /// dtor
  delete fOneEventData;
  delete fAccumulatedData;
	if (fIsOwnerOfRawReader) delete fRawReader;
}

//_____________________________________________________________________________
Bool_t
AliMUONTrackerRawDataMaker::Add(const AliMUONTrackerRawDataMaker& other) 
{
  /// Adds other to this
  
//  AliRawReader* fRawReader; //!< reader of the data (owner or not)
//  Bool_t fIsOwnerOfRawReader; //!< whether we must delete rawReader or not
//  AliMUONVTrackerData* fAccumulatedData; ///< data (owner)
//  AliMUONVStore* fOneEventData; ///< data for one event (owner)
//  TString fSource; ///< where the data comes from
//  Bool_t fIsRunning; ///< whether we are running or are paused
//  Int_t fNumberOfEvents; ///< number of events seen
//  Int_t fRunNumber; ///< run number of the data
//  Bool_t fIsEventByEvent; ///< we only keep one event's data (no accumulation)
//  Bool_t fUseHPDecoder; ///< whether to use high performance decoder or not
//  static Int_t fgkCounter; ///< to count the number of instances
  
  if (!fAccumulatedData) return kFALSE;
  
  if ( fIsEventByEvent )
  {
    AliError("Cannot add event by event objects !");
    return kFALSE;
  }
  
  if ( fRunNumber != other.fRunNumber ) fRunNumber = -1;
  
  fSource += "\n";
  fSource += other.fSource;
  
  fNumberOfEvents += other.fNumberOfEvents;
  
  TList list;
  list.Add(other.fAccumulatedData);
  
  fAccumulatedData->Merge(&list);
  
  return kTRUE;
}

//_____________________________________________________________________________
Long64_t
AliMUONTrackerRawDataMaker::Merge(TCollection* list)
{
  /// Merge objects in collection
  
  if (!list) return 0;
  
  if ( list->IsEmpty() ) return NumberOfEvents();
  
  TIter next(list);
  const TObject* o(0x0);
  
  while ( ( o = next() ) )
  {
    const AliMUONTrackerRawDataMaker* data = dynamic_cast<const AliMUONTrackerRawDataMaker*>(o);
    if (!o)
    {
      AliError(Form("Object named %s is not an AliMUONTrackerRawDataMaker ! Skipping it",
                    o->GetName()));
    }
    else
    {
      Bool_t ok = Add(*data);
      if (!ok)
      {
        AliError("Got incompatible objects");
      }
    }
  }
  
  return NumberOfEvents();
}

//_____________________________________________________________________________
Bool_t 
AliMUONTrackerRawDataMaker::NextEvent()
{
  /// Read and process next event
 
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

	if ( ProcessEvent() )
	{
		++ngood;
	}
	
	AliDebug(1,Form("n %10d nphysics %10d ngood %10d",fNumberOfEvents,nphysics,ngood));
	
	return kTRUE;
}

//_____________________________________________________________________________
Bool_t 
AliMUONTrackerRawDataMaker::ProcessEvent()
{
	/// Process current event
	
  AliMUONVRawStreamTracker* stream = 0x0;
  
  if ( fUseHPDecoder ) 
  {
    stream = new AliMUONRawStreamTrackerHP(fRawReader);
  }
  else
  {
    stream = new AliMUONRawStreamTracker(fRawReader);
  }
    
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
  
	Bool_t good(kFALSE);
	
  if ( !stream->IsErrorMessage() )
  {
    good = kTRUE;
    fAccumulatedData->Add(*fOneEventData);
  }

  delete stream;
  
  return good;
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

//_____________________________________________________________________________
void 
AliMUONTrackerRawDataMaker::SetRawReader(AliRawReader* rawReader)
{
  /// Points to another raw reader
	
	if ( fIsOwnerOfRawReader ) 
	{
    AliFatal("Improper use of this class ! Cannot change raw reader in this case");
	}
	fRawReader = rawReader;
}
