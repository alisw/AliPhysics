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

/// \class AliMUONTrackerDataMaker
/// 
/// Implementation of VTrackerDataMaker to read raw data and 
/// calibrate it (if required)
/// 
/// \author Laurent Aphecetche, Subatech

#include "AliMUONTrackerDataMaker.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCodeTimer.h"
#include "AliDAQ.h"
#include "AliLog.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONLogger.h"
#include "AliMUONRawStreamTrackerHP.h"
#include "AliMUONTrackerData.h"
#include "AliMpDDLStore.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawReader.h"
#include "Riostream.h"

/// \cond CLASSIMP
ClassImp(AliMUONTrackerDataMaker)
/// \endcond

Int_t AliMUONTrackerDataMaker::fgkCounter(0);

//_____________________________________________________________________________
AliMUONTrackerDataMaker::AliMUONTrackerDataMaker(TRootIOCtor*) 
: 
AliMUONVTrackerDataMaker(),
fRawReader(0x0),
fAccumulatedData(0x0),
fIsOwnerOfAccumulatedData(kTRUE),
fOneEventData(0x0),
fDigitCalibrator(0x0),
fCalibrationData(0x0), 
fSource(""),
fOCDBPath(""),
fNumberOfEvents(0),
fRunNumber(0),
fIsRunning(kFALSE),
fIsOwnerOfRawReader(kFALSE),
fIsEventByEvent(kFALSE),
fLogger(0x0),
fLastEventWasEmpty(kFALSE),
fNumberOfPhysicsEvents(0),
fNumberOfGoodPhysicsEvents(0),
fTryRecover(kFALSE),
fFirstEvent(-1),
fLastEvent(-1)
{
/// Root IO ctor
}

//_____________________________________________________________________________
AliMUONTrackerDataMaker::AliMUONTrackerDataMaker(const AliMUONRecoParam* recoParam,
                                                 Int_t runNumber,
                                                 AliRawReader* rawReader,
                                                 const char* cdbPath,
                                                 const char* calibMode,
                                                 Bool_t histogram,
                                                 Double_t xmin,
                                                 Double_t xmax)
:
AliMUONVTrackerDataMaker(),
fRawReader(rawReader),
fAccumulatedData(0x0),
fIsOwnerOfAccumulatedData(kTRUE),
fOneEventData(new AliMUON2DMap(true)),
fDigitCalibrator(0x0),
fCalibrationData(0x0), 
fSource(""),
fOCDBPath(cdbPath),
fNumberOfEvents(0),
fRunNumber(runNumber),
fIsRunning(kFALSE),
fIsOwnerOfRawReader(kFALSE),
fIsEventByEvent(kFALSE),
fLogger(0x0),
fLastEventWasEmpty(kFALSE),
fNumberOfPhysicsEvents(0),
fNumberOfGoodPhysicsEvents(0),
fTryRecover(kFALSE),
fFirstEvent(-1),
fLastEvent(-1)
{
  /// Ctor in which this object will NOT be the owner of the reader
  /// and can NOT apply rewind to it, nor use Next on it. 
  Ctor(recoParam,runNumber,calibMode,histogram,xmin,xmax);
}


//_____________________________________________________________________________
AliMUONTrackerDataMaker::AliMUONTrackerDataMaker(const AliMUONRecoParam* recoParam,
                                                 AliRawReader* rawReader,
                                                 const char* cdbPath,
                                                 const char* calibMode,
                                                 Bool_t histogram,
                                                 Double_t xmin,
                                                 Double_t xmax)
:
AliMUONVTrackerDataMaker(),
fRawReader(rawReader),
fAccumulatedData(0x0),
fIsOwnerOfAccumulatedData(kTRUE),
fOneEventData(new AliMUON2DMap(true)),
fDigitCalibrator(0x0),
fCalibrationData(0x0), 
fSource(""),
fOCDBPath(cdbPath),
fNumberOfEvents(0),
fRunNumber(0),
fIsRunning(kFALSE),
fIsOwnerOfRawReader(kTRUE),
fIsEventByEvent(kFALSE),
fLogger(0x0),
fLastEventWasEmpty(kFALSE),
fNumberOfPhysicsEvents(0),
fNumberOfGoodPhysicsEvents(0),
fTryRecover(kFALSE),
fFirstEvent(-1),
fLastEvent(-1)
{
  /// Ctor in which we take the ownership of the rawReader, so we can rewind
  /// and advance it as we wish
  
  if ( fRawReader && fRawReader->NextEvent() ) 
  {
    fRunNumber = fRawReader->GetRunNumber();
    fRawReader->RewindEvents();
  }
  
  Ctor(recoParam,fRunNumber,calibMode,histogram,xmin,xmax);
}

//_____________________________________________________________________________
AliMUONTrackerDataMaker::AliMUONTrackerDataMaker(AliRawReader* rawReader, Bool_t histogram)
:
AliMUONVTrackerDataMaker(),
fRawReader(rawReader),
fAccumulatedData(0x0),
fIsOwnerOfAccumulatedData(kTRUE),
fOneEventData(new AliMUON2DMap(true)),
fDigitCalibrator(0x0),
fCalibrationData(0x0), 
fSource(""),
fOCDBPath(""),
fNumberOfEvents(0),
fRunNumber(0),
fIsRunning(kFALSE),
fIsOwnerOfRawReader(kTRUE),
fIsEventByEvent(kFALSE),
fLogger(0x0),
fLastEventWasEmpty(kFALSE),
fNumberOfPhysicsEvents(0),
fNumberOfGoodPhysicsEvents(0),
fTryRecover(kFALSE),
fFirstEvent(-1),
fLastEvent(-1)
{
  /// Ctor from raw data reader
  if ( fRawReader && fRawReader->NextEvent() ) 
  {
    fRunNumber = fRawReader->GetRunNumber();
    fRawReader->RewindEvents();
  }
  
  Ctor(0x0,fRunNumber,"",histogram);
  
}

//_____________________________________________________________________________
void 
AliMUONTrackerDataMaker::Ctor(const AliMUONRecoParam* recoParam,
                              Int_t runNumber,
                              const char* calibMode,
                              Bool_t histogram,
                              Double_t xmin, Double_t xmax)
{
  /// "designated constructor"

  Bool_t calibrate = ( strlen(calibMode) > 0 );
  
  TString name;
  TString type("RAW");
  
  if ( calibrate ) 
  {
    TString scalib(calibMode);
    scalib.ToUpper();
    if ( scalib == "GAIN" ) type = "CALC";
    if ( scalib == "NOGAIN" ) type = "CALZ";
    if ( scalib == "GAINCONSTANTCAPA") type = "CALG";
    if ( scalib == "INJECTIONGAIN" ) type = "CALE";
  }
  
  if ( !fRunNumber ) 
  {
    ++fgkCounter;
    name = Form("%s%s_%d",(histogram?"H":""),type.Data(),fgkCounter);
  }
  else
  {
    name = Form("%s%s%d",(histogram?"H":""),type.Data(),fRunNumber);
  }
  
  fAccumulatedData = new AliMUONTrackerData(name.Data(),"charge values",1);
  fAccumulatedData->SetDimensionName(0,(calibrate ? "Calibrated charge" : "Raw charge"));
  if (histogram)
  {
    fAccumulatedData->MakeHistogramForDimension(0,kTRUE,xmin,xmax);
  }
  
  if ( calibrate ) 
  {
    fCalibrationData = new AliMUONCalibrationData(runNumber);
    
    // force the reading of calibration NOW
    // FIXME: not really elegant and error prone (as we have the list of calib data twice, 
    // once here and once in the digitcalibrator class, hence the change of them getting
    // out of sync)
    // But with the current CDBManager implementation, I don't know how to solve
    // this better (e.g. to avoid clearing cache messages and so on).

    AliCDBStorage* storage(0x0);
    
    if ( fOCDBPath.Length() > 0 )
    {
      storage = AliCDBManager::Instance()->GetDefaultStorage();

      if ( storage && ( storage->GetURI() != fOCDBPath.Data() ) )
      {
        AliCDBManager::Instance()->SetDefaultStorage(fOCDBPath.Data());
      }
    }
    
    fCalibrationData->Pedestals();
    fCalibrationData->Gains();
    fCalibrationData->Neighbours();
    fCalibrationData->HV();
    fCalibrationData->Capacitances();
    
    if ( storage && ( storage->GetURI() != fOCDBPath.Data() ) )
    {
      AliCDBManager::Instance()->SetDefaultStorage(storage);
    }
    
    fDigitCalibrator = new AliMUONDigitCalibrator(*fCalibrationData,recoParam);
    //FIXME: get the reco param from GUI and/or from OCDB if not used from the QA code ?
  }
}

//_____________________________________________________________________________
AliMUONTrackerDataMaker::~AliMUONTrackerDataMaker()
{
/// dtor

  delete fOneEventData;
  if ( fIsOwnerOfAccumulatedData ) delete fAccumulatedData;
  if ( fIsOwnerOfRawReader ) delete fRawReader;
  delete fCalibrationData;
  delete fDigitCalibrator;
}

//_____________________________________________________________________________
Bool_t 
AliMUONTrackerDataMaker::Add(const AliMUONTrackerDataMaker& other)
{
  /// Adds other to this
    
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
  fNumberOfPhysicsEvents += other.fNumberOfPhysicsEvents;
  fNumberOfGoodPhysicsEvents += other.fNumberOfGoodPhysicsEvents;

  TList list;
  list.Add(other.fAccumulatedData);
  
  fAccumulatedData->Merge(&list);
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t 
AliMUONTrackerDataMaker::NextEvent()
{
  /// Read and process next event
  
  if ( !fIsOwnerOfRawReader ) 
  {
    AliError("I'm not the owner of the raw reader. Cannot use NextEvent");
    return kFALSE;
  }
  
  AliCodeTimerAuto("",0);
  
  if ( !IsRunning() ) return kTRUE;
  
  Bool_t ok(kTRUE);
  
  if ( fLastEvent >= fFirstEvent && fLastEvent > 0 ) // do we have an event range to consider ?
  {
    // skip up to first event
    
    while ( (fNumberOfEvents-1) < fFirstEvent && ( ok = fRawReader->NextEvent() ) ) 
    {
      ++fNumberOfEvents; 
    }
    
    if ( ok && (fNumberOfEvents-1) <= fLastEvent ) 
    {
      ok = fRawReader->NextEvent();
    }
    else
    {
      fNumberOfEvents=fLastEvent+1;
      return kFALSE;
    }
  }
  else
  {
    // no event range, just proceed...
    ok = fRawReader->NextEvent();
  }
  
  if (!ok) 
  {
    return kFALSE;
  }
  
	ProcessEvent();
	
	return kTRUE;  
}

//_____________________________________________________________________________
Bool_t AliMUONTrackerDataMaker::ProcessEvent()
{
  /// Process current event 
  /// 
  /// Note that in case of calibration, we do not simply reuse the 
  /// AliMUONDigitCalibrator::Calibrate(AliMUONVDigitStore&) method, 
  /// as this would require filling first a digitStore, and then calibrate it,
  /// and then convert it into a VStore, all this taking too much time.
  /// But we *do* reuse the AliMUONDigitCalibrator::CalibrateDigit in order not to 
  /// duplicate this critical piece of calibration code !
  ///
  
  ++fNumberOfEvents;
  
  Int_t eventType = fRawReader->GetType();
  
  if (eventType != AliRawEventHeaderBase::kPhysicsEvent ) 
  {
    return kTRUE; // for the moment
  }
  
  ++fNumberOfPhysicsEvents;
  
  fLastEventWasEmpty = kFALSE;
  
  AliCodeTimerAuto("",0);
  
  AliMUONRawStreamTrackerHP stream(fRawReader);

  stream.DisableWarnings();
  stream.DisableRawReaderErrorLogger();
  stream.DisableMUONErrorLogger();

  if ( fTryRecover ) 
  {
    stream.TryRecover(kTRUE);
  }
  else
  {
    stream.TryRecover(kFALSE);
  }
  
  if (fLogger)
  {
    stream.EnableMUONErrorLogger();  
    stream.SetMUONErrorLogger(fLogger);    
    stream.SetLoggingDetailLevel(AliMUONRawStreamTrackerHP::kMediumErrorDetail);
  }
  
  const Int_t nddls = AliDAQ::NumberOfDdls("MUONTRK");
  TArrayI nevents(nddls);
  
  for ( Int_t i = 0; i < nddls; ++i ) 
  {
    nevents[i] = 0;
  }
  
  fOneEventData->Clear();
  
  Int_t buspatchId;
  UShort_t  manuId;
  UChar_t manuChannel;
  UShort_t adc;
  
  stream.First();
  
  while ( stream.Next(buspatchId,manuId,manuChannel,adc,kTRUE) )
  {    
    Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(buspatchId);
    
    Int_t ddl = AliMpDDLStore::Instance()->GetDDLfromBus(buspatchId);
    
    nevents[ddl] = 1;
    
    AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(fOneEventData->FindObject(detElemId,manuId));
    if (!param)
    {
      param = new AliMUONCalibParamND(1,64,detElemId,manuId,
                                      AliMUONVCalibParam::InvalidFloatValue());
      fOneEventData->Add(param);
    }
    
    Double_t charge(adc);
    
    if ( fDigitCalibrator ) 
    {
      if ( fDigitCalibrator->IsValidDigit(detElemId, manuId, manuChannel) )
      {
        charge = fDigitCalibrator->CalibrateDigit(detElemId, manuId, manuChannel,adc);
      }
      else
      {
        charge = 0.0;
      }
    }
    
    if (charge > 0.0 ) 
    {
      param->SetValueAsDouble(manuChannel,0,charge);
    }
  }
  
	Bool_t badEvent = stream.HasPaddingError() || stream.HasGlitchError();

	if (!badEvent)
  {
    fAccumulatedData->Add(*fOneEventData,&nevents);  
    if ( fOneEventData->GetSize() == 0 ) fLastEventWasEmpty = kTRUE;
    ++fNumberOfGoodPhysicsEvents;
  }
    
  AliDebug(1,Form("n %10d nphysics %10d ngood %10d",fNumberOfEvents,fNumberOfPhysicsEvents,fNumberOfGoodPhysicsEvents));
	
  return !badEvent;
}


//_____________________________________________________________________________
void 
AliMUONTrackerDataMaker::Print(Option_t*) const
{
  /// Printout
  
  cout << "Source=" << Source() << " Running=" << ( IsRunning() ? "YES" : "NO")
  << endl;
}

//_____________________________________________________________________________
void AliMUONTrackerDataMaker::Rewind()
{
  /// Rewind events
  if ( fIsOwnerOfRawReader ) 
  {
    fRawReader->RewindEvents();
    fNumberOfEvents=0;  
    fNumberOfPhysicsEvents=0;
    fNumberOfGoodPhysicsEvents=0;
  }
  else
  {
    AliError("Wrong usage of this class : cannot rewind as I am not owner of the raw reader !");
  }
}

//_____________________________________________________________________________
Long64_t AliMUONTrackerDataMaker::Merge(TCollection* list)
{
  /// Merge objects in collection
  
  if (!list) return 0;
  
  if ( list->IsEmpty() ) return NumberOfEvents();
  
  TIter next(list);
  const TObject* o(0x0);
  
  while ( ( o = next() ) )
  {
    const AliMUONTrackerDataMaker* data = dynamic_cast<const AliMUONTrackerDataMaker*>(o);
    if (!data)
    {
      AliError(Form("Object named %s is not an AliMUONTrackerDataMaker ! Skipping it",
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
void 
AliMUONTrackerDataMaker::SetRawReader(AliRawReader* rawReader)
{
  /// Change the rawreader (only works if isowner=true)
  
  if ( fIsOwnerOfRawReader ) 
	{
    AliFatal("Improper use of this class ! Cannot change raw reader in this case");
	}
	
  fRawReader = rawReader;
  
}
