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

#include "AliMUONTrackerCalibratedDataMaker.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCodeTimer.h"
#include "AliDAQ.h"
#include "AliLog.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONRawStreamTracker.h"
#include "AliMUONRawStreamTrackerHP.h"
#include "AliMUONRecoParam.h"
#include "AliMUONReconstructor.h"
#include "AliMUONTrackerData.h"
#include "AliMpDDLStore.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawReader.h"
#include <Riostream.h>

///\class AliMUONTrackerCalibratedDataMaker
///
/// Creator of AliMUONVTrackerData from AliRawReader
/// 
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONTrackerCalibratedDataMaker)
///\endcond

Int_t AliMUONTrackerCalibratedDataMaker::fgkCounter(0);

//_____________________________________________________________________________
AliMUONTrackerCalibratedDataMaker::AliMUONTrackerCalibratedDataMaker(TRootIOCtor*)
: AliMUONVTrackerDataMaker(),
fRawReader(0x0),
fIsOwnerOfRawReader(kFALSE),
fAccumulatedData(0x0),
fOneEventData(0x0),
fSource(""),
fIsRunning(kFALSE),
fDigitCalibrator(0x0),
fCalibrationData(0x0),
fCDBPath(""),
fNumberOfEvents(0),
fUseHPDecoder(kTRUE)
{
  /// Root IO ctor
}

//_____________________________________________________________________________
AliMUONTrackerCalibratedDataMaker::AliMUONTrackerCalibratedDataMaker(const AliMUONRecoParam* recoParam,
                                                                     Int_t runNumber,
                                                                     AliRawReader* reader, 
                                                                     const char* cdbpath,
                                                                     const char* calibMode,
                                                                     Bool_t histogram,
                                                                     Double_t xmin,
                                                                     Double_t xmax,                                                                     
                                                                     Bool_t useHPdecoder)
: AliMUONVTrackerDataMaker(),
fRawReader(reader),
fIsOwnerOfRawReader(kFALSE),
fAccumulatedData(0x0),
fOneEventData(new AliMUON2DMap(true)),
fSource("unspecified"),
fIsRunning(kFALSE),
fDigitCalibrator(0x0),
fCalibrationData(0x0),
fCDBPath(cdbpath),
fNumberOfEvents(0),
fUseHPDecoder(useHPdecoder)
{
  /// Ctor in which this object will NOT be owner of the reader,
  /// and can NOT apply rewind to it, nor use Next on it
  
  Ctor(recoParam,runNumber,calibMode,histogram,xmin,xmax);
}

//_____________________________________________________________________________
AliMUONTrackerCalibratedDataMaker::AliMUONTrackerCalibratedDataMaker(const AliMUONRecoParam* recoParam,
                                                                     AliRawReader* reader,
                                                                     const char* cdbpath,
                                                                     const char* calibMode,
                                                                     Bool_t histogram,
                                                                     Double_t xmin,
                                                                     Double_t xmax,
                                                                     Bool_t useHPDecoder)
: AliMUONVTrackerDataMaker(),
fRawReader(reader),
fIsOwnerOfRawReader(kTRUE),
fAccumulatedData(0x0),
fOneEventData(new AliMUON2DMap(true)),
fSource("unspecified"),
fIsRunning(kFALSE),
fDigitCalibrator(0x0),
fCalibrationData(0x0),
fCDBPath(cdbpath),
fNumberOfEvents(0),
fUseHPDecoder(useHPDecoder)
{
  /// Ctor, in which we are the owner of the reader, so we can rewind and advance it
  /// as we wish

  Int_t runNumber(0);
  
  if ( fRawReader ) 
  {
    fRawReader->NextEvent(); // to be sure to get run number available
    runNumber = reader->GetRunNumber();
    fRawReader->RewindEvents();
  }
  
  Ctor(recoParam,runNumber,calibMode,histogram,xmin,xmax);
}

//_____________________________________________________________________________
void
AliMUONTrackerCalibratedDataMaker::Ctor(const AliMUONRecoParam* recoParam, 
                                        Int_t runNumber, const char* calibMode,
                                        Bool_t histogram, Double_t xmin, Double_t xmax)
{
  /// "designated" constructor.

  ++fgkCounter;
  
  Bool_t calibrate = ( fCDBPath.Length() > 0 );
  TString name;
  TString basename("RAW");
  
  if ( calibrate ) 
  {
    TString scalib(calibMode);
    scalib.ToUpper();
    if ( scalib == "GAIN" ) basename = "CALC";
    if ( scalib == "NOGAIN" ) basename = "CALZ";
    if ( scalib == "GAINCONSTANTCAPA" ) basename = "CALG";
  }
  
  if (!runNumber)
  {
    name = Form("%s%s_%d",
                (histogram ? "H" : ""),
                basename.Data(),
                fgkCounter);
  }
  else
  {
    name = Form("%s%s%d",
                (histogram ? "H" : ""),
                basename.Data(),
                runNumber);
  }
  
  fAccumulatedData = new AliMUONTrackerData(name.Data(),"charge values",1);
  fAccumulatedData->SetDimensionName(0,(calibrate ? "Calibrated charge" : "Raw charge"));
  if ( histogram ) 
  {
    fAccumulatedData->MakeHistogramForDimension(0,kTRUE,xmin,xmax);
    AliInfo(Form("Will histogram between %e and %e",xmin,xmax));
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
    
    AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
    
    if ( storage->GetURI() != fCDBPath.Data() ) 
    {
      AliCDBManager::Instance()->SetDefaultStorage(fCDBPath.Data());
    }
    
    fCalibrationData->Pedestals();
    fCalibrationData->Gains();
    fCalibrationData->Neighbours();
    fCalibrationData->HV();
    fCalibrationData->Capacitances();
    
    if ( storage->GetURI() != fCDBPath.Data() ) 
    {
      AliCDBManager::Instance()->SetDefaultStorage(storage);
    }
    
    fDigitCalibrator = new AliMUONDigitCalibrator(*fCalibrationData,recoParam,calibMode);
		//FIXME: get the reco param from GUI and/or from OCDB if not used from the QA code ?
  }
}

//_____________________________________________________________________________
AliMUONTrackerCalibratedDataMaker::~AliMUONTrackerCalibratedDataMaker()
{
  /// dtor
  delete fOneEventData;
  delete fAccumulatedData;
  if ( fIsOwnerOfRawReader ) delete fRawReader;
  delete fCalibrationData;
  delete fDigitCalibrator;
}

//_____________________________________________________________________________
Long64_t 
AliMUONTrackerCalibratedDataMaker::Merge(TCollection*)
{
  /// Merge
  AliError("Not implemented yet");
  return 0;
}

//_____________________________________________________________________________
Bool_t 
AliMUONTrackerCalibratedDataMaker::NextEvent()
{
  /// Read and process next event
 
  if ( !fIsOwnerOfRawReader ) 
  {
    AliError("I'm not the owner of the raw reader. Cannot use NextEvent");
    return kFALSE;
  }
  
  AliCodeTimerAuto("");
  
  static Int_t nphysics(0);
  static Int_t ngood(0);

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

	Bool_t pok = ProcessEvent();
	
	if ( pok ) ++ngood;
	
	AliDebug(1,Form("n %10d nphysics %10d ngood %10d",fNumberOfEvents,nphysics,ngood));
	
	return kTRUE;
}

//_____________________________________________________________________________
Bool_t
AliMUONTrackerCalibratedDataMaker::ProcessEvent()
{
	/// Process current event 
  /// Note that we do not simply reuse the AliMUONDigitCalibrator::Calibrate(AliMUONVDigitStore&)
  /// method, as this would require filling first a digitStore, and then calibrate it, and
  /// then convert it into a VStore, all this taking too much time.
  ///
  /// But we *do* reuse the AliMUONDigitCalibrator::CalibrateDigit in order not to 
  /// duplicate this critical piece of calibration code !
  ///
  
  AliCodeTimerAuto("");

  AliMUONVRawStreamTracker* stream = 0x0;
  
  if ( fUseHPDecoder ) 
  {
    stream = new AliMUONRawStreamTrackerHP(fRawReader);
  }
  else
  {
    stream = new AliMUONRawStreamTracker(fRawReader);
  }
  
  stream->EnabbleErrorLogger();
  
  stream->First();
  
  Int_t buspatchId;
  UShort_t manuId;
  UChar_t manuChannel;
	UShort_t adc;
  const Int_t nddls = AliDAQ::NumberOfDdls("MUONTRK");
  TArrayI nevents(nddls);

  for ( Int_t i = 0; i < nddls; ++i ) 
  {
    nevents[i] = 0;
  }
  
  fOneEventData->Clear();
  
  while ( stream->Next(buspatchId,manuId,manuChannel,adc) )
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

    Bool_t ok = fDigitCalibrator->IsValidDigit(detElemId, manuId, manuChannel);

    if ( ok ) 
    {
      Float_t charge = fDigitCalibrator->CalibrateDigit(detElemId, manuId, manuChannel,adc,3.0);
    
      if (charge > 0.0 ) 
      {
        param->SetValueAsDouble(manuChannel,0,charge);
      }
    }
  }    
  
	Bool_t good(kFALSE);
	
  if ( !stream->IsErrorMessage() )
  {
    good = kTRUE;
    fAccumulatedData->Add(*fOneEventData,&nevents);    
  }
  
  delete stream;
  
  return good;  
}

//_____________________________________________________________________________
void
AliMUONTrackerCalibratedDataMaker::Print(Option_t*) const
{
  /// Printout
  
  cout << "Source=" << Source() << " Running=" << ( IsRunning() ? "YES" : "NO")
  << endl;
  
}

//_____________________________________________________________________________
void 
AliMUONTrackerCalibratedDataMaker::Rewind()
{
  /// Rewind events
  fRawReader->RewindEvents();
  fNumberOfEvents=0;
}

//_____________________________________________________________________________
void 
AliMUONTrackerCalibratedDataMaker::SetRawReader(AliRawReader* rawReader)
{
  /// Points to another raw reader
	
	if ( fIsOwnerOfRawReader ) 
	{
    AliFatal("Improper use of this class ! Cannot change raw reader in this case");
	}
	fRawReader = rawReader;
}
