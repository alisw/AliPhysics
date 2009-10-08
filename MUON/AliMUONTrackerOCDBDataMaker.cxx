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

#include "AliMUONTrackerOCDBDataMaker.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONPadStatusMapMaker.h"
#include "AliMUONTrackerData.h"
#include "AliMUONVStore.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpDEManager.h"
#include "AliMpDCSNamer.h"
#include "AliMpManuIterator.h"
#include "Riostream.h"
#include <TClass.h>
#include <TMap.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>

///\class AliMUONTrackerOCDBDataMaker
///
/// Producer of AliMUONVTrackerData from OCDB data
///
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONTrackerOCDBDataMaker)
///\endcond

//_____________________________________________________________________________
AliMUONTrackerOCDBDataMaker::AliMUONTrackerOCDBDataMaker(const char* ocdbPath,
                                                           Int_t runNumber,
                                                           const char* type)
: AliMUONVTrackerDataMaker(),
  fIsValid(kTRUE),
  fData(0x0),
  fSource(Form("%s-%010d-%s",ocdbPath,runNumber,type))
{
	/// Ctor
	AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
	
	AliCDBManager::Instance()->SetDefaultStorage(ocdbPath);
	
	AliMUONVStore* store(0x0);
	
	TString stype(type);
	stype.ToUpper();
	Bool_t isSingleEvent(kTRUE);
	Int_t startOfValidity(0);
	
	if ( stype == "PEDESTALS" )
	{
		store = AliMUONCalibrationData::CreatePedestals(runNumber,&startOfValidity);
		fData = CreateDataPedestals(startOfValidity);
	}
	else if ( stype == "OCCUPANCY" )
	{
		store = AliMUONCalibrationData::CreateOccupancyMap(runNumber,&startOfValidity);
    
    if (store)
    {
      fData = new AliMUONTrackerData(Form("OCC%d",runNumber),"OccupancyMap",*store);
      fData->SetDimensionName(0,"One");
      fData->SetDimensionName(1,"Zero");
    }
	}
	else if ( stype == "GAINS" ) 
	{
    AliMUONVStore* gains = AliMUONCalibrationData::CreateGains(runNumber,&startOfValidity);
    store = PatchGainStore(*gains);
    delete gains;
		fData = CreateDataGains(startOfValidity);
	}
	else if ( stype == "CAPACITANCES" )
	{
		store = AliMUONCalibrationData::CreateCapacitances(runNumber,&startOfValidity);
		fData = CreateDataCapacitances(startOfValidity);
	}
	else if ( stype == "HV" )
	{
		TMap* m = AliMUONCalibrationData::CreateHV(runNumber,&startOfValidity);
		fData = new AliMUONTrackerData(Form("HV%d",startOfValidity),"High Voltages",1,!isSingleEvent);
		fData->SetDimensionName(0,"HV");
		store = CreateHVStore(*m);
		delete m;
	}
  else if ( stype == "STATUSMAP" )
  {
    fData = new AliMUONTrackerData(Form("STATUSMAP%d",runNumber),"Status map",2,kTRUE);
    fData->SetDimensionName(0,"Bits");
    fData->SetDimensionName(1,"Dead");
    store = CreateStatusMapStore(runNumber);
  }
  else if ( stype == "STATUS" )
  {
    fData = new AliMUONTrackerData(Form("STATUS%d",runNumber),"Status",1,kTRUE);
    fData->SetDimensionName(0,"Bits");
    store = CreateStatusStore(runNumber);
  }

	AliCDBManager::Instance()->SetDefaultStorage(storage);
	
	if (!store)
	{
		fIsValid = kFALSE;
		delete fData;
		fData = 0x0;
		AliError("Could not create store");
		return;
	}
	
  if ( stype != "OCCUPANCY" )
  {
    fData->Add(*store);
	}
  
	delete store;
}

//_____________________________________________________________________________
AliMUONTrackerOCDBDataMaker::~AliMUONTrackerOCDBDataMaker()
{
  /// dtor
  delete fData;
}

//_____________________________________________________________________________
AliMUONVTrackerData*
AliMUONTrackerOCDBDataMaker::CreateDataCapacitances(Int_t runNumber)
{
  /// Create data to hold capa values
  
  AliMUONVTrackerData* data = new AliMUONTrackerData(Form("CAPA%d",runNumber),"Capacitances",2,kTRUE);
  data->SetDimensionName(0,"Capa");
  data->SetDimensionName(1,"Injection gain");
  return data;
}

//_____________________________________________________________________________
AliMUONVTrackerData*
AliMUONTrackerOCDBDataMaker::CreateDataGains(Int_t runNumber)
{
  /// Create data to hold gains values
  
  AliMUONVTrackerData* data = new AliMUONTrackerData(Form("GAIN%d",runNumber),"Gains",7,kTRUE);
  data->SetDimensionName(0,"gain");
  data->SetDimensionName(1,"a1");
  data->SetDimensionName(2,"a2");
  data->SetDimensionName(3,"thres");
  data->SetDimensionName(4,"qual1");
  data->SetDimensionName(5,"qual2");
  data->SetDimensionName(6,"sat");
  return data;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONTrackerOCDBDataMaker::PatchGainStore(AliMUONVStore& gains)
{
  /// Polish the gain store : 
  /// a) adding a dimension, computed from a1, and called gain = 1/a1/0.2 
  ///     where 0.2 is internal capa in pF, and gain is then in mV/fC
  /// b) splitting the quality in two
  
  AliMUONVStore* store = gains.Create();
  
  TIter next(gains.CreateIterator());
  AliMUONVCalibParam* param;
  
  while ( ( param = static_cast<AliMUONVCalibParam*>(next()) ) ) 
  {
    AliMUONVCalibParam* nd = new AliMUONCalibParamND(param->Dimension()+2,
                                                     param->Size(),
                                                     param->ID0(),
                                                     param->ID1());
    for ( Int_t i = 0; i < param->Size(); ++i ) 
    {

      Int_t qual = param->ValueAsInt(i,3);
			Int_t q1 = (qual & 0xF0) >> 4;  // linear fit quality
			Int_t q2 = qual & 0xF;		// parabolic fit quality
			Double_t gain = 0.0;
      
      if ( param->ValueAsFloat(i,0) > 1E-9 ) gain = 1.0/param->ValueAsFloat(i,0)/0.2;
			
      nd->SetValueAsDouble(i,0,gain); // gain
      nd->SetValueAsDouble(i,1,param->ValueAsFloat(i,0)); // a1
      nd->SetValueAsDouble(i,2,param->ValueAsFloat(i,1)); // a2
      nd->SetValueAsInt(i,3,param->ValueAsInt(i,2)); // thres
      nd->SetValueAsInt(i,4,q1); // qual1
      nd->SetValueAsInt(i,5,q2); // qual2
      nd->SetValueAsInt(i,6,param->ValueAsInt(i,4)); // sat
    }
    store->Add(nd);
  }
  
  return store;
}

//_____________________________________________________________________________
AliMUONVTrackerData*
AliMUONTrackerOCDBDataMaker::CreateDataPedestals(Int_t runNumber)
{
  /// Create data to hold pedestal values
  
  AliMUONVTrackerData* data  = new AliMUONTrackerData(Form("PED%d",runNumber),"Pedestals",2,kTRUE);
  data->SetDimensionName(0,"Mean");
  data->SetDimensionName(1,"Sigma");
  return data;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONTrackerOCDBDataMaker::CreateStatusStore(Int_t runNumber)
{
  /// Get the status store
  
  AliMUONDigitCalibrator calibrator(runNumber);
  
  AliMUONVStore* sm = new AliMUON2DMap(kTRUE);
  
  AliMpManuIterator it;
  Int_t detElemId, manuId;
  
  while (it.Next(detElemId,manuId))
  {
    AliMUONVCalibParam* np = new AliMUONCalibParamNI(1,AliMpConstants::ManuNofChannels(),detElemId,manuId);
    for ( Int_t i = 0; i < np->Size(); ++i ) 
    {
      Int_t value = calibrator.PadStatus(detElemId,manuId,i);
      np->SetValueAsInt(i,0,value); // "raw" value of the status
    }
    sm->Add(np);
  }
  
  return sm;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONTrackerOCDBDataMaker::CreateStatusMapStore(Int_t runNumber)
{
  /// Get the status map, and polish it a bit for representation purposes

  AliMUONDigitCalibrator calibrator(runNumber);
  
  AliMUONVStore* sm = new AliMUON2DMap(kTRUE);
  
  AliMpManuIterator it;
  Int_t detElemId, manuId;
  
  while (it.Next(detElemId,manuId))
  {
    AliMUONVCalibParam* np = new AliMUONCalibParamNI(2,AliMpConstants::ManuNofChannels(),detElemId,manuId);
    for ( Int_t i = 0; i < np->Size(); ++i ) 
    {
      Int_t value = calibrator.StatusMap(detElemId,manuId,i);
      Int_t channelIsDead = ( value & AliMUONPadStatusMapMaker::SelfDeadMask() );
      np->SetValueAsInt(i,0,value); // "raw" value of the status map
      np->SetValueAsInt(i,1,channelIsDead); // simple 0 or 1 for this channel
    }
    sm->Add(np);
  }
  
  return sm;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONTrackerOCDBDataMaker::CreateHVStore(TMap& m)
{
  /// Create a store from hv values
  
  AliMUONVStore* store = new AliMUON2DMap(kTRUE);
  
  TIter next(&m);
  TObjString* s;
  AliMpDCSNamer hvNamer("TRACKER");
  
  while ( ( s = static_cast<TObjString*>(next()) ) )
  {
    TString name(s->String());

    Int_t detElemId = hvNamer.DetElemIdFromDCSAlias(name.Data());
    
    if ( !AliMpDEManager::IsValidDetElemId(detElemId) )
    {
      AliErrorClass(Form("Got an invalid DE = %d from alias = %s",
                         detElemId,name.Data()));
      continue;
    }
    
    Int_t nindex = 1;
    Int_t hvIndex = hvNamer.DCSIndexFromDCSAlias(name.Data());
    
    if ( hvIndex > 0 && detElemId >= 500 ) 
    {
      AliFatalClass("FIXME"); // there's now switch aliases which should be taken into account
    }

    if ( hvIndex == -2 ) // we should consider switch alias there...
    {
      nindex = hvNamer.NumberOfPCBs(detElemId);
      hvIndex = 0;
    }

    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

    for ( int i = 0 ; i < nindex; ++i )
    {
      Int_t index = hvIndex + i ;
      
      const AliMpArrayI* manus = de->ManusForHV(index);
      
      TPair* p = static_cast<TPair*>(m.FindObject(name.Data()));
      TObjArray* a = static_cast<TObjArray*>(p->Value());
      TIter n2(a);
      AliDCSValue* v;
      Float_t hvValue(0);
      Int_t n(0);
      while ( ( v = static_cast<AliDCSValue*>(n2()) ) )
      {
        hvValue += v->GetFloat();  
        ++n;
      }
      if ( n ) hvValue /= n;
      
      Int_t nofChannels(AliMpConstants::ManuNofChannels());
      
      for ( Int_t k = 0 ; k < manus->GetSize(); ++k )
      {
        Int_t manuId = manus->GetValue(k);
        AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(store->FindObject(detElemId,manuId));
        if ( ! param ) 
        {
          param = new AliMUONCalibParamND(1,nofChannels,detElemId,manuId,0);
          store->Add(param);
        }
        for ( Int_t j = 0 ; j < nofChannels; ++j )
        {
          param->SetValueAsDouble(j,0,hvValue);
        }
      }
    }
  }
  
  return store;
  
}

//_____________________________________________________________________________
Long64_t 
AliMUONTrackerOCDBDataMaker::Merge(TCollection*)
{
  /// Merge
  AliError("Not implemented. Does it have sense ?");
  return 0;
}
