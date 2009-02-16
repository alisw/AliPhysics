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
#include "AliMUONCalibrationData.h"
#include "AliMUONTrackerData.h"
#include "AliMUONVStore.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpDEManager.h"
#include "AliMpDCSNamer.h"
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
	else if ( stype == "GAINS" ) 
	{
		AliMUONVStore* gains = AliMUONCalibrationData::CreateGains(runNumber,&startOfValidity);
		fData = CreateDataGains(startOfValidity);
		store = SplitQuality(*gains);
		delete gains;
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
	
	AliCDBManager::Instance()->SetDefaultStorage(storage);
	
	if (!store)
	{
		fIsValid = kFALSE;
		delete fData;
		fData = 0x0;
		AliError("Could not create store");
		return;
	}
	
	fData->Add(*store);
	
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
  
  AliMUONVTrackerData* data = new AliMUONTrackerData(Form("GAIN%d",runNumber),"Gains",6,kTRUE);
  data->SetDimensionName(0,"a1");
  data->SetDimensionName(1,"a2");
  data->SetDimensionName(2,"thres");
  data->SetDimensionName(3,"qual1");
  data->SetDimensionName(4,"qual2");
  data->SetDimensionName(5,"sat");
  return data;
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

//_____________________________________________________________________________
AliMUONVStore*
AliMUONTrackerOCDBDataMaker::SplitQuality(const AliMUONVStore& gains)
{
  /// Create a new store, identical to source gain store, except that qual 
  /// dimension is "decompacted" in two separated values
  
  AliMUONVStore* store = gains.Create();
  
  TIter next(gains.CreateIterator());
  AliMUONVCalibParam* param;
  
  while ( ( param = static_cast<AliMUONVCalibParam*>(next()) ) ) 
  {
    AliMUONVCalibParam* nd = new AliMUONCalibParamND(param->Dimension()+1,
                                                      param->Size(),
                                                      param->ID0(),
                                                      param->ID1());
    for ( Int_t i = 0; i < param->Size(); ++i ) 
    {
      for ( Int_t k = 0; k < param->Dimension(); ++k ) 
      {
        if ( k == 3 ) continue;
        Int_t m = ( k < 3 ? k : k+1 ) ;
        nd->SetValueAsDouble(i,m,param->ValueAsFloat(i,k));
      }
      Int_t qual = param->ValueAsInt(i,3);
			
			Int_t q1 = (qual & 0xF0) >> 4;  // linear fit quality
			Int_t q2 = qual & 0xF;		// parabolic fit quality
			
      nd->SetValueAsInt(i,3,q1);
      nd->SetValueAsInt(i,4,q2);
    }
    store->Add(nd);
  }
  return store;
}
