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

#include "AliMUONTrackerConditionDataMaker.h"

///\class AliMUONTrackerConditionDataMaker
///
/// Producer of AliMUONVTrackerData from OCDB or Ascii file condition data
///
/// \author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONTrackerConditionDataMaker)
///\endcond

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliMpManuIterator.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONCalibParamNF.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONPadStatusMapMaker.h"
#include "AliMUONTrackerData.h"
#include "AliMUONTrackerIO.h"
#include "AliMpArrayI.h"
#include "AliMpConstants.h"
#include "AliMpDCSNamer.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "TClass.h"
#include "TMap.h"
#include "TObjString.h"
#include "Riostream.h"
#include "TString.h"
#include <sstream>
#include "TSystem.h"

//_____________________________________________________________________________
AliMUONTrackerConditionDataMaker::AliMUONTrackerConditionDataMaker():
AliMUONVTrackerDataMaker(),
fData(0x0),
fSource("")
{
  /// default ctor to be able to stream
}

//_____________________________________________________________________________
AliMUONTrackerConditionDataMaker::AliMUONTrackerConditionDataMaker(Int_t runNumber, const char* ocdbPath, const char* type):
AliMUONVTrackerDataMaker(),
fData(0x0),
fSource(Form("%s-%010d-%s",ocdbPath,runNumber,type))
{
  /// ctor from OCDB

  AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
	
	AliCDBManager::Instance()->SetDefaultStorage(ocdbPath);

  Int_t startOfValidity;
  AliMUONVStore* store = CreateStore(runNumber,ocdbPath,type,startOfValidity);
  AliDebug(1,Form("runNumber=%d ocdbPath=%s type=%s startOfValidity=%d store=%p",
                  runNumber,ocdbPath,type,startOfValidity,store));
  if ( store )
  {
    fData = CreateData(type,*store,startOfValidity);
  }  

  delete store;
  
  AliCDBManager::Instance()->SetDefaultStorage(storage);
}

//_____________________________________________________________________________
AliMUONTrackerConditionDataMaker::AliMUONTrackerConditionDataMaker(const char* filename, const char* type):
AliMUONVTrackerDataMaker(),
fData(0x0),
fSource(Form("%s-%s",filename,type))
{
  /// ctor from an ASCII file
  
  TString sFilename(gSystem->ExpandPathName(filename));
  
  std::ifstream in(sFilename.Data());
  if (in.good()) 
  {
    std::ostringstream stream;
    char line[1024];
    while ( in.getline(line,1024) )
    {
      stream << line << "\n";    
    }
  
    in.close();
  
    Int_t dummy;
    
    AliMUONVStore* store = CreateStore(-1,stream.str().c_str(),type,dummy);
    
    if ( store )
    {
      fData = CreateData(type,*store,dummy);
    }
    delete store;
  }
}

//_____________________________________________________________________________
AliMUONTrackerConditionDataMaker::AliMUONTrackerConditionDataMaker(const char* data, const char* type, Bool_t) :
AliMUONVTrackerDataMaker(),
fData(0x0),
fSource(Form("direct-%s",type))
{
  /// ctor from a string containing the ASCII data
  /// the last parameter is there just to distinguish this ctor from the previous one
  
  Int_t dummy;
  
  AliMUONVStore* store = CreateStore(-1,data,type,dummy);
  
  if ( store )
  {
    fData = CreateData(type,*store,dummy);
  }
  delete store;
  
}

//_____________________________________________________________________________
AliMUONTrackerConditionDataMaker::~AliMUONTrackerConditionDataMaker()
{
  /// dtor
  delete fData;
}


//_____________________________________________________________________________
AliMUONVTrackerData*
AliMUONTrackerConditionDataMaker::CreateData(const char* type, AliMUONVStore& store, Int_t startOfValidity)
{
  /// Create the data source 
  AliMUONVTrackerData* data(0x0);
  
  TString stype(type);
  stype.ToUpper();
  
  if ( stype == "CAPACITANCES" )
  {    
    data = new AliMUONTrackerData(Form("CAPA%d",startOfValidity),"Capacitances",2,kTRUE);
    data->SetDimensionName(0,"Capa");
    data->SetDimensionName(1,"Injection gain");    
  }
  else if ( stype == "CONFIG" ) 
  {
    data = new AliMUONTrackerData(Form("CONFIG%d",startOfValidity),"Configuration",1);
    data->SetDimensionName(0,"there");
    data->DisableChannelLevel();
  }
  else if ( stype == "GAINS" ) 
  {
    data = new AliMUONTrackerData(Form("GAIN%d",startOfValidity),"Gains",7,kTRUE);
    data->SetDimensionName(0,"gain");
    data->SetDimensionName(1,"a1");
    data->SetDimensionName(2,"a2");
    data->SetDimensionName(3,"thres");
    data->SetDimensionName(4,"qual1");
    data->SetDimensionName(5,"qual2");
    data->SetDimensionName(6,"sat");    
  }
  else if ( stype == "HV" ) 
  {
    data = new AliMUONTrackerData(Form("HV%d",startOfValidity),"High Voltages",1); //,!isSingleEvent);
		data->SetDimensionName(0,"HV");
  }
  else if ( stype == "OCCUPANCY" ) 
  {
    data = new AliMUONTrackerData(Form("OCC%d",startOfValidity),"OccupancyMap",store);
    data->SetDimensionName(0,"One");
    return data; // important to return now to avoid the data->Add(store) later on...
  }
  else if ( stype == "PEDESTALS" ) 
  {
    data  = new AliMUONTrackerData(Form("PED%d",startOfValidity),"Pedestals",2,kTRUE);
    data->SetDimensionName(0,"Mean");
    data->SetDimensionName(1,"Sigma");    
  }
  else if ( stype == "STATUS" ) 
  {
    data = new AliMUONTrackerData(Form("STATUS%d",startOfValidity),"Status",1,kTRUE);
    data->SetDimensionName(0,"Bits");
  }
  else if ( stype == "STATUSMAP" ) 
  {
    data = new AliMUONTrackerData(Form("STATUSMAP%d",startOfValidity),"Status map",2,kTRUE);
    data->SetDimensionName(0,"Bits");
    data->SetDimensionName(1,"Dead");
  }

  if (!data)
  {
    AliErrorClass(Form("Could not create data for type=%s",type));
    return 0x0;
  }
  
  data->Add(store);
  
  return data;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONTrackerConditionDataMaker::CreateHVStore(TMap& m)
{
  /// Create a store from hv values
  
  AliMUONVStore* store = new AliMUON2DMap(kTRUE);
  
  TIter next(&m);
  TObjString* s;
  AliMpDCSNamer hvNamer("TRACKER");
  
  while ( ( s = static_cast<TObjString*>(next()) ) )
  {
    TString name(s->String());
    
    Int_t hvIndex = hvNamer.DCSIndexFromDCSAlias(name.Data());

    Int_t detElemId = hvNamer.DetElemIdFromDCSAlias(name.Data());
    
    if ( hvIndex >= 0 && detElemId < 0 )
    {
      // skip switches
      continue;      
    }
    
    if ( !AliMpDEManager::IsValidDetElemId(detElemId) )
    {
      AliErrorClass(Form("Got an invalid DE = %d from alias = %s",
                         detElemId,name.Data()));
      continue;
    }

    Int_t nPCBs = hvNamer.NumberOfPCBs(detElemId);
    Int_t indexMin = nPCBs ? 0 : hvIndex;
    Int_t indexMax = nPCBs ? nPCBs : hvIndex+1;
    
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    
    for ( int i = indexMin ; i < indexMax; ++i )
    {
      Float_t switchValue(1.0);
      
      if ( nPCBs ) 
      {
        TString switchName(hvNamer.DCSSwitchName(detElemId,i));

        TPair* p = static_cast<TPair*>(m.FindObject(switchName.Data()));
        TObjArray* a = static_cast<TObjArray*>(p->Value());
        
        switchValue = AliMUONPadStatusMaker::SwitchValue(*a);                                           
      }
      
      const AliMpArrayI* manus = de->ManusForHV(i);
      
      if (!manus) continue;
      
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
      hvValue *= switchValue;  
      
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
AliMUONVStore*
AliMUONTrackerConditionDataMaker::CreateStatusStore(Int_t runNumber)
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
AliMUONTrackerConditionDataMaker::CreateStatusMapStore(Int_t runNumber)
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
AliMUONTrackerConditionDataMaker::CreateStore(Int_t runNumber, 
                                              const char* source, 
                                              const char* type, 
                                              Int_t& startOfValidity)
{
  /// Create the store by reading it from OCDB or from an ASCII file
  
  TString stype(type);
  stype.ToUpper();
  
  AliMUONVStore* store(0x0);
  
  startOfValidity = 0;
  
  Bool_t ocdb = (runNumber>=0);
  
  if ( stype == "CAPACITANCES" )
  {    
    if ( ocdb ) 
    {
      store = AliMUONCalibrationData::CreateCapacitances(runNumber,&startOfValidity);    
    }
    else
    {
      store = new AliMUON2DMap(20000);
      AliMUONTrackerIO::DecodeCapacitances(source,*store);
    }
  }
  else if ( stype == "CONFIG" ) 
  {
    AliMUONVStore* tmp(0x0);
    if ( ocdb ) 
    {
      tmp = AliMUONCalibrationData::CreateConfig(runNumber,&startOfValidity);
    }
    else
    {
      tmp = new AliMUON2DMap(kTRUE);
      AliMUONTrackerIO::DecodeConfig(source,*tmp);
    }
    if ( tmp ) 
    {
      store = ExpandConfig(*tmp);      
    }
    delete tmp;
  }
  else if ( stype == "GAINS" ) 
  {
    AliMUONVStore* gains(0x0);
    if ( ocdb ) 
    {
      gains = AliMUONCalibrationData::CreateGains(runNumber,&startOfValidity);
    }
    else
    {
      gains = new AliMUON2DMap(kTRUE);
      TString comment;
      AliMUONTrackerIO::DecodeGains(source,*gains,comment);
    }
    store = PatchGainStore(*gains);
    delete gains;
  }
  else if ( stype == "OCCUPANCY" ) 
  {
    if ( ocdb ) 
    {
      store = AliMUONCalibrationData::CreateOccupancyMap(runNumber,&startOfValidity);
    }
    else
    {
      store = new AliMUON2DMap(kTRUE);
      AliMUONTrackerIO::DecodeOccupancy(source,*store);
    }
  }
  else if ( stype == "PEDESTALS" ) 
  {
    if ( ocdb ) 
    {
      store = AliMUONCalibrationData::CreatePedestals(runNumber,&startOfValidity);
    }
    else
    {
      store = new AliMUON2DMap(kTRUE);
      AliMUONTrackerIO::DecodePedestals(source,*store);
    }
  }
  
  /// Below are source that can only be accessed from OCDB
  if (!store && !ocdb) 
  {
    return 0x0;
  }
  
  if ( stype == "HV" ) 
  {
    TMap* m = AliMUONCalibrationData::CreateHV(runNumber,&startOfValidity);
		store = CreateHVStore(*m);
    delete m;
  }
  else if ( stype == "STATUS" ) 
  {
    store = CreateStatusStore(runNumber);
  }
  else if ( stype == "STATUSMAP" ) 
  {
    store = CreateStatusMapStore(runNumber);
  }
  
  return store;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONTrackerConditionDataMaker::ExpandConfig(const AliMUONVStore& manuConfig)
{
  /// Convert the config from manu level to channel level (just to
  /// be able to add it correctly to the trackerdata...)
  
  AliMUONVStore* store = manuConfig.Create();
  
  TIter next(manuConfig.CreateIterator());
  AliMUONVCalibParam* p;
  
  while ( ( p = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    AliMUONVCalibParam* c = new AliMUONCalibParamNF(1,AliMpConstants::ManuNofChannels(),p->ID0(),p->ID1(),0.0);
    
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(p->ID0());
    
    for ( Int_t i = 0; i < c->Size(); ++i ) 
    {
      if ( de->IsExistingChannel(p->ID1(),i) )
      {
        c->SetValueAsFloat(i,0,1.0);        
      }
    }
    
    store->Add(c);
  }
  return store;
}

//_____________________________________________________________________________
Long64_t 
AliMUONTrackerConditionDataMaker::Merge(TCollection*)
{
  /// Merge
  AliError("Not implemented. Does it have sense ?");
  return 0;
}

//_____________________________________________________________________________
AliMUONVStore*
AliMUONTrackerConditionDataMaker::PatchGainStore(const AliMUONVStore& gains)
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

