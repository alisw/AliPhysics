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
#include "AliMpArrayI.h"
#include "AliMpConstants.h"
#include "AliMpDCSNamer.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpManuIterator.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONCalibParamNF.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONPadStatusMapMaker.h"
#include "AliMUONRejectList.h"
#include "AliMUONTrackerData.h"
#include "AliMUONTrackerDataSourceTypes.h"
#include "AliMUONTrackerIO.h"
#include "Riostream.h"
#include "TClass.h"
#include "TMap.h"
#include "TObjString.h"
#include "TString.h"
#include "TSystem.h"
#include <sstream>

//_____________________________________________________________________________
AliMUONTrackerConditionDataMaker::AliMUONTrackerConditionDataMaker():
AliMUONVTrackerDataMaker(),
fData(0x0),
fSource(""),
fIsOwnerOfData(kTRUE)
{
  /// default ctor to be able to stream
}

//_____________________________________________________________________________
AliMUONTrackerConditionDataMaker::AliMUONTrackerConditionDataMaker(Int_t runNumber, const char* ocdbPath, const char* type):
AliMUONVTrackerDataMaker(),
fData(0x0),
fSource(Form("%s-%010d-%s",ocdbPath,runNumber,type)),
fIsOwnerOfData(kTRUE)
{
  /// ctor from OCDB

  AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();

	AliCDBManager::Instance()->SetDefaultStorage(ocdbPath);

  Int_t startOfValidity;

  if ( AliMUONTrackerDataSourceTypes::IsRejectList(type) )
  {
    AliMUONRejectList* rl = AliMUONCalibrationData::CreateRejectList(runNumber,&startOfValidity);

    if (rl)
    {
      fData = new AliMUONTrackerData(Form("RL%d",startOfValidity),"RejectList",*rl);
    }

    delete rl;
  }
  else
  {
    AliMUONVStore* store = CreateStore(runNumber,ocdbPath,type,startOfValidity);

    AliDebug(1,Form("runNumber=%d ocdbPath=%s type=%s startOfValidity=%d store=%p",
                    runNumber,ocdbPath,type,startOfValidity,store));
    if ( store )
    {
      fData = CreateData(type,*store,startOfValidity);
    }

    AliDebug(1,Form("runNumber=%d ocdbPath=%s type=%s startOfValidity=%d store=%p",
                    runNumber,ocdbPath,type,startOfValidity,store));

    delete store;
  }

  if ( fData )
  {
    TString shortName(fData->GetName());
    TString cdbPath(ocdbPath);

    shortName = type;

    shortName += Form("%d",startOfValidity);

    shortName += "(";

    if ( cdbPath.Contains("cvmfs/alice") )
    {
      shortName += "cvmfs";
    }
    else if ( cdbPath.BeginsWith("alien") && cdbPath.Contains("/alice/data") )
    {
      shortName += "alien";
    }
    else if ( cdbPath.BeginsWith("alien") && cdbPath.Contains("user") )
    {
      shortName.ReplaceAll("/alice.cern.ch/user/","...");
    }
    else
    {
      shortName += cdbPath;
    }

    shortName += ")";


    fData->SetName(shortName);
  }

  AliCDBManager::Instance()->SetDefaultStorage(storage);
}

//_____________________________________________________________________________
AliMUONTrackerConditionDataMaker::AliMUONTrackerConditionDataMaker(const char* filename, const char* type):
AliMUONVTrackerDataMaker(),
fData(0x0),
fSource(Form("%s-%s",filename,type)),
fIsOwnerOfData(kTRUE)
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
fSource(Form("direct-%s",type)),
fIsOwnerOfData(kTRUE)

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
  if ( fIsOwnerOfData ) delete fData;
}


//_____________________________________________________________________________
AliMUONVTrackerData*
AliMUONTrackerConditionDataMaker::CreateData(const char* type, AliMUONVStore& store, Int_t startOfValidity)
{
  /// Create the data source
  AliMUONVTrackerData* data(0x0);

  if ( AliMUONTrackerDataSourceTypes::IsConfig(type ) )
  {
    data = new AliMUONTrackerData(Form("%s%d",AliMUONTrackerDataSourceTypes::ShortNameForConfig(),startOfValidity),"Configuration",1);
    data->SetDimensionName(0,"there");
    data->DisableChannelLevel();
  }
  else if ( AliMUONTrackerDataSourceTypes::IsHV(type) )
  {
    data = new AliMUONTrackerData(Form("%s%d",AliMUONTrackerDataSourceTypes::ShortNameForHV(),startOfValidity),"High Voltages",1); //,!isSingleEvent);
		data->SetDimensionName(0,"HV");
  }
  else if ( AliMUONTrackerDataSourceTypes::IsLV(type) )
  {
    data = new AliMUONTrackerData(Form("%s%d",AliMUONTrackerDataSourceTypes::ShortNameForLV(),startOfValidity),"Low Voltages",3); //,!isSingleEvent);
		data->SetDimensionName(0,"ann"); // analog negative
    data->SetDimensionName(1,"dig"); // digital
    data->SetDimensionName(2,"anp"); // analog positive
  }
  else if ( AliMUONTrackerDataSourceTypes::IsOccupancy(type) )
  {
    data = new AliMUONTrackerData(Form("%s%d",AliMUONTrackerDataSourceTypes::ShortNameForOccupancy(),startOfValidity),"OccupancyMap",store);
    data->SetDimensionName(0,"One");
    return data; // important to return now to avoid the data->Add(store) later on...
  }
  else if ( AliMUONTrackerDataSourceTypes::IsPedestals(type) )
  {
    data  = new AliMUONTrackerData(Form("%s%d",AliMUONTrackerDataSourceTypes::ShortNameForPedestals(),startOfValidity),"Pedestals",2,kTRUE);
    data->SetDimensionName(0,"Mean");
    data->SetDimensionName(1,"Sigma");
  }
  else if ( AliMUONTrackerDataSourceTypes::IsStatus(type) )
  {
    data = new AliMUONTrackerData(Form("%s%d",AliMUONTrackerDataSourceTypes::ShortNameForStatus(),startOfValidity),"Status",1,kTRUE);
    data->SetDimensionName(0,"Bits");
  }
  else if ( AliMUONTrackerDataSourceTypes::IsStatusMap(type) )
  {
    data = new AliMUONTrackerData(Form("%s%d",AliMUONTrackerDataSourceTypes::ShortNameForStatusMap(),startOfValidity),"Status map",2,kTRUE);
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
        TString switchName(hvNamer.DCSSwitchAliasName(detElemId,i));

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
      Int_t noff(0);

      while ( ( v = static_cast<AliDCSValue*>(n2()) ) )
      {
        hvValue += v->GetFloat();
        if ( v->GetFloat() < AliMpDCSNamer::TrackerHVOFF() ) ++noff;
        ++n;
      }
      hvValue *= switchValue;

      if ( n ) hvValue /= n;

      if (noff>0 && noff<n)
      {
        // that's a trip
        hvValue = -1.0;
      }

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
AliMUONTrackerConditionDataMaker::CreateLVStore(TMap& m)
{
  /// Create a store from low voltage values

  AliMUONVStore* store = new AliMUON2DMap(kTRUE);

  TIter next(&m);
  TObjString* s;
  AliMpDCSNamer hvNamer("TRACKER");

  while ( ( s = static_cast<TObjString*>(next()) ) )
  {
    TString name(s->String());

    Int_t index(0);

    if ( name.Contains("ann") ) index = 0;
    if ( name.Contains("dig") ) index = 1;
    if ( name.Contains("anp") ) index = 2;

    Int_t* detElemId(0x0);
    Int_t numberOfDetectionElements;
    AliMp::PlaneType planeType;

    hvNamer.DecodeDCSMCHLVAlias(name.Data(), detElemId, numberOfDetectionElements, planeType);

    if (planeType != AliMp::kBendingPlane)
    {
      std::cout << "cool " << detElemId[0] << std::endl;
    }
    // compute the value to be associated with that alias : either the mean
    // or zero if a trip is detected
    Float_t lvValue(0.0);

    TPair* p = static_cast<TPair*>(m.FindObject(name.Data()));
    TObjArray* a = static_cast<TObjArray*>(p->Value());
    TIter n2(a);
    AliDCSValue* v;
    Int_t n(0);
    Int_t noff(0);

    while ( ( v = static_cast<AliDCSValue*>(n2()) ) )
    {
      lvValue += v->GetFloat();
      if ( v->GetFloat() < AliMpDCSNamer::TrackerLVOFF() ) ++noff;
      ++n;
    }

    if ( n ) lvValue /= n;

    if (noff>0 && noff<n)
    {
      // that's a trip
      lvValue = -1.0;
    }

    Int_t nofChannels(AliMpConstants::ManuNofChannels());

    // now assign this LV value to all the manus of the relevant plane(s)
    // of the detection element(s)

    for ( int i = 0; i < numberOfDetectionElements; ++i )
    {
      AliMpManuIterator manuIterator;
      Int_t de, manuId;

      //std::cout << Form("DE %04d planeType %d value %g",detElemId[i],planeType,lvValue) << std::endl;

      while (manuIterator.Next(de,manuId))
      {
        if ( de != detElemId[i] ) continue;
        if ( de < 500 )
        {
          if ( planeType == AliMp::kNonBendingPlane &&
             ( ( manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane) ) == 0 ) ) continue;
             if ( planeType == AliMp::kBendingPlane &&
                ( ( manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane) ) != 0 ) ) continue;
        }
        AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(store->FindObject(de,manuId));
        if ( ! param )
        {
          param = new AliMUONCalibParamND(3,nofChannels,de,manuId,0);
          store->Add(param);
        }
        for ( Int_t j = 0 ; j < nofChannels; ++j )
        {
          param->SetValueAsDouble(j,index,lvValue);
        }
      }

    }
    delete[] detElemId;
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

  AliMUONVStore* store(0x0);

  startOfValidity = 0;

  Bool_t ocdb = (runNumber>=0);

  if ( AliMUONTrackerDataSourceTypes::IsConfig(type) )
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
  else if ( AliMUONTrackerDataSourceTypes::IsOccupancy(type) )
  {
    if ( ocdb )
    {
      store = AliMUONCalibrationData::CreateOccupancyMap(runNumber,&startOfValidity);
      if (store) store = static_cast<AliMUONVStore*>(store->Clone());
    }
    else
    {
      store = new AliMUON2DMap(kTRUE);
      AliMUONTrackerIO::DecodeOccupancy(source,*store);
    }
  }
  else if ( AliMUONTrackerDataSourceTypes::IsPedestals(type) )
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

  if ( AliMUONTrackerDataSourceTypes::IsHV(type) )
  {
    TMap* m = AliMUONCalibrationData::CreateHV(runNumber,&startOfValidity);
		store = CreateHVStore(*m);
    delete m;
  }
  if ( AliMUONTrackerDataSourceTypes::IsLV(type) )
  {
    TMap* m = AliMUONCalibrationData::CreateLV(runNumber,&startOfValidity);
		store = CreateLVStore(*m);
    delete m;
  }
  else if ( AliMUONTrackerDataSourceTypes::IsStatus(type) )
  {
    store = CreateStatusStore(runNumber);
  }
  else if ( AliMUONTrackerDataSourceTypes::IsStatusMap(type) )
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
