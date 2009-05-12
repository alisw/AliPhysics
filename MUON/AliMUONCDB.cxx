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

/* $Id$ */

//-----------------------------------------------------------------------------
/// \class AliMUONCDB
///
/// Helper class to experience the OCDB
/// It allows to generate dummy (but complete) containers for all the
/// calibration data types we have for tracker and trigger, and to write
/// them into OCDB.
///
/// For more information, please see READMECDB
///
// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMUONCDB.h"

#include "AliMUON1DArray.h"
#include "AliMUON1DMap.h"
#include "AliMUON2DMap.h"
#include "AliMUON2DStoreValidator.h"
#include "AliMUONCalibParamNF.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONConstants.h"
#include "AliMUONTrackerIO.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerLut.h"
#include "AliMUONVStore.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONGlobalCrateConfig.h"
#include "AliMUONRegionalTriggerConfig.h"

#include "AliMpCDB.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpManuStore.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpFiles.h"
#include "AliMpDCSNamer.h"
#include "AliMpManuIterator.h"
#include "AliMpSegmentation.h"
#include "AliMpStationType.h"
#include "AliMpVSegmentation.h"

#include "AliCodeTimer.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliDCSValue.h"
#include "AliLog.h"

#include <Riostream.h>
#include <TArrayI.h>
#include <TClass.h>
#include <TH1F.h>
#include <TList.h>
#include <TMap.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TMath.h>


/// \cond CLASSIMP
ClassImp(AliMUONCDB)
/// \endcond

namespace
{
  //_____________________________________________________________________________
AliMUONVStore* Create2DMap()
{
  return new AliMUON2DMap(true);
}

  //_____________________________________________________________________________
void getBoundaries(const AliMUONVStore& store, Int_t dim,
                   Float_t* xmin, Float_t* xmax)
{
  /// Assuming the store contains AliMUONVCalibParam objects, compute the
  /// limits of the value contained in the VCalibParam, for each of its dimensions
  /// xmin and xmax must be of dimension dim
  
  for ( Int_t i = 0; i < dim; ++i ) 
  {
    xmin[i]=1E30;
    xmax[i]=-1E30;
  }
  
  TIter next(store.CreateIterator());
  AliMUONVCalibParam* value;
  
  while ( ( value = dynamic_cast<AliMUONVCalibParam*>(next() ) ) )
  {
    Int_t detElemId = value->ID0();
    Int_t manuId = value->ID1();
    
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
        
    for ( Int_t manuChannel = 0; manuChannel < value->Size(); ++manuChannel )
    {
      AliMpPad pad = seg->PadByLocation(manuId,manuChannel,kFALSE);
      if (!pad.IsValid()) continue;
      
      for ( Int_t i = 0; i < dim; ++i ) 
      {
        Float_t x0 = value->ValueAsFloat(manuChannel,i);
      
        xmin[i] = TMath::Min(xmin[i],x0);
        xmax[i] = TMath::Max(xmax[i],x0);
      }
    }
  }

  for ( Int_t i = 0; i < dim; ++i ) 
  {
    if ( TMath::Abs(xmin[i]-xmax[i]) < 1E-3 ) 
    {
      xmin[i] -= 1;
      xmax[i] += 1;
    }
  }
}

//_____________________________________________________________________________
Double_t GetRandom(Double_t mean, Double_t sigma, Bool_t mustBePositive)
{
  Double_t x(-1);
  if ( mustBePositive ) 
  {
    while ( x < 0 ) 
    {
      x = gRandom->Gaus(mean,sigma);
    }
  }
  else
  {
    x = gRandom->Gaus(mean,sigma);
  }
  return x;
}

}

//_____________________________________________________________________________
AliMUONCDB::AliMUONCDB(const char* cdbpath)
: TObject(),
  fCDBPath(cdbpath),
  fMaxNofChannelsToGenerate(-1)
{
  /// ctor
    // Load mapping
    if ( ! AliMpCDB::LoadDDLStore() ) {
      AliFatal("Could not access mapping from OCDB !");
    }

    if ( ! AliMpCDB::LoadManuStore() ) {
      AliFatal("Could not access run-dependent mapping from OCDB !");
    }
}

//_____________________________________________________________________________
AliMUONCDB::~AliMUONCDB()
{
  /// dtor
}

//_____________________________________________________________________________
AliMUONVStore* 
AliMUONCDB::Diff(AliMUONVStore& store1, AliMUONVStore& store2, 
                 const char* opt)
{
  /// creates a store which contains store1-store2
  /// if opt="abs" the difference is absolute one,
  /// if opt="rel" then what is stored is (store1-store2)/store1
  /// if opt="percent" then what is stored is rel*100
  ///
  /// WARNING Works only for stores which holds AliMUONVCalibParam objects
  
  TString sopt(opt);
  sopt.ToUpper();
  
  if ( !sopt.Contains("ABS") && !sopt.Contains("REL") && !sopt.Contains("PERCENT") )
  {
    AliErrorClass(Form("opt %s not supported. Only ABS, REL, PERCENT are",opt));
    return 0x0;
  }
  
  AliMUONVStore* d = static_cast<AliMUONVStore*>(store1.Clone());
  
  TIter next(d->CreateIterator());
  
  AliMUONVCalibParam* param;
  
  while ( ( param = dynamic_cast<AliMUONVCalibParam*>(next() ) ) )
  {
    Int_t detElemId = param->ID0();
    Int_t manuId = param->ID1();
    
    AliMUONVCalibParam* param2 = dynamic_cast<AliMUONVCalibParam*>(store2.FindObject(detElemId,manuId));
    //FIXME: this might happen. Handle it.
    if (!param2) 
    {
      cerr << "param2 is null : FIXME : this might happen !" << endl;
      delete d;
      return 0;
    }
    
    for ( Int_t i = 0; i < param->Size(); ++i )
    {
      for ( Int_t j = 0; j < param->Dimension(); ++j )
      {
        Float_t value(0);
        if ( sopt.Contains("ABS") )
        {
          value = param->ValueAsFloat(i,j) - param2->ValueAsFloat(i,j);
        }
        else if ( sopt.Contains("REL") || sopt.Contains("PERCENT") )
        {
          if ( param->ValueAsFloat(i,j) ) 
          {
            value = (param->ValueAsFloat(i,j) - param2->ValueAsFloat(i,j))/param->ValueAsFloat(i,j);
          }
          else 
          {
            continue;
          }
          if ( sopt.Contains("PERCENT") ) value *= 100.0;
        }
        param->SetValueAsFloat(i,j,value);
      }      
    }
  }
  return d;
}

//_____________________________________________________________________________
void 
AliMUONCDB::Plot(const AliMUONVStore& store, const char* name, Int_t nbins)
{
  /// Make histograms of each dimension of the AliMUONVCalibParam
  /// contained inside store.
  /// It produces histograms named name_0, name_1, etc...
  
  TIter next(store.CreateIterator());
  AliMUONVCalibParam* param;
  Int_t n(0);
  const Int_t kNStations = AliMpConstants::NofTrackingChambers()/2;
  Int_t* nPerStation = new Int_t[kNStations];
  TH1** h(0x0);
  
  for ( Int_t i = 0; i < kNStations; ++i ) nPerStation[i]=0;
  
  while ( ( param = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    if (!h)
    {
      Int_t dim = param->Dimension();
      h = new TH1*[dim];
      Float_t* xmin = new Float_t[dim];
      Float_t* xmax = new Float_t[dim];
      getBoundaries(store,dim,xmin,xmax);
      
      for ( Int_t i = 0; i < dim; ++i ) 
      {
        h[i] = new TH1F(Form("%s_%d",name,i),Form("%s_%d",name,i),
                            nbins,xmin[i],xmax[i]);
        AliInfo(Form("Created histogram %s",h[i]->GetName()));
      }
    }
    
    Int_t detElemId = param->ID0();
    Int_t manuId = param->ID1();
    Int_t station = AliMpDEManager::GetChamberId(detElemId)/2;
    
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
    
    for ( Int_t manuChannel = 0; manuChannel < param->Size(); ++manuChannel )
    {
      AliMpPad pad = seg->PadByLocation(manuId,manuChannel,kFALSE);
      if (!pad.IsValid()) continue;

      ++n;
      ++nPerStation[station];
      
      for ( Int_t dim = 0; dim < param->Dimension(); ++dim ) 
      {
        h[dim]->Fill(param->ValueAsFloat(manuChannel,dim));
      }
    }
  } 
  
  for ( Int_t i = 0; i < kNStations; ++i )
  {
    AliInfo(Form("Station %d %d ",(i+1),nPerStation[i]));
  }

  AliInfo(Form("Number of channels = %d",n));
  
  delete[] nPerStation;
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeHVStore(TMap& aliasMap, Bool_t defaultValues)
{
  /// Create a HV store
  
  AliMpDCSNamer hvNamer("TRACKER");
  
  TObjArray* aliases = hvNamer.GenerateAliases();
  
  Int_t nSwitch(0);
  Int_t nChannels(0);
  
  for ( Int_t i = 0; i < aliases->GetEntries(); ++i ) 
  {
    TObjString* alias = static_cast<TObjString*>(aliases->At(i));
    TString& aliasName = alias->String();
    if ( aliasName.Contains("sw") ) 
    {
      // HV Switch (St345 only)
      TObjArray* valueSet = new TObjArray;
      valueSet->SetOwner(kTRUE);
      
      Bool_t value = kTRUE;
      
      if (!defaultValues)
      {
        Float_t r = gRandom->Uniform();
        if ( r < 0.007 ) value = kFALSE;      
      } 
      
      for ( UInt_t timeStamp = 0; timeStamp < 60*3; timeStamp += 60 )
      {
        AliDCSValue* dcsValue = new AliDCSValue(value,timeStamp);
        valueSet->Add(dcsValue);
      }
      aliasMap.Add(new TObjString(*alias),valueSet);
      ++nSwitch;
    }
    else
    {
      TObjArray* valueSet = new TObjArray;
      valueSet->SetOwner(kTRUE);
      for ( UInt_t timeStamp = 0; timeStamp < 60*15; timeStamp += 120 )
      {
        Float_t value = 1500;
        if (!defaultValues) value = GetRandom(1750,62.5,true);
        AliDCSValue* dcsValue = new AliDCSValue(value,timeStamp);
        valueSet->Add(dcsValue);
      }
      aliasMap.Add(new TObjString(*alias),valueSet);
      ++nChannels;
    }
  }
  
  delete aliases;
  
  AliInfo(Form("%d HV channels and %d switches",nChannels,nSwitch));
  
  return nChannels+nSwitch;
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeTriggerDCSStore(TMap& aliasMap, Bool_t defaultValues)
{
  /// Create a Trigger HV and Currents store
  
  AliMpDCSNamer triggerDCSNamer("TRIGGER");
  
  TObjArray* aliases = triggerDCSNamer.GenerateAliases();
  
  Int_t nChannels[2] = {0, 0};
  
  for ( Int_t i = 0; i < aliases->GetEntries(); ++i ) 
  {
    TObjString* alias = static_cast<TObjString*>(aliases->At(i));
    TString& aliasName = alias->String();

    TObjArray* valueSet = new TObjArray;
    valueSet->SetOwner(kTRUE);
    Int_t measureType = triggerDCSNamer.DCSvariableFromDCSAlias(aliasName.Data());
    for ( UInt_t timeStamp = 0; timeStamp < 60*15; timeStamp += 120 )
    {
      Float_t value = 
	(measureType == AliMpDCSNamer::kDCSI) ? 2. : 8000.;
      if (!defaultValues) {
	switch (measureType){
	case AliMpDCSNamer::kDCSI:
	  value = GetRandom(2.,0.4,true);
	  break;
	case AliMpDCSNamer::kDCSHV:
	  value = GetRandom(8000.,16.,true);
	  break;
	}
      }
      AliDCSValue* dcsValue = new AliDCSValue(value,timeStamp);
      valueSet->Add(dcsValue);
    }
    aliasMap.Add(new TObjString(*alias),valueSet);
    ++nChannels[measureType];
  }
  
  delete aliases;
  
  AliInfo(Form("Trigger channels I -> %i   HV -> %i",nChannels[0], nChannels[1]));
  
  return nChannels[0] + nChannels[1];
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakePedestalStore(AliMUONVStore& pedestalStore, Bool_t defaultValues)
{
  /// Create a pedestal store. if defaultValues=true, ped.mean=ped.sigma=1,
  /// otherwise mean and sigma are from a gaussian (with parameters
  /// defined below by the kPedestal* constants)

  AliCodeTimerAuto("");
  
  Int_t nchannels(0);
  Int_t nmanus(0);
  
  const Int_t kChannels(AliMpConstants::ManuNofChannels());
  
  // bending
  const Float_t kPedestalMeanMeanB(200.);
  const Float_t kPedestalMeanSigmaB(10.);
  const Float_t kPedestalSigmaMeanB(1.);
  const Float_t kPedestalSigmaSigmaB(0.2);
  
  // non bending
  const Float_t kPedestalMeanMeanNB(200.);
  const Float_t kPedestalMeanSigmaNB(10.);
  const Float_t kPedestalSigmaMeanNB(1.);
  const Float_t kPedestalSigmaSigmaNB(0.2);
  
  const Float_t kFractionOfDeadManu(0.); // within [0.,1.]

  Int_t detElemId;
  Int_t manuId;
    
  AliMpManuIterator it;
  
  while ( it.Next(detElemId,manuId) )
  {
    // skip a given fraction of manus
    if (kFractionOfDeadManu > 0. && gRandom->Uniform() < kFractionOfDeadManu) continue;
    
    ++nmanus;

    AliMUONVCalibParam* ped = 
      new AliMUONCalibParamNF(2,kChannels,detElemId,manuId,AliMUONVCalibParam::InvalidFloatValue());

    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    
    for ( Int_t manuChannel = 0; manuChannel < kChannels; ++manuChannel )
    {
      if ( ! de->IsConnectedChannel(manuId,manuChannel) ) continue;
      
      ++nchannels;
      
      Float_t meanPedestal;
      Float_t sigmaPedestal;
      
      if ( defaultValues ) 
      {
        meanPedestal = 0.0;
        sigmaPedestal = 1.0;
      }
      else
      {
        Bool_t positive(kTRUE);
        meanPedestal = 0.0;
	
	if ( manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane) ) { // manu in non bending plane
	  
	  while ( meanPedestal == 0.0 ) // avoid strict zero 
	  {
	    meanPedestal = GetRandom(kPedestalMeanMeanNB,kPedestalMeanSigmaNB,positive);
	  }
	  sigmaPedestal = GetRandom(kPedestalSigmaMeanNB,kPedestalSigmaSigmaNB,positive);
	  
	} else { // manu in bending plane
	  
	  while ( meanPedestal == 0.0 ) // avoid strict zero 
	  {
	    meanPedestal = GetRandom(kPedestalMeanMeanB,kPedestalMeanSigmaB,positive);
	  }
	  sigmaPedestal = GetRandom(kPedestalSigmaMeanB,kPedestalSigmaSigmaB,positive);
	  
	}
	
      }
      
      ped->SetValueAsFloat(manuChannel,0,meanPedestal);
      ped->SetValueAsFloat(manuChannel,1,sigmaPedestal);
      
    }
    Bool_t ok = pedestalStore.Add(ped);
    if (!ok)
    {
      AliError(Form("Could not set DetElemId=%d manuId=%d",detElemId,manuId));
    }
    if ( fMaxNofChannelsToGenerate > 0 && nchannels >= fMaxNofChannelsToGenerate ) break;
  }
  
  AliInfo(Form("%d Manus and %d channels.",nmanus,nchannels));
  return nchannels;
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeKillMapStore(AliMUONVStore& killMapStore)
{
  /// Create a kill map.
  
  AliCodeTimerAuto("");
  
  Int_t nchannels(0);
  Int_t nmanus(0);
  
  const Int_t kChannels(AliMpConstants::ManuNofChannels());
  
  const Float_t kFractionOfDeadManu(0.1); // within [0.,1.]
  
  if ( kFractionOfDeadManu == 0.0 ) return 0.0;
  
  Int_t detElemId;
  Int_t manuId;
  
  AliMpManuIterator it;
  
  while ( it.Next(detElemId,manuId) )
  {
    // skip a given fraction of manus
    if ( gRandom->Uniform() > kFractionOfDeadManu) continue;
    
    ++nmanus;
    
    AliMUONVCalibParam* kill = new AliMUONCalibParamNI(1,kChannels,detElemId,manuId,0);
    
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    
    for ( Int_t manuChannel = 0; manuChannel < kChannels; ++manuChannel )
    {
      if ( ! de->IsConnectedChannel(manuId,manuChannel) ) continue;
      
      ++nchannels;
      
      kill->SetValueAsInt(manuChannel,0,1);
      
    }
    Bool_t ok = killMapStore.Add(kill);
    if (!ok)
    {
      AliError(Form("Could not set DetElemId=%d manuId=%d",detElemId,manuId));
    }
    if ( fMaxNofChannelsToGenerate > 0 && nchannels >= fMaxNofChannelsToGenerate ) break;
  }
  
  AliInfo(Form("%d Manus and %d channels.",nmanus,nchannels));
  return nchannels;
}

//_____________________________________________________________________________
Int_t
AliMUONCDB::MakeCapacitanceStore(AliMUONVStore& capaStore, const char* file)
{
  /// Read the capacitance values from file and append them to the capaStore
  
  return AliMUONTrackerIO::ReadCapacitances(file,capaStore);
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeCapacitanceStore(AliMUONVStore& capaStore, Bool_t defaultValues)
{
  /// Create a capacitance store. if defaultValues=true, all capa are 1.0,
  /// otherwise they are from a gaussian with parameters defined in the
  /// kCapa* constants below.

  AliCodeTimerAuto("");
  
  Int_t nchannels(0);
  Int_t nmanus(0);
  Int_t nmanusOK(0); // manus for which we got the serial number
    
  const Float_t kCapaMean(0.3);
  const Float_t kCapaSigma(0.1);
  const Float_t kInjectionGainMean(3);
  const Float_t kInjectionGainSigma(1);

  Int_t detElemId;
  Int_t manuId;
  
  AliMpManuIterator it;
  
  while ( it.Next(detElemId,manuId) )
  {
    ++nmanus;
    
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId); 
    Int_t serialNumber = AliMpManuStore::Instance()->GetManuSerial(detElemId, manuId);
      
    if ( serialNumber <= 0 ) continue;
    
    ++nmanusOK;
    
    AliMUONVCalibParam* capa = static_cast<AliMUONVCalibParam*>(capaStore.FindObject(serialNumber));
    
    if (!capa)
    {
      capa = new AliMUONCalibParamNF(2,AliMpConstants::ManuNofChannels(),serialNumber,0,1.0);
      Bool_t ok = capaStore.Add(capa);
      if (!ok)
      {
        AliError(Form("Could not set serialNumber=%d manuId=%d",serialNumber,manuId));
      }      
    }
    
    for ( Int_t manuChannel = 0; manuChannel < capa->Size(); ++manuChannel )
    {
      if ( ! de->IsConnectedChannel(manuId,manuChannel) ) continue;
      
      ++nchannels;
      
      Float_t capaValue;
      Float_t injectionGain;
      
      if ( defaultValues ) 
      {
        capaValue = 1.0;
        injectionGain = 1.0;
      }
      else
      {
        capaValue = GetRandom(kCapaMean,kCapaSigma,kTRUE);
        injectionGain = GetRandom(kInjectionGainMean,kInjectionGainSigma,kTRUE);
      }
      capa->SetValueAsFloat(manuChannel,0,capaValue);
      capa->SetValueAsFloat(manuChannel,1,injectionGain);
    }
  }
  
  Float_t percent = 0;
  if ( nmanus ) percent = 100*nmanusOK/nmanus;
  AliInfo(Form("%5d manus with serial number (out of %5d manus = %3.0f%%)",
               nmanusOK,nmanus,percent));
  AliInfo(Form("%5d channels",nchannels));
  if ( percent < 100 ) 
  {
    AliWarning("Did not get all serial numbers. capaStore is incomplete !!!!");
  }
  return nchannels;
  
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeGainStore(AliMUONVStore& gainStore, Bool_t defaultValues)
{  
  /// Create a gain store. if defaultValues=true, all gains set so that
  /// charge = (adc-ped)
  ///
  /// otherwise parameters are taken from gaussians with parameters 
  /// defined in the k* constants below.
  
  AliCodeTimerAuto("");
  
  Int_t nchannels(0);
  Int_t nmanus(0);
    
  const Int_t kSaturation(3000);
  const Double_t kA0Mean(1.2);
  const Double_t kA0Sigma(0.1);
  const Double_t kA1Mean(1E-5);
  const Double_t kA1Sigma(1E-6);
  const Double_t kQualMean(0xFF);
  const Double_t kQualSigma(0x10);
  const Int_t kThresMean(1600);
  const Int_t kThresSigma(100);
  
  Int_t detElemId;
  Int_t manuId;
  
  AliMpManuIterator it;
  
  while ( it.Next(detElemId,manuId) )
  {
    ++nmanus;

    AliMUONVCalibParam* gain = 
      new AliMUONCalibParamNF(5,AliMpConstants::ManuNofChannels(),
                              detElemId,
                              manuId,
                              AliMUONVCalibParam::InvalidFloatValue());

    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

    for ( Int_t manuChannel = 0; manuChannel < gain->Size(); ++manuChannel )
    {
      if ( ! de->IsConnectedChannel(manuId,manuChannel) ) continue;
      
      ++nchannels;
      
      if ( defaultValues ) 
      {
        gain->SetValueAsFloat(manuChannel,0,1.0);
        gain->SetValueAsFloat(manuChannel,1,0.0);
        gain->SetValueAsInt(manuChannel,2,4095); 
        gain->SetValueAsInt(manuChannel,3,1);
        gain->SetValueAsInt(manuChannel,4,kSaturation);
      }
      else
      {
        Bool_t positive(kTRUE);
        gain->SetValueAsFloat(manuChannel,0,GetRandom(kA0Mean,kA0Sigma,positive));
        gain->SetValueAsFloat(manuChannel,1,GetRandom(kA1Mean,kA1Sigma,!positive));
        gain->SetValueAsInt(manuChannel,2,(Int_t)TMath::Nint(GetRandom(kThresMean,kThresSigma,positive)));
        gain->SetValueAsInt(manuChannel,3,(Int_t)TMath::Nint(GetRandom(kQualMean,kQualSigma,positive)));
        gain->SetValueAsInt(manuChannel,4,kSaturation);        
      }
      
    }
    Bool_t ok = gainStore.Add(gain);
    if (!ok)
    {
      AliError(Form("Could not set DetElemId=%d manuId=%d",detElemId,manuId));
    }
    if ( fMaxNofChannelsToGenerate > 0 && nchannels >= fMaxNofChannelsToGenerate ) break;
  }
  
  AliInfo(Form("%d Manus and %d channels.",nmanus,nchannels));
  return nchannels;
}

//_____________________________________________________________________________
Int_t
AliMUONCDB::MakeLocalTriggerMaskStore(AliMUONVStore& localBoardMasks) const
{
  /// Generate local trigger masks store. All masks are set to FFFF
  
  AliCodeTimerAuto("");
  
  Int_t ngenerated(0);
  // Generate fake mask values for all localboards and put that into
  // one single container (localBoardMasks)
  for ( Int_t i = 1; i <= AliMpConstants::TotalNofLocalBoards(); ++i )
  {
    AliMUONVCalibParam* localBoard = new AliMUONCalibParamNI(1,8,i,0,0);
    for ( Int_t x = 0; x < 2; ++x )
    {
      for ( Int_t y = 0; y < 4; ++y )
      {
        Int_t index = x*4+y;
        localBoard->SetValueAsInt(index,0,0xFFFF);
        ++ngenerated;
      }
    }
    localBoardMasks.Add(localBoard);
  }
  return ngenerated;
}

//_____________________________________________________________________________
Int_t
AliMUONCDB::MakeRegionalTriggerConfigStore(AliMUONRegionalTriggerConfig& rtm) const
{
  /// Make a regional trigger config store. Mask is set to FFFF for each local board (Ch.F.)
  
  AliCodeTimerAuto("");
  
  if ( ! rtm.ReadData(AliMpFiles::LocalTriggerBoardMapping()) ) {
    AliErrorStream() << "Error when reading from mapping file" << endl;
    return 0;
  }
    
  return rtm.GetNofTriggerCrates();  
}


//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeGlobalTriggerConfigStore(AliMUONGlobalCrateConfig& gtm) const
{
  /// Make a global trigger config store. All masks (disable) set to 0x00 for each Darc board (Ch.F.)
  
  AliCodeTimerAuto("");
  
  return gtm.ReadData(AliMpFiles::GlobalTriggerBoardMapping());
}


//_____________________________________________________________________________
AliMUONTriggerLut* 
AliMUONCDB::MakeTriggerLUT(const char* file) const
{
  /// Make a triggerlut object, from a file.
  
  AliCodeTimerAuto("");
  
  AliMUONTriggerLut* lut = new AliMUONTriggerLut;
  lut->ReadFromFile(file);
  return lut;
}

//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells*
AliMUONCDB::MakeTriggerEfficiency(const char* file) const
{
  /// Make a trigger efficiency object from a file.
  
  AliCodeTimerAuto("");
  
  return new AliMUONTriggerEfficiencyCells(file);
}

//_____________________________________________________________________________
void 
AliMUONCDB::WriteToCDB(const char* calibpath, TObject* object, 
                       Int_t startRun, Int_t endRun, 
                       const char* filename)
{
  /// Write a given object to OCDB
  
  TString comment(gSystem->ExpandPathName(filename));
  
  WriteToCDB(object, calibpath, startRun, endRun, comment.Data());
}

//_____________________________________________________________________________
void 
AliMUONCDB::WriteToCDB(const char* calibpath, TObject* object, 
                       Int_t startRun, Int_t endRun, Bool_t defaultValues)
{
  /// Write a given object to OCDB
  
  TString comment;
  if ( defaultValues ) comment += "Test with default values";
  else comment += "Test with random values";
  
  WriteToCDB(object, calibpath, startRun, endRun, comment.Data());
}

//_____________________________________________________________________________
void
AliMUONCDB::WriteToCDB(TObject* object, const char* calibpath, Int_t startRun, Int_t endRun,
		       const char* comment, const char* responsible)
{
  /// Write a given object to OCDB
  
  AliCDBId id(calibpath,startRun,endRun);
  AliCDBMetaData md;
  md.SetAliRootVersion(gROOT->GetVersion());
  md.SetComment(comment);
  md.SetResponsible(responsible);
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) man->SetDefaultStorage(fCDBPath);
  man->Put(object,id,&md);
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeNeighbourStore(AliMUONVStore& neighbourStore)
{
  /// Fill the neighbours store with, for each channel, a TObjArray of its
  /// neighbouring pads (including itself)
  
  AliCodeTimerAuto("");
  
  AliInfo("Generating NeighbourStore. This will take a while. Please be patient.");
  
  Int_t nchannels(0);
  
  TObjArray tmp;

  Int_t detElemId;
  Int_t manuId;
  
  AliMpManuIterator it;
  
  while ( it.Next(detElemId,manuId) )
  {
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
    
    AliMUONVCalibParam* calibParam = static_cast<AliMUONVCalibParam*>(neighbourStore.FindObject(detElemId,manuId));
    if (!calibParam)
    {
      Int_t dimension(11);
      Int_t size(AliMpConstants::ManuNofChannels());
      Int_t defaultValue(-1);
      Int_t packingFactor(size);
      
      calibParam = new AliMUONCalibParamNI(dimension,size,detElemId,manuId,defaultValue,packingFactor);
      Bool_t ok = neighbourStore.Add(calibParam);
      if (!ok)
      {
        AliError(Form("Could not set DetElemId=%d manuId=%d",detElemId,manuId));
        return -1;
      }      
    }
    
    for ( Int_t manuChannel = 0; manuChannel < AliMpConstants::ManuNofChannels(); ++manuChannel )
    {
      AliMpPad pad = seg->PadByLocation(manuId,manuChannel,kFALSE);
      
      if (pad.IsValid()) 
      {
        ++nchannels;

        seg->GetNeighbours(pad,tmp,true,true);
        Int_t nofPadNeighbours = tmp.GetEntriesFast();
            
        for ( Int_t i = 0; i < nofPadNeighbours; ++i )
        {
          AliMpPad* p = static_cast<AliMpPad*>(tmp.UncheckedAt(i));
          Int_t x;
//          Bool_t ok =
              calibParam->PackValues(p->GetManuId(),p->GetManuChannel(),x);
//          if (!ok)
//          {
//            AliError("Could not pack value. Something is seriously wrong. Please check");
//            StdoutToAliError(pad->Print(););
//            return -1;
//          }
          calibParam->SetValueAsInt(manuChannel,i,x);
        }
      }
    }
    }
  
  return nchannels;
}

//_____________________________________________________________________________
void
AliMUONCDB::SetMaxNofChannelsToGenerate(Int_t n)
{
  /// Set the maximum number of channels to generate (used for testing only)
  /// n < 0 means no limit
  fMaxNofChannelsToGenerate = n;
}

//_____________________________________________________________________________
void
AliMUONCDB::WriteLocalTriggerMasks(Int_t startRun, Int_t endRun)
{  
  /// Write local trigger masks to OCDB
  
  AliMUONVStore* ltm = new AliMUON1DArray(AliMpConstants::TotalNofLocalBoards()+1);
  Int_t ngenerated = MakeLocalTriggerMaskStore(*ltm);
  AliInfo(Form("Ngenerated = %d",ngenerated));
  if (ngenerated>0)
  {
    WriteToCDB("MUON/Calib/LocalTriggerBoardMasks",ltm,startRun,endRun,true);
  }
  delete ltm;
}

//_____________________________________________________________________________
void
AliMUONCDB::WriteRegionalTriggerConfig(Int_t startRun, Int_t endRun)
{  
  /// Write regional trigger masks to OCDB
  
  AliMUONRegionalTriggerConfig* rtm = new AliMUONRegionalTriggerConfig();
  Int_t ngenerated = MakeRegionalTriggerConfigStore(*rtm);
  AliInfo(Form("Ngenerated = %d",ngenerated));
  if (ngenerated>0)
  {
    WriteToCDB("MUON/Calib/RegionalTriggerConfig",rtm,startRun,endRun,true);
  }
  delete rtm;
}


//_____________________________________________________________________________
void
AliMUONCDB::WriteGlobalTriggerConfig(Int_t startRun, Int_t endRun)
{  
  /// Write global trigger masks to OCDB
  
  AliMUONGlobalCrateConfig* gtm = new AliMUONGlobalCrateConfig();

  Int_t ngenerated = MakeGlobalTriggerConfigStore(*gtm);
  AliInfo(Form("Ngenerated = %d",ngenerated));
  if (ngenerated>0)
  {
    WriteToCDB("MUON/Calib/GlobalTriggerCrateConfig",gtm,startRun,endRun,true);
  }
  delete gtm;
}


//_____________________________________________________________________________
void
AliMUONCDB::WriteTriggerLut(Int_t startRun, Int_t endRun)
{  
  /// Write trigger LUT to OCDB
  
  AliMUONTriggerLut* lut = MakeTriggerLUT();
  if (lut)
  {
    WriteToCDB("MUON/Calib/TriggerLut",lut,startRun,endRun,true);
  }
  delete lut;
}

//_____________________________________________________________________________
void
AliMUONCDB::WriteTriggerEfficiency(Int_t startRun, Int_t endRun)
{  
  /// Write trigger efficiency to OCDB
  
  AliMUONTriggerEfficiencyCells* eff = MakeTriggerEfficiency();
  if (eff)
  {
    WriteToCDB("MUON/Calib/TriggerEfficiency",eff,startRun,endRun,true);
  }
  delete eff;
}

//_____________________________________________________________________________
void 
AliMUONCDB::WriteNeighbours(Int_t startRun, Int_t endRun)
{
  /// Write neighbours to OCDB
  
  AliMUONVStore* neighbours = Create2DMap();
  Int_t ngenerated = MakeNeighbourStore(*neighbours);
  AliInfo(Form("Ngenerated = %d",ngenerated));
  if (ngenerated>0)
  {
    WriteToCDB("MUON/Calib/Neighbours",neighbours,startRun,endRun,true);
  }
  delete neighbours;
}

//_____________________________________________________________________________
void 
AliMUONCDB::WriteHV(Bool_t defaultValues,
                    Int_t startRun, Int_t endRun)
{
  /// generate HV values (either cste = 1500 V) if defaultValues=true or random
  /// if defaultValues=false, see makeHVStore) and
  /// store them into CDB located at cdbpath, with a validity period
  /// ranging from startRun to endRun
  
  TMap* hvStore = new TMap;
  Int_t ngenerated = MakeHVStore(*hvStore,defaultValues);
  AliInfo(Form("Ngenerated = %d",ngenerated));
  if (ngenerated>0)
  {
    WriteToCDB("MUON/Calib/HV",hvStore,startRun,endRun,defaultValues);
  }
  delete hvStore;
}

//_____________________________________________________________________________
void 
AliMUONCDB::WriteTriggerDCS(Bool_t defaultValues,
                    Int_t startRun, Int_t endRun)
{
  /// generate Trigger HV and current values (either const if defaultValues=true or random
  /// if defaultValues=false, see makeTriggerDCSStore) and
  /// store them into CDB located at cdbpath, with a validity period
  /// ranging from startRun to endRun
  
  TMap* triggerDCSStore = new TMap;
  Int_t ngenerated = MakeTriggerDCSStore(*triggerDCSStore,defaultValues);
  AliInfo(Form("Ngenerated = %d",ngenerated));
  if (ngenerated>0)
  {
    WriteToCDB("MUON/Calib/TriggerDCS",triggerDCSStore,startRun,endRun,defaultValues);
  }
  delete triggerDCSStore;
}

//_____________________________________________________________________________
void 
AliMUONCDB::WritePedestals(Bool_t defaultValues,
                           Int_t startRun, Int_t endRun)
{
  /// generate pedestal values (either 0 if defaultValues=true or random
  /// if defaultValues=false, see makePedestalStore) and
  /// store them into CDB located at cdbpath, with a validity period
  /// ranging from startRun to endRun
  
  AliMUONVStore* pedestalStore = Create2DMap();
  Int_t ngenerated = MakePedestalStore(*pedestalStore,defaultValues);
  AliInfo(Form("Ngenerated = %d",ngenerated));
  WriteToCDB("MUON/Calib/Pedestals",pedestalStore,startRun,endRun,defaultValues);
  delete pedestalStore;
}

//_____________________________________________________________________________
void 
AliMUONCDB::WriteKillMap(Bool_t defaultValues,
                         Int_t startRun, Int_t endRun)
{
  /// generate kill map values (either empty one if defaultValues=true, or
  /// random one, see MakeKillMapStore) and
  /// store them into CDB located at cdbpath, with a validity period
  /// ranging from startRun to endRun
  
  AliMUONVStore* killMapStore = Create2DMap();
  if ( !defaultValues )
  {
    Int_t ngenerated = MakeKillMapStore(*killMapStore);
    AliInfo(Form("Ngenerated = %d",ngenerated));
  }
  WriteToCDB("MUON/Calib/KillMap",killMapStore,startRun,endRun,defaultValues);
  delete killMapStore;
}


//_____________________________________________________________________________
void 
AliMUONCDB::WriteGains(Bool_t defaultValues,
                       Int_t startRun, Int_t endRun)
{
  /// generate gain values (either 1 if defaultValues=true or random
  /// if defaultValues=false, see makeGainStore) and
  /// store them into CDB located at cdbpath, with a validity period
  /// ranging from startRun to endRun
  
  AliMUONVStore* gainStore = Create2DMap();
  Int_t ngenerated = MakeGainStore(*gainStore,defaultValues);
  AliInfo(Form("Ngenerated = %d",ngenerated));  
  WriteToCDB("MUON/Calib/Gains",gainStore,startRun,endRun,defaultValues);
  delete gainStore;
}

//_____________________________________________________________________________
void 
AliMUONCDB::WriteCapacitances(const char* filename,
                              Int_t startRun, Int_t endRun)
{
  /// read manu capacitance and injection gain values from file 
  /// and store them into CDB located at cdbpath, with a validity period
  /// ranging from startRun to endRun
  
  AliMUONVStore* capaStore = new AliMUON1DMap(16828);
  Int_t ngenerated = MakeCapacitanceStore(*capaStore,filename);
  AliInfo(Form("Ngenerated = %d",ngenerated));
  if ( ngenerated > 0 ) 
  {
    WriteToCDB("MUON/Calib/Capacitances",capaStore,startRun,endRun,filename);
  }
  delete capaStore;
}

//_____________________________________________________________________________
void 
AliMUONCDB::WriteCapacitances(Bool_t defaultValues,
                              Int_t startRun, Int_t endRun)
{
  /// generate manu capacitance values (either 1 if defaultValues=true or random
  /// if defaultValues=false, see makeCapacitanceStore) and
  /// store them into CDB located at cdbpath, with a validity period
  /// ranging from startRun to endRun
  
  AliMUONVStore* capaStore = new AliMUON1DMap(16828);
  Int_t ngenerated = MakeCapacitanceStore(*capaStore,defaultValues);
  AliInfo(Form("Ngenerated = %d",ngenerated));
  WriteToCDB("MUON/Calib/Capacitances",capaStore,startRun,endRun,defaultValues);
  delete capaStore;
}

//_____________________________________________________________________________
void
AliMUONCDB::WriteTrigger(Bool_t defaultValues, Int_t startRun, Int_t endRun)
{
  /// Writes all Trigger related calibration to CDB
  WriteTriggerDCS(defaultValues,startRun,endRun);
  WriteLocalTriggerMasks(startRun,endRun);
  WriteRegionalTriggerConfig(startRun,endRun);
  WriteGlobalTriggerConfig(startRun,endRun);
  WriteTriggerLut(startRun,endRun);
  WriteTriggerEfficiency(startRun,endRun);
}

//_____________________________________________________________________________
void
AliMUONCDB::WriteTracker(Bool_t defaultValues, Int_t startRun, Int_t endRun)
{
  /// Writes all Tracker related calibration to CDB
  WriteHV(defaultValues,startRun,endRun);
  WritePedestals(defaultValues,startRun,endRun);
  WriteGains(defaultValues,startRun,endRun);
  WriteCapacitances(defaultValues,startRun,endRun);
  WriteNeighbours(startRun,endRun);
  WriteKillMap(startRun,endRun,defaultValues);
}

