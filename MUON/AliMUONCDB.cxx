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
/// \namespace AliMUONCDB
///
/// Helper functions to experience the OCDB
///
/// They allow to read magnetic field, mapping and recoParam from OCDB
///
/// And also to generate dummy (but complete) containers for all the
/// calibration data types we have for tracker and trigger, and to write
/// them into OCDB.
///
/// For more information, please see READMEcalib
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMUONCDB.h"

#include "AliMUON1DArray.h"
#include "AliMUON1DMap.h"
#include "AliMUON2DMap.h"
#include "AliMUON2DStoreValidator.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONCalibParamNF.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONConstants.h"
#include "AliMUONGlobalCrateConfig.h"
#include "AliMUONLogger.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONPadStatusMapMaker.h"
#include "AliMUONRecoParam.h"
#include "AliMUONRegionalTriggerConfig.h"
#include "AliMUONRejectList.h"
#include "AliMUONTrackerData.h"
#include "AliMUONTrackerIO.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerLut.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVStore.h"

#include "AliMpCDB.h"
#include "AliMpConstants.h"
#include "AliMpDEStore.h"
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
#include "AliGRPManager.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliMpBusPatch.h"

#include <Riostream.h>
#include <TArrayI.h>
#include <TClass.h>
#include <TFile.h>
#include <TH1F.h>
#include <TList.h>
#include <TMap.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TMath.h>
#include <TGeoGlobalMagField.h>
#include <TClonesArray.h>
#include <sstream>
#include <set>

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
        
    if (!seg) continue;
    
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
Bool_t AliMUONCDB::CheckOCDB(Bool_t pathOnly)
{
  /// Check that OCDB path and run number are properly set
  
  AliCDBManager* man = AliCDBManager::Instance();
  
  // first OCDB path
  if (!man->IsDefaultStorageSet()) {
    AliErrorGeneral("AliMUONCDB", "OCDB path must be properly set");
    return kFALSE;
  }
  
  // then run number if required
  if (pathOnly) return kTRUE;
  if (man->GetRun() < 0) {
    AliErrorGeneral("AliMUONCDB", "Run number must be properly set");
    return kFALSE;
  }
  
  return kTRUE;
  
}

//_____________________________________________________________________________
Bool_t AliMUONCDB::CheckMapping(Bool_t segmentationOnly)
{
  /// Check that the mapping has been loaded
  
  // first the segmentation
  if (!AliMpSegmentation::Instance(false)) {
    AliErrorGeneral("AliMUONCDB", "Mapping segmentation must be loaded first");
    return kFALSE;
  }
  
  // then the others if required
  if (segmentationOnly) return kTRUE;
  if (!AliMpDDLStore::Instance(false) || !AliMpDEStore::Instance(false) || !AliMpManuStore::Instance(false)) {
    AliErrorGeneral("AliMUONCDB", "Full mapping must be loaded first");
    return kFALSE;
  }
  
  return kTRUE;
  
}

//_____________________________________________________________________________
Bool_t AliMUONCDB::LoadField()
{
  /// Load magnetic field (existing field will be deleted).
  /// OCDB path and run number are supposed to be set.
  
  AliInfoGeneral("AliMUONCDB","Loading field map from GRP...");
  
  if (!AliMUONCDB::CheckOCDB()) return kFALSE;
  
  AliGRPManager grpMan;
  
  // make sure the old field is deleted even if it is locked
  if(TGeoGlobalMagField::Instance()->IsLocked()) delete TGeoGlobalMagField::Instance();
  
  if (!grpMan.ReadGRPEntry() || !grpMan.SetMagField()) {
    AliErrorGeneral("AliMUONCDB", "failed to load magnetic field from OCDB");
    return kFALSE;
  }
  
  return kTRUE;
  
}

//_____________________________________________________________________________
Bool_t AliMUONCDB::LoadMapping(Bool_t segmentationOnly)
{
  /// Load mapping (existing mapping will be unloaded).
  /// OCDB path and run number are supposed to be set.
  
  AliInfoGeneral("AliMUONCDB", "Loading mapping from OCDB...");
  
  if (!AliMUONCDB::CheckOCDB()) return kFALSE;
  
  // in case it has already been set
  AliMpCDB::UnloadAll();
  
  if (segmentationOnly) {
    
    if (!AliMpCDB::LoadMpSegmentation(kTRUE)){
      AliErrorGeneral("AliMUONCDB", "failed to load segmentation from OCDB");
      return kFALSE;
    }
    
  } else {
    
    if (!AliMpCDB::LoadAll(kTRUE)) {
      AliErrorGeneral("AliMUONCDB", "failed to load mapping from OCDB");
      return kFALSE;
    }
    
  }
  
  return kTRUE;
  
}

//_____________________________________________________________________________
AliMUONRecoParam* AliMUONCDB::LoadRecoParam()
{
  /// Load and return reconstruction parameters.
  /// OCDB path is supposed to be set.
  
  AliInfoGeneral("AliMUONCDB", "Loading RecoParam from OCDB...");
  
  if (!AliMUONCDB::CheckOCDB()) return 0x0;
  
  AliMUONRecoParam* recoParam = 0x0;
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("MUON/Calib/RecoParam");
  
  if(entry) {
    
    // load recoParam according OCDB content (single or array)
    if (!(recoParam = dynamic_cast<AliMUONRecoParam*>(entry->GetObject()))) {
      
      TObjArray* recoParamArray = static_cast<TObjArray*>(entry->GetObject());
//      recoParamArray->SetOwner(kTRUE); // FIXME: this should be done, but is causing a problem at the end of the reco... investigate why...
      
      for(Int_t i = 0; i < recoParamArray->GetEntriesFast(); i++) {
	recoParam = static_cast<AliMUONRecoParam*>(recoParamArray->UncheckedAt(i));
	if (recoParam->IsDefault()) break;
	recoParam = 0x0;
      }
      
    }
    
  }
  
  if (!recoParam) AliErrorGeneral("AliMUONCDB", "failed to load RecoParam from OCDB");
  
  return recoParam;
  
}

//_____________________________________________________________________________
TClonesArray* AliMUONCDB::LoadAlignmentData()
{
  /// Load and return the array of alignment objects.
  
  AliInfoGeneral("AliMUONCDB", "Loading Alignemnt from OCDB...");
  
  if (!AliMUONCDB::CheckOCDB()) return 0x0;
  
  TClonesArray* alignmentArray = 0x0;
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("MUON/Align/Data");
  
  if (entry) {
    // load alignement array
    alignmentArray = dynamic_cast<TClonesArray*>(entry->GetObject());
  }
  
  if (!alignmentArray) { 
    AliErrorGeneral("AliMUONCDB", "failed to load Alignemnt from OCDB");
  }  
  
  return alignmentArray;
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
    AliErrorGeneral("AliMUONCDB", Form("opt %s not supported. Only ABS, REL, PERCENT are",opt));
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
TH1** 
AliMUONCDB::Plot(const AliMUONVStore& store, const char* name, Int_t nbins)
{
  /// Make histograms of each dimension of the AliMUONVCalibParam
  /// contained inside store.
  /// It produces histograms named name_0, name_1, etc...
  
  if (!AliMUONCDB::CheckMapping(kTRUE)) return 0x0;
  
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
        AliInfoGeneral("AliMUONCDB", Form("Created histogram %s",h[i]->GetName()));
      }
      delete [] xmin;
      delete [] xmax;
    }
    
    Int_t detElemId = param->ID0();
    Int_t manuId = param->ID1();
    Int_t station = AliMpDEManager::GetChamberId(detElemId)/2;
    
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
    
    if (!seg) continue;
    
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
    AliInfoGeneral("AliMUONCDB", Form("Station %d %d ",(i+1),nPerStation[i]));
  }

  AliInfoGeneral("AliMUONCDB", Form("Number of channels = %d",n));
  
  delete[] nPerStation;
  
  return h;
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeHVStore(TMap& aliasMap, Bool_t defaultValues)
{
  /// Create a HV store
  
  if (!AliMUONCDB::CheckMapping()) return 0;
  
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
  
  AliInfoGeneral("AliMUONCDB", Form("%d HV channels and %d switches",nChannels,nSwitch));
  
  return nChannels+nSwitch;
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeTriggerDCSStore(TMap& aliasMap, Bool_t defaultValues)
{
  /// Create a Trigger HV and Currents store
  
  if (!AliMUONCDB::CheckMapping()) return 0;
  
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
    if ( measureType < 0 ) {
      AliErrorGeneralStream("AliMUONCDB") 
        << "Failed to get DCS variable from an alias (trigger): "
        << aliasName.Data() << endl;
      return 0;
    }
        
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
  
  AliInfoGeneral("AliMUONCDB", Form("Trigger channels I -> %i   HV -> %i",nChannels[0], nChannels[1]));
  
  return nChannels[0] + nChannels[1];
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakePedestalStore(AliMUONVStore& pedestalStore, Bool_t defaultValues)
{
  /// Create a pedestal store. if defaultValues=true, ped.mean=ped.sigma=1,
  /// otherwise mean and sigma are from a gaussian (with parameters
  /// defined below by the kPedestal* constants)

  AliCodeTimerAutoGeneral("",0);
  
  if (!AliMUONCDB::CheckMapping()) return 0;
  
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
      AliErrorGeneral("AliMUONCDB", Form("Could not set DetElemId=%d manuId=%d",detElemId,manuId));
    }
  }
  
  AliInfoGeneral("AliMUONCDB", Form("%d Manus and %d channels.",nmanus,nchannels));
  return nchannels;
}

//_____________________________________________________________________________
AliMUONRejectList* 
AliMUONCDB::MakeRejectListStore(Bool_t defaultValues)
{
  /// Create a reject list
  
  AliCodeTimerAutoGeneral("",0);

  AliMUONRejectList* rl = new AliMUONRejectList;
  
  if (!defaultValues)
  {
    rl->SetDetectionElementProbability(510);
    rl->SetDetectionElementProbability(508);
    return rl;
  }
  
  return rl;
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeOccupancyMapStore(AliMUONVStore& occupancyMapStore, Bool_t defaultValues)
{
  /// Create an occupancy map.
  
  AliCodeTimerAutoGeneral("",0);
  
  if (!AliMUONCDB::CheckMapping()) return 0;
  
  Int_t nmanus(0);
  
  Int_t detElemId;
  Int_t manuId;
  
  AliMpManuIterator it;
  
  Int_t nevents(1000);
  
  while ( it.Next(detElemId,manuId) )
  {
    ++nmanus;
    
    AliMUONVCalibParam* occupancy = new AliMUONCalibParamND(5,1,detElemId,manuId,0);
    
    Double_t occ = 0.0;

    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    
    Int_t numberOfChannelsInManu = de->NofChannelsInManu(manuId);
    
    if (!defaultValues) occ = gRandom->Rndm(1);

    Double_t sumn = occ*nevents;
    
    occupancy->SetValueAsFloat(0,0,sumn); 
    occupancy->SetValueAsFloat(0,1,sumn);
    occupancy->SetValueAsFloat(0,2,sumn);
    occupancy->SetValueAsInt(0,3,numberOfChannelsInManu);
    occupancy->SetValueAsInt(0,4,nevents);
    
    Bool_t ok = occupancyMapStore.Add(occupancy);
    if (!ok)
    {
      AliErrorGeneral("AliMUONCDB", Form("Could not set DetElemId=%d manuId=%d",detElemId,manuId));
    }
  }
  
  return nmanus;
}

//_____________________________________________________________________________
Int_t
AliMUONCDB::MakeCapacitanceStore(AliMUONVStore& capaStore, const char* file)
{
  /// Read the capacitance values from file and append them to the capaStore
  
  if (!AliMUONCDB::CheckMapping()) return 0;
  
  return AliMUONTrackerIO::ReadCapacitances(file,capaStore);
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeCapacitanceStore(AliMUONVStore& capaStore, Bool_t defaultValues)
{
  /// Create a capacitance store. if defaultValues=true, all capa are 1.0,
  /// otherwise they are from a gaussian with parameters defined in the
  /// kCapa* constants below.

  AliCodeTimerAutoGeneral("",0);
  
  if (!AliMUONCDB::CheckMapping()) return 0;
  
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
        AliErrorGeneral("AliMUONCDB", Form("Could not set serialNumber=%d manuId=%d",serialNumber,manuId));
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
  AliInfoGeneral("AliMUONCDB", Form("%5d manus with serial number (out of %5d manus = %3.0f%%)",
               nmanusOK,nmanus,percent));
  AliInfoGeneral("AliMUONCDB", Form("%5d channels",nchannels));
  if ( percent < 100 ) 
  {
    AliWarningGeneral("AliMUONCDB", "Did not get all serial numbers. capaStore is incomplete !!!!");
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
  
  AliCodeTimerAutoGeneral("",0);
  
  if (!AliMUONCDB::CheckMapping()) return 0;
  
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
      AliErrorGeneral("AliMUONCDB", Form("Could not set DetElemId=%d manuId=%d",detElemId,manuId));
    }
  }
  
  AliInfoGeneral("AliMUONCDB", Form("%d Manus and %d channels.",nmanus,nchannels));
  return nchannels;
}

//_____________________________________________________________________________
Int_t
AliMUONCDB::MakeLocalTriggerMaskStore(AliMUONVStore& localBoardMasks)
{
  /// Generate local trigger masks store. All masks are set to FFFF
  
  AliCodeTimerAutoGeneral("",0);
  
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
AliMUONCDB::MakeRegionalTriggerConfigStore(AliMUONRegionalTriggerConfig& rtm)
{
  /// Make a regional trigger config store. Mask is set to FFFF for each local board (Ch.F.)
  
  AliCodeTimerAutoGeneral("",0);
  
  if ( ! rtm.ReadData(AliMpFiles::LocalTriggerBoardMapping()) ) {
    AliErrorGeneral("AliMUONCDB", "Error when reading from mapping file");
    return 0;
  }
    
  return rtm.GetNofTriggerCrates();  
}


//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeGlobalTriggerConfigStore(AliMUONGlobalCrateConfig& gtm)
{
  /// Make a global trigger config store. All masks (disable) set to 0x00 for each Darc board (Ch.F.)
  
  AliCodeTimerAutoGeneral("",0);
  
  return gtm.ReadData(AliMpFiles::GlobalTriggerBoardMapping());
}


//_____________________________________________________________________________
AliMUONTriggerLut* 
AliMUONCDB::MakeTriggerLUT(const char* file)
{
  /// Make a triggerlut object, from a file.
  
  AliCodeTimerAutoGeneral("",0);
  
  AliMUONTriggerLut* lut = new AliMUONTriggerLut;
  lut->ReadFromFile(file);
  return lut;
}

//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells*
AliMUONCDB::MakeTriggerEfficiency(const char* file)
{
  /// Make a trigger efficiency object from a file.
  
  AliCodeTimerAutoGeneral("",0);
  
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
  
  if (!AliMUONCDB::CheckOCDB(kTRUE)) return;
  
  AliCDBId id(calibpath,startRun,endRun);
  AliCDBMetaData md;
  md.SetAliRootVersion(gROOT->GetVersion());
  md.SetComment(comment);
  md.SetResponsible(responsible);
  AliCDBManager::Instance()->Put(object,id,&md);
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeNeighbourStore(AliMUONVStore& neighbourStore)
{
  /// Fill the neighbours store with, for each channel, a TObjArray of its
  /// neighbouring pads (including itself)
  
  AliCodeTimerAutoGeneral("",0);
  
  if (!AliMUONCDB::CheckMapping()) return 0;
  
  AliInfoGeneral("AliMUONCDB", "Generating NeighbourStore. This will take a while. Please be patient.");
  
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
        AliErrorGeneral("AliMUONCDB", Form("Could not set DetElemId=%d manuId=%d",detElemId,manuId));
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
AliMUONCDB::WriteLocalTriggerMasks(Int_t startRun, Int_t endRun)
{  
  /// Write local trigger masks to OCDB
  
  AliMUONVStore* ltm = new AliMUON1DArray(AliMpConstants::TotalNofLocalBoards()+1);
  Int_t ngenerated = MakeLocalTriggerMaskStore(*ltm);
  AliInfoGeneral("AliMUONCDB", Form("Ngenerated = %d",ngenerated));
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
  AliInfoGeneral("AliMUONCDB", Form("Ngenerated = %d",ngenerated));
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
  AliInfoGeneral("AliMUONCDB", Form("Ngenerated = %d",ngenerated));
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
  AliInfoGeneral("AliMUONCDB", Form("Ngenerated = %d",ngenerated));
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
  AliInfoGeneral("AliMUONCDB", Form("Ngenerated = %d",ngenerated));
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
  AliInfoGeneral("AliMUONCDB", Form("Ngenerated = %d",ngenerated));
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
  AliInfoGeneral("AliMUONCDB", Form("Ngenerated = %d",ngenerated));
  WriteToCDB("MUON/Calib/Pedestals",pedestalStore,startRun,endRun,defaultValues);
  delete pedestalStore;
}

//_____________________________________________________________________________
void 
AliMUONCDB::WriteOccupancyMap(Bool_t defaultValues,
                              Int_t startRun, Int_t endRun)
{
  /// generate occupancy map values (either empty one if defaultValues=true, or
  /// random one, see MakeOccupancyMapStore) and
  /// store them into CDB located at cdbpath, with a validity period
  /// ranging from startRun to endRun
  
  AliMUONVStore* occupancyMapStore = Create2DMap();
  Int_t ngenerated = MakeOccupancyMapStore(*occupancyMapStore,defaultValues);
  AliInfoGeneral("AliMUONCDB", Form("Ngenerated = %d",ngenerated));
  WriteToCDB("MUON/Calib/OccupancyMap",occupancyMapStore,startRun,endRun,defaultValues);
  delete occupancyMapStore;
}

//_____________________________________________________________________________
void 
AliMUONCDB::WriteRejectList(Bool_t defaultValues,
                              Int_t startRun, Int_t endRun)
{
  /// generate reject list values (either empty one if defaultValues=true, or
  /// random one, see MakeRejectListStore) and
  /// store them into CDB located at cdbpath, with a validity period
  /// ranging from startRun to endRun
  
  AliMUONRejectList* rl = MakeRejectListStore(defaultValues);
  WriteToCDB("MUON/Calib/RejectList",rl,startRun,endRun,defaultValues);
  delete rl;
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
  AliInfoGeneral("AliMUONCDB", Form("Ngenerated = %d",ngenerated));  
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
  AliInfoGeneral("AliMUONCDB", Form("Ngenerated = %d",ngenerated));
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
  AliInfoGeneral("AliMUONCDB", Form("Ngenerated = %d",ngenerated));
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
AliMUONCDB::WriteConfig(Int_t startRun, Int_t endRun)
{
  /// Write complete tracker configuration to OCDB
  ostringstream lines;
  TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  AliMpBusPatch* bp;
  while ( ( bp = static_cast<AliMpBusPatch*>(next()) ) )
  {
    for (Int_t imanu = 0; imanu < bp->GetNofManus(); ++imanu) 
    {
      lines << bp->GetId() << " " << bp->GetManuId(imanu) << endl;
    }
  }
  
  AliMUON2DMap config(kTRUE);
  
  AliMUONTrackerIO::DecodeConfig(lines.str().c_str(),config);
  
  WriteToCDB("MUON/Calib/Config",&config,startRun,endRun,kTRUE);
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
  WriteOccupancyMap(defaultValues,startRun,endRun);
  WriteRejectList(defaultValues,startRun,endRun);
  WriteConfig(startRun,endRun);
}

//_____________________________________________________________________________
void 
AliMUONCDB::ShowCapacitances()
{
  /// Show briefly the number of capa values we have in the OCDB,
  /// and the list of manu that are actually in the config and for which
  /// we miss the capa (if any).
  
  if (!AliMUONCDB::CheckOCDB()) return;
  
  AliMUONCDB::LoadMapping();
  
  if (!AliMUONCDB::CheckMapping()) return;
  
  AliCDBEntry* e = AliCDBManager::Instance()->Get("MUON/Calib/Config");
  
  if (!e) return ;
  
  AliMUONVStore* config = static_cast<AliMUONVStore*>(e->GetObject());
  
  e = AliCDBManager::Instance()->Get("MUON/Calib/Capacitances");
  
  if (!e) return;
  
  AliMUONVStore* capacitances = static_cast<AliMUONVStore*>(e->GetObject());
  
  AliInfoGeneral("ShowCapacitances",Form("%d capacitances are in OCDB",capacitances->GetSize()));
  
  TIter nextManu(config->CreateIterator());
  AliMUONVCalibParam* param;
  
  while ( ( param = static_cast<AliMUONVCalibParam*>(nextManu()) ) )
  {
    Int_t detElemId = param->ID0();
    Int_t manuId = param->ID1();
    
    Int_t serialNumber 
    = AliMpManuStore::Instance()->GetManuSerial(detElemId, manuId);
    
    if (serialNumber<0)
    {
      AliErrorGeneral("ShowCapacitances",Form("Did not find serial for DE %04d MANUID %04d",detElemId,manuId));
    }
    else
    {
      AliMUONVCalibParam* capa = static_cast<AliMUONVCalibParam*>(capacitances->FindObject(serialNumber));
      if (!capa)
      {
        AliErrorGeneral("ShowCapacitances",Form("Did not find capacitance for DE %04d MANUID %04d SERIAL %d",detElemId,manuId,serialNumber));
      }
    }
  }
  
}

//_____________________________________________________________________________
void 
AliMUONCDB::ShowConfig(Bool_t withStatusMap)
{  
  /// Dumps the current tracker configuration, i.e. number and identity of missing buspatches
  /// If statusMap is true, will also take into account the status map to report the number
  /// of good channels
  
  if (!AliMUONCDB::CheckOCDB()) return;
  
  AliMUONCDB::LoadMapping();
  
  if (!AliMUONCDB::CheckMapping()) return;
  
  AliCDBEntry* e = AliCDBManager::Instance()->Get("MUON/Calib/Config");
  
  if (!e) return ;
  
  AliMUONVStore* config = static_cast<AliMUONVStore*>(e->GetObject());
  
  AliMUONPadStatusMapMaker* statusMapMaker(0x0);
  AliMUONCalibrationData* cd(0x0);
  AliMUONPadStatusMaker* statusMaker(0x0);

  if ( withStatusMap ) 
  {
    cd = new AliMUONCalibrationData(AliCDBManager::Instance()->GetRun());
  
    statusMaker = new AliMUONPadStatusMaker(*cd);

    AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  
    statusMaker->SetLimits(*recoParam);
  
    UInt_t mask = recoParam->PadGoodnessMask();

    delete recoParam;
    
    const Bool_t deferredInitialization = kFALSE;
  
    statusMapMaker = new AliMUONPadStatusMapMaker(*cd,mask,deferredInitialization);
  }
  
  TIter nextManu(config->CreateIterator());
  AliMUONVCalibParam* param;
  
  AliMpExMap buspatches;
  
  while ( ( param = static_cast<AliMUONVCalibParam*>(nextManu()) ) )
  {
    Int_t detElemId = param->ID0();
    Int_t manuId = param->ID1();
    Int_t busPatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);
    if ( buspatches.GetValue(busPatchId) == 0x0 ) 
    {      
      buspatches.Add(busPatchId,new TObjString(Form("BP%04d",busPatchId)));
    }
  }

  TArrayI removed(buspatches.GetSize());

  TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  AliMpBusPatch* bp;
  Int_t n(0);
  Int_t nok(0);
  Int_t nremoved(0);
  
  // accounting of bus patches first
  
  while ( ( bp = static_cast<AliMpBusPatch*>(next())))
  {
    if ( buspatches.GetValue(bp->GetId()) )
    {
      ++nok;
    }
    else
    {
      removed.SetAt(bp->GetId(),nremoved++);
    }
  }
  
  // accounting of channels
  
  AliMpManuIterator it;

  Int_t totalNumberOfChannels(0);
  Int_t removedChannels(0);
  Int_t badChannels(0);
  Int_t badAndRemovedChannels(0);
  Int_t badOrRemovedChannels(0);
  
  Int_t detElemId, manuId;
  
  while ( it.Next(detElemId,manuId) )
  {
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    for ( Int_t i = 0; i < AliMpConstants::ManuNofChannels(); ++i )
    {
      Int_t busPatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);
      
      if ( de->IsConnectedChannel(manuId,i) )
      {
        ++totalNumberOfChannels;
        Bool_t badBusPatch = ( buspatches.GetValue(busPatchId) == 0x0 );
        
        if ( withStatusMap ) 
        {
          Bool_t badChannel = ( ( statusMapMaker->StatusMap(detElemId,manuId,i) & AliMUONPadStatusMapMaker::SelfDeadMask() ) != 0);
          if ( badChannel ) ++badChannels;
          if ( badBusPatch && badChannel ) ++badAndRemovedChannels;
          if ( badBusPatch || badChannel ) ++badOrRemovedChannels;          
        }
      
        if ( badBusPatch) ++removedChannels;
      }
    }
  }
  
  
  Int_t* indices = new Int_t[nremoved];
  
  TMath::Sort(nremoved,removed.GetArray(),indices,kFALSE);
  
  for ( Int_t i = 0; i < nremoved; ++i ) 
  {
    Int_t busPatchId = removed[indices[i]];
    bp = AliMpDDLStore::Instance()->GetBusPatch(busPatchId);
    bp->Print();
  }
  
  delete[] indices;  
  
  cout << endl;
  cout << Form("Bus patches n=%3d nok=%3d nremoved=%3d",n,nok,nremoved) << endl;

  cout << Form("Channels n=%6d nremoved=%6d bad=%6d bad and removed=%6d bad or removed=%6d",
               totalNumberOfChannels,removedChannels,badChannels,badAndRemovedChannels,badOrRemovedChannels) << endl;
  
  if (totalNumberOfChannels>0)
  {
    cout << Form("Percentage of readout channels %5.1f %%",removedChannels*100.0/totalNumberOfChannels) << endl;
    if ( withStatusMap )
    {
      cout << Form("Percentage of non useable channels (bad or removed) %5.1f %%",
                   badOrRemovedChannels*100.0/totalNumberOfChannels) << endl;
    }
  }
  
  
  delete statusMapMaker;
  delete cd;
  delete statusMaker;  
}

//______________________________________________________________________________
void AliMUONCDB::ReadIntegers(const char* filename, std::vector<int>& integers)
{
  /// Read integers from filename, where integers are either
  /// separated by "," or by return carriage
  ifstream in(gSystem->ExpandPathName(filename));
  int i;
  
  std::set<int> runset;
  
  char line[10000];
  
  in.getline(line,10000,'\n');
  
  TString sline(line);
  
  if (sline.Contains(","))
  {
    TObjArray* a = sline.Tokenize(",");
    TIter next(a);
    TObjString* s;
    while ( ( s = static_cast<TObjString*>(next()) ) )
    {
      runset.insert(s->String().Atoi());
    }
    delete a;
  }
  else
  {
    runset.insert(sline.Atoi());
    
    while ( in >> i )
    {
      runset.insert(i);
    }
  }
  
  for ( std::set<int>::const_iterator it = runset.begin(); it != runset.end(); ++it ) 
  {
    integers.push_back((*it)); 
  }
  
  std::sort(integers.begin(),integers.end());
}

//______________________________________________________________________________
void AliMUONCDB::ShowFaultyBusPatches(const char* runlist, double occLimit,
                                      const char* outputBaseName,
                                      const char* ocdbPath)
{
  /// Shows the list of bus patches above a given occupancy limit,
  /// for each run in the runlist
  
  AliLog::GetRootLogger()->SetGlobalLogLevel(AliLog::kError);
  
  //  AliLog::SetPrintType(AliLog::kInfo,kFALSE);
  //  AliLog::SetPrintType(AliLog::kWarning,kFALSE);
  //  gErrorIgnoreLevel=kError; // to avoid all the TAlienFile::Open messages...
  
  AliCDBManager* man = AliCDBManager::Instance();
  
  man->SetDefaultStorage(ocdbPath);
  
  Bool_t first(kTRUE);
  
  std::vector<int> runnumbers;
  
  ReadIntegers(runlist,runnumbers);
  
  AliMUON2DMap bpValues(kFALSE);
  
  std::ofstream outfile(Form("%s.txt",outputBaseName));
  
  for ( unsigned int i = 0 ; i < runnumbers.size(); ++i )
  {
    int runNumber = runnumbers[i];
    
    man->SetRun(runNumber);
    
    if ( first ) 
    {
      AliMpCDB::LoadAll();  
      first = kFALSE;
    }
    
    AliCDBEntry* e = man->Get("MUON/Calib/OccupancyMap",runNumber);
    
    if (!e)
    {
      AliErrorGeneral("AliMUONCDB::ShowFaultyBusPatches",
                      Form("Could not get OccupancyMap for run %09d",runNumber));
      continue;
    }
    
    AliMUONVStore* occmap = static_cast<AliMUONVStore*>(e->GetObject());
    
    AliMUONTrackerData td("occ","occ",*occmap);
    
    TIter nextBP(AliMpDDLStore::Instance()->CreateBusPatchIterator());
    AliMpBusPatch* bp;
    std::set<int> buspatches;
    Double_t sumn = 1000.0;
    
    while ( ( bp = static_cast<AliMpBusPatch*>(nextBP()) ) )
    {      
      Double_t occ = td.BusPatch(bp->GetId(),2);
      
      if (occ>occLimit) 
      {
        buspatches.insert(bp->GetId());
        
        AliMUONVCalibParam* param = static_cast<AliMUONVCalibParam*>(bpValues.FindObject(bp->GetId()));
        
        if (!param)
        {
          param = new AliMUONCalibParamND(5, 1, bp->GetId(), 0);
          bpValues.Add(param);
          
          Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(bp->GetId());
          AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
          
          Int_t nchannels(0);
          
          for ( Int_t imanu = 0; imanu < bp->GetNofManus(); ++imanu ) 
          {
            Int_t manuId = bp->GetManuId(imanu);
            nchannels += de->NofChannelsInManu(manuId);
          }
          
          param->SetValueAsDouble(0,2,sumn);
          param->SetValueAsDouble(0,3,nchannels);
          param->SetValueAsDouble(0,4,1);          
        }
        
        Double_t sumw = sumn*(param->ValueAsDouble(0)/sumn+1.0/runnumbers.size());
        Double_t sumw2 = 0.0; //(sumn-1)*ey*ey+sumw*sumw/sumn;
        
        param->SetValueAsDouble(0,0,sumw);
        param->SetValueAsDouble(0,1,sumw2);
        
      }
    }
    
    outfile << Form("RUN %09d",runNumber);
    
    for ( std::set<int>::const_iterator bit = buspatches.begin(); bit != buspatches.end(); ++bit )
    {
      outfile << Form(" %4d",*bit);
    }
    outfile << endl;
  }
  
  outfile.close();
  
  const char* name = "BPfailureRate";
  
  AliMUONTrackerData* mpData = new AliMUONTrackerData(name,name,bpValues,2);
  mpData->SetDimensionName(0,name);
  
  TFile f(Form("%s.root",outputBaseName),"recreate");
  mpData->Write();
  f.Close();
  
  cout << Form("Results are in %s.txt and %s.root",outputBaseName,outputBaseName) << endl;
  
  gSystem->Exec(Form("cat %s.txt",outputBaseName));
  
}

//______________________________________________________________________________
void AliMUONCDB::CheckHV(Int_t runNumber, Int_t verbose)
{
  /// Check the HV values in OCDB for a given run
  
  TList messages;
  messages.SetOwner(kTRUE);
  
  Bool_t patched(kTRUE);
  
  if (!AliCDBManager::Instance()->IsDefaultStorageSet())
  {
    AliCDBManager::Instance()->SetDefaultStorage("raw://");
  }

  AliCDBManager::Instance()->SetRun(runNumber);

  LoadMapping();
  
  AliMUONCalibrationData::CreateHV(runNumber,0,patched,&messages);
  
  AliMUONCalibrationData cd(runNumber,true);
  
  AliMUONPadStatusMaker statusMaker(cd);
  
  AliMUONRecoParam* rp = AliMUONCDB::LoadRecoParam();
  
  if (!rp)
  {
    AliErrorGeneral("AliMUONCDB::CheckHV","Could not get RecoParam !!!");
    return;
  }
  
  statusMaker.SetLimits(*rp);
  
  TIter next(&messages);
  TObjString* s;
  AliMpDCSNamer hvNamer("TRACKER");
  AliMUONLogger log;
  
  while ( ( s = static_cast<TObjString*>(next()) ) )
  {
    TObjArray* a = s->String().Tokenize(":");
    
    TString name(static_cast<TObjString*>(a->At(0))->String());
    
    if ( name.Contains("sw") || name.Contains("SUMMARY") ) continue;
    
    Int_t index = hvNamer.DCSIndexFromDCSAlias(name.Data());
    
    Int_t detElemId = hvNamer.DetElemIdFromDCSAlias(name.Data());
    
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
    
    if (!de)
    {
      AliErrorGeneral("AliMUONCDB::CheckHV",Form("Could not get detElemId from dcsAlias %s",name.Data()));
      continue;
    }
    
    Int_t manuId;
    
    if ( index >= 0 )
    {
      const AliMpArrayI* array = de->ManusForHV(index);
      manuId = array->GetValue(0);
    }
    else
      
    {
      AliMpBusPatch* bp = AliMpDDLStore::Instance()->GetBusPatch(de->GetBusPatchId(0));
      
      manuId = bp->GetManuId(0);
    }
    
    Int_t status = statusMaker.HVStatus(detElemId,manuId);
    
    log.Log(AliMUONPadStatusMaker::AsString(status).Data());
    
    s->String() += Form(" (DE %4d) ",detElemId);
    s->String() += AliMUONPadStatusMaker::AsString(status).Data();
    
    delete a;
  }    
  
  TIter nextMessage(&messages);
  TObjString* msg;
  
  while ( ( msg = static_cast<TObjString*>(nextMessage()) ) )
  {
    if ( verbose > 0 || msg->String().Contains("SUMMARY") )
    {
      AliInfoGeneral("AliMUONCDB::CheckHV",Form("RUN %09d HVchannel %s",runNumber,msg->String().Data()));
    }
  }
  
  TString lmsg;
  Int_t occurance;
  TString totalLog;
  
  while (log.Next(lmsg,occurance))
  {
    totalLog += Form("%s(%d)",lmsg.Data(),occurance);
    totalLog += " | ";
  }

  AliInfoGeneral("AliMUONCDB::CheckHV",Form("RUN %09d %s",runNumber,totalLog.Data()));

  // one last loop to get the list of problematic HV channels
  nextMessage.Reset();
  
  while ( ( msg = static_cast<TObjString*>(nextMessage()) ) )
  {
    if ( msg->String().Contains("HV ") )
    {
      AliInfoGeneral("AliMUONCDB::CheckHV",Form("     Problem at %s",msg->String().Data()));      
    }
  }
   
  AliCDBManager::Instance()->ClearCache();
}


