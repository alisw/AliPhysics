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

#include "AliMUONCDB.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliDCSValue.h"
#include "AliLog.h"
#include "AliMUON1DArray.h"
#include "AliMUON1DMap.h"
#include "AliMUON2DMap.h"
#include "AliMUON2DStoreValidator.h"
#include "AliMUONCalibParamNF.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONConstants.h"
#include "AliMUONHVNamer.h"
#include "AliMUONObjectPair.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerLut.h"
#include "AliMUONV2DStore.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVDataIterator.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpManuList.h"
#include "AliMpSegmentation.h"
#include "AliMpStationType.h"
#include "AliMpVSegmentation.h"
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

/// \cond CLASSIMP
ClassImp(AliMUONCDB)
/// \endcond

const Int_t AliMUONCDB::fgkMaxNofChannelsPerManu=64;

namespace
{
//_____________________________________________________________________________
void getBoundaries(const AliMUONV2DStore& store,
                   Float_t& x0min, Float_t& x0max,
                   Float_t& x1min, Float_t& x1max)
{
  x0min=1E30;
  x0max=-1E30;
  x1min=1E30;
  x1max=-1E30;
  
  AliMUONVDataIterator* it = store.Iterator();
  
  AliMUONObjectPair* p;
  
  while ( ( p = dynamic_cast<AliMUONObjectPair*>(it->Next() ) ) )
  {
    AliMpIntPair* dm = dynamic_cast<AliMpIntPair*>(p->Key());
    AliMUONVCalibParam* value = dynamic_cast<AliMUONVCalibParam*>(p->Value());
    
    Int_t detElemId = dm->GetFirst();
    Int_t manuId = dm->GetSecond();
    
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
        
    for ( Int_t manuChannel = 0; manuChannel < value->Size(); ++manuChannel )
    {
      AliMpPad pad = seg->PadByLocation(AliMpIntPair(manuId,manuChannel),kFALSE);
      if (!pad.IsValid()) continue;
      
      Float_t x0 = value->ValueAsFloat(manuChannel,0);
      
      x0min = TMath::Min(x0min,x0);
      x0max = TMath::Max(x0max,x0);
      if ( value->Dimension()>1 )
      {
        Float_t x1 = value->ValueAsFloat(manuChannel,1);
        x1min = TMath::Min(x1min,x1);
        x1max = TMath::Max(x1max,x1);
      }
    }
    if (it->IsOwner()) delete p;
  }  
}

}

//_____________________________________________________________________________
AliMUONCDB::AliMUONCDB(const char* cdbpath)
: TObject(),
  fCDBPath(cdbpath),
  fManuList(AliMpManuList::ManuList())
{
    /// ctor
}

//_____________________________________________________________________________
AliMUONCDB::~AliMUONCDB()
{
  /// dtor
  delete fManuList;
}

//_____________________________________________________________________________
AliMUONV2DStore* 
AliMUONCDB::Diff(AliMUONV2DStore& store1, AliMUONV2DStore& store2, 
                 const char* opt)
{
  /// creates a store which contains store1-store2
  /// if opt="abs" the difference is absolute one,
  /// if opt="rel" then what is stored is (store1-store2)/store1
  /// WARNING Works only for stores which holds AliMUONVCalibParam objects
  
  TString sopt(opt);
  sopt.ToUpper();
  
  if ( !sopt.Contains("ABS") && !sopt.Contains("REL") )
  {
    AliErrorClass(Form("opt %s not supported. Only ABS or REL are",opt));
    return 0x0;
  }
  
  AliMUONV2DStore* d = static_cast<AliMUONV2DStore*>(store1.Clone());
  
  AliMUONVDataIterator* it = d->Iterator();
  
  AliMUONObjectPair* p;
  
  while ( ( p = dynamic_cast<AliMUONObjectPair*>(it->Next() ) ) )
  {
    AliMpIntPair* dm = dynamic_cast<AliMpIntPair*>(p->Key());
    //FIXME: this might happen (if a full manu is missing, for instance)
    //handle it.
    if (!dm) 
    {
      cerr << "dm is null : FIXME: this might happen !" << endl;
      delete d; 
      return 0;
    }
    AliMUONVCalibParam* param = dynamic_cast<AliMUONVCalibParam*>(p->Value());
    //FIXMENOT: this should *not* happen
    if (!param) 
    {
      cerr << "param is null" << endl;
      delete d;
      return 0;
    }
    
    AliMUONVCalibParam* param2 = dynamic_cast<AliMUONVCalibParam*>(store2.Get(dm->GetFirst(),dm->GetSecond()));
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
        else if ( sopt.Contains("REL") )
        {
          if ( param->ValueAsFloat(i,j) ) 
          {
            value = (param->ValueAsFloat(i,j) - param2->ValueAsFloat(i,j))/param->ValueAsFloat(i,j);
          }
          else 
          {
            continue;
          }
        }
        param->SetValueAsFloat(i,j,value);
      }      
    }
    if (it->IsOwner()) delete p;
  }
  return d;
}

//_____________________________________________________________________________
void 
AliMUONCDB::Plot(const AliMUONV2DStore& store, const char* name, Int_t nbins)
{
  /// Make a plot of the first 1 or 2 dimensions of the AliMUONVCalibParam
  /// contained inside store.
  /// It produces histograms named name_0, name_1, etc...
  
  Float_t x0min, x0max, x1min, x1max;
  
  getBoundaries(store,x0min,x0max,x1min,x1max);
  
  if ( x0min > x0max ) 
  {
    cerr << Form("Something is wrong with boundaries : x0(min,max)=%e,%e",
                 x0min,x0max) << endl;
    return;
  }

  if ( TMath::Abs(x0min-x0max) < 1E-3 ) 
  {
    x0min -= 1;
    x0max += 1;
  }
  
  TH1* h0 = new TH1F(Form("%s_0",name),Form("%s_0",name),
                    nbins,x0min,x0max);
  
  TH1* h1(0);
  
  if ( x1max > x1min )
  {
    h1 = new TH1F(Form("%s_1",name),Form("%s_1",name),
                  nbins,x1min,x1max);
  }
  
  TIter next(fManuList);
  AliMpIntPair* p;
  Int_t n(0);
  Int_t nPerStation[7];
  
  for ( Int_t i = 0; i < 7; ++i ) nPerStation[i]=0;
  
  while ( ( p = (AliMpIntPair*)next() ) )
  {
    Int_t detElemId = p->GetFirst();
    Int_t manuId = p->GetSecond();
    Int_t station = AliMpDEManager::GetChamberId(detElemId);
    
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
    
    AliMUONVCalibParam* value = 
      dynamic_cast<AliMUONVCalibParam*>(store.Get(detElemId,manuId));
    
    if (value)
    {
      for ( Int_t manuChannel = 0; manuChannel < value->Size(); ++manuChannel )
      {
        AliMpPad pad = seg->PadByLocation(AliMpIntPair(manuId,manuChannel),kFALSE);
        if (!pad.IsValid()) continue;

        ++n;
        ++nPerStation[station];
        Float_t x = value->ValueAsFloat(manuChannel,0);
        if ( x>1E4 ) 
        {
          AliInfo(Form("DE %d Manu %d Ch %d x=%e",detElemId,manuId,manuChannel,x));
        }
        h0->Fill(x);
        if (h1)
        {
          h1->Fill(value->ValueAsFloat(manuChannel,1));
        }
      }
    }
    else
    {
      AliWarning(Form("Got a null value for DE=%d manuId=%d",detElemId,manuId));
    }
  }
  
  AliInfo(Form("Number of channels = %d",n));
  for ( Int_t i = 0; i < 7; ++i )
  {
    AliInfo(Form("Station %d %d ",(i+1),nPerStation[i]));
  }
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeHVStore(TMap& aliasMap, Bool_t defaultValues)
{
  /// Create a HV store
  
  AliMUONHVNamer hvNamer;
  TRandom random;
  
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
        Float_t r = random.Uniform();
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
        if (!defaultValues) value = random.Gaus(1750,62.5);
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
AliMUONCDB::MakePedestalStore(AliMUONV2DStore& pedestalStore, Bool_t defaultValues)
{
  /// Create a pedestal store. if defaultValues=true, ped.mean=ped.sigma=1,
  /// otherwise mean and sigma are from a gaussian (with parameters
  /// defined below by the kPedestal* constants)
  
  TIter next(fManuList);
  
  AliMpIntPair* p;
  
  Int_t nchannels(0);
  Int_t nmanus(0);
  
  Bool_t replace = kFALSE;
  
  const Int_t kChannels(64);
  const Float_t kPedestalMeanMean(200);
  const Float_t kPedestalMeanSigma(10);
  const Float_t kPedestalSigmaMean(1.0);
  const Float_t kPedestalSigmaSigma(0.2);
  
  while ( ( p = (AliMpIntPair*)next() ) )
  {
    ++nmanus;
    AliMUONVCalibParam* ped = new AliMUONCalibParamNF(2,kChannels,AliMUONVCalibParam::InvalidFloatValue());
    
    Int_t detElemId = p->GetFirst();
        
    Int_t manuId = p->GetSecond();
    
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
    
    for ( Int_t manuChannel = 0; manuChannel < kChannels; ++manuChannel )
    {
      AliMpPad pad = seg->PadByLocation(AliMpIntPair(manuId,manuChannel),kFALSE);
      if (!pad.IsValid()) continue;
      
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
        meanPedestal = -1;
        while ( meanPedestal < 0 )
        {
          meanPedestal = gRandom->Gaus(kPedestalMeanMean,kPedestalMeanSigma);
        }
        sigmaPedestal = -1;
        while ( sigmaPedestal < 0 )
        {
          sigmaPedestal = gRandom->Gaus(kPedestalSigmaMean,kPedestalSigmaSigma);
        }
      }
      ped->SetValueAsFloat(manuChannel,0,meanPedestal);
      ped->SetValueAsFloat(manuChannel,1,sigmaPedestal);
      
    }
    Bool_t ok = pedestalStore.Set(detElemId,manuId,ped,replace);
    if (!ok)
    {
      AliError(Form("Could not set DetElemId=%d manuId=%d",detElemId,manuId));
    }
  }
  
  AliInfo(Form("%d Manus and %d channels.",nmanus,nchannels));
  return nchannels;
  
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeCapacitanceStore(AliMUONV1DStore& capaStore, Bool_t defaultValues)
{
  /// Create a capacitance store. if defaultValues=true, all capa are 1.0,
  /// otherwise they are from a gaussian with parameters defined in the
  /// kCapa* constants below.
  
  TIter next(fManuList);
  
  AliMpIntPair* p;
  
  Int_t nchannels(0);
  Int_t nmanus(0);
  Int_t nmanusOK(0); // manus for which we got the serial number
  
  Bool_t replace = kFALSE;
  
  const Float_t kCapaMean(1.0);
  const Float_t kCapaSigma(0.5);
  
  while ( ( p = (AliMpIntPair*)next() ) )
  {
    ++nmanus;
    
    Int_t detElemId = p->GetFirst();
    Int_t manuId = p->GetSecond();
    
    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId); 
    Int_t serialNumber = de->GetManuSerialFromId(manuId);
      
    if ( serialNumber > 0 ) ++nmanusOK;
    
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);

    AliMUONVCalibParam* capa = static_cast<AliMUONVCalibParam*>(capaStore.Get(serialNumber));
    
    if (!capa)
    {
      capa = new AliMUONCalibParamNF(1,fgkMaxNofChannelsPerManu,1.0);
      Bool_t ok = capaStore.Set(serialNumber,capa,replace);
      if (!ok)
      {
        AliError(Form("Could not set serialNumber=%d manuId=%d",serialNumber,manuId));
      }      
    }
    
    for ( Int_t manuChannel = 0; manuChannel < capa->Size(); ++manuChannel )
    {
      AliMpPad pad = seg->PadByLocation(AliMpIntPair(manuId,manuChannel),kFALSE);
      if (!pad.IsValid()) continue;
      
      ++nchannels;
      
      Float_t capaValue;
      
      if ( defaultValues ) 
      {
        capaValue = 1.0;
      }
      else
      {
        capaValue = -1;
        while ( capaValue < 0 )
        {
          capaValue = gRandom->Gaus(kCapaMean,kCapaSigma);
        }
      }
      capa->SetValueAsFloat(manuChannel,0,capaValue);
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
AliMUONCDB::MakeGainStore(AliMUONV2DStore& gainStore, Bool_t defaultValues)
{  
  /// Create a gain store. if defaultValues=true, all gain are 1.0,
  /// otherwise they are from a gaussian with parameters defined in the
  /// kGain* constants below.
  
  TIter next(fManuList);
  
  AliMpIntPair* p;
  
  Int_t nchannels(0);
  Int_t nmanus(0);
  
  Bool_t replace = kFALSE;
  
  const Double_t kSaturation(3000);
  const Double_t kGainMean(1.0);
  const Double_t kGainSigma(0.05);
  
  while ( ( p = (AliMpIntPair*)next() ) )
  {
    ++nmanus;
    AliMUONVCalibParam* gain = 
      new AliMUONCalibParamNF(2,fgkMaxNofChannelsPerManu,
                              AliMUONVCalibParam::InvalidFloatValue());

    Int_t detElemId = p->GetFirst();
    Int_t manuId = p->GetSecond();

    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);

    for ( Int_t manuChannel = 0; manuChannel < gain->Size(); ++manuChannel )
    {
      AliMpPad pad = seg->PadByLocation(AliMpIntPair(manuId,manuChannel),kFALSE);
      if (!pad.IsValid()) continue;
      
      ++nchannels;
      
      Float_t meanGain;
      Float_t saturation(kSaturation);
    
      if ( defaultValues ) 
      {
        meanGain = 1.0;
      }
      else
      {
        meanGain = -1;
        while ( meanGain < 0 )
        {
          meanGain = gRandom->Gaus(kGainMean,kGainSigma);
        }
      }
      gain->SetValueAsFloat(manuChannel,0,meanGain);
      gain->SetValueAsFloat(manuChannel,1,saturation);
      
    }
    Bool_t ok = gainStore.Set(detElemId,manuId,gain,replace);
    if (!ok)
    {
      AliError(Form("Could not set DetElemId=%d manuId=%d",detElemId,manuId));
    }
  }
  
  AliInfo(Form("%d Manus and %d channels.",nmanus,nchannels));
  return nchannels;
}

//_____________________________________________________________________________
Int_t
AliMUONCDB::MakeLocalTriggerMaskStore(AliMUONV1DStore& localBoardMasks) const
{
  /// Generate local trigger masks store. All masks are set to FFFF
  
  Int_t ngenerated(0);
  // Generate fake mask values for 234 localboards and put that into
  // one single container (localBoardMasks)
  for ( Int_t i = 1; i <= 234; ++i )
  {
    AliMUONVCalibParam* localBoard = new AliMUONCalibParamNI(1,8);
    for ( Int_t x = 0; x < 2; ++x )
    {
      for ( Int_t y = 0; y < 4; ++y )
      {
        Int_t index = x*4+y;
        localBoard->SetValueAsInt(index,0,0xFFFF);
        ++ngenerated;
      }
    }
    localBoardMasks.Set(i,localBoard,kFALSE);
  }
  return ngenerated;
}

//_____________________________________________________________________________
Int_t
AliMUONCDB::MakeRegionalTriggerMaskStore(AliMUONV1DStore& rtm) const
{
  /// Make a regional trigger masks store. All masks are set to 3F
  
  Int_t ngenerated(0);
  for ( Int_t i = 0; i < 16; ++i )
  {
    AliMUONVCalibParam* regionalBoard = new AliMUONCalibParamNI(1,16);
    for ( Int_t j = 0; j < 16; ++j )
    {
      regionalBoard->SetValueAsInt(j,0,0x3F);
      ++ngenerated;
    }
    rtm.Set(i,regionalBoard,kFALSE);
  }
  
  return ngenerated;
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeGlobalTriggerMaskStore(AliMUONVCalibParam& gtm) const
{
  /// Make a global trigger masks store. All masks set to FFF
  
  Int_t ngenerated(0);
  
  for ( Int_t j = 0; j < 16; ++j )
  {
    gtm.SetValueAsInt(j,0,0xFFF);
    ++ngenerated;
  }
  return ngenerated;
}

//_____________________________________________________________________________
AliMUONTriggerLut* 
AliMUONCDB::MakeTriggerLUT(const char* file) const
{
  /// Make a triggerlut object, from a file.
  
  AliMUONTriggerLut* lut = new AliMUONTriggerLut;
  lut->ReadFromFile(file);
  return lut;
}

//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells*
AliMUONCDB::MakeTriggerEfficiency(const char* file) const
{
  /// Make a trigger efficiency object from a file.
  
  return new AliMUONTriggerEfficiencyCells(file);
}

//_____________________________________________________________________________
void 
AliMUONCDB::WriteToCDB(const char* calibpath, TObject* object, 
                       Int_t startRun, Int_t endRun, Bool_t defaultValues)
{
  /// Write a given object to OCDB
  
  AliCDBId id(calibpath,startRun,endRun);
  AliCDBMetaData md;
  md.SetAliRootVersion(gROOT->GetVersion());
  if ( defaultValues )
  {
    md.SetComment("Test with default values");
  }
  else
  {
    md.SetComment("Test with random values");
  }
  md.SetResponsible("AliMUONCDB tester class");
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(fCDBPath);
  man->Put(object,id,&md);
}

//_____________________________________________________________________________
Int_t 
AliMUONCDB::MakeNeighbourStore(AliMUONV2DStore& neighbourStore)
{
  /// Fill the neighbours store with, for each channel, a TObjArray of its
  /// neighbouring pads (including itself)
  
  TStopwatch timer;
  
  timer.Start(kTRUE);
  
  TIter next(fManuList);
  
  AliMpIntPair* p;
  
  Int_t nchannels(0);
  
  TObjArray tmp;
  
  while ( ( p = (AliMpIntPair*)next() ) )
  {
    Int_t detElemId = p->GetFirst();
    Int_t manuId = p->GetSecond();
    
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
    
    AliMUONVCalibParam* calibParam = static_cast<AliMUONVCalibParam*>(neighbourStore.Get(detElemId,manuId));
    if (!calibParam)
    {
      Int_t dimension(11);
      Int_t size(64);
      Int_t defaultValue(-1);
      Int_t packingFactor(64);
      
      calibParam = new AliMUONCalibParamNI(dimension,size,defaultValue,packingFactor);
      Bool_t ok = neighbourStore.Set(detElemId,manuId,calibParam,kFALSE);
      if (!ok)
      {
        AliError(Form("Could not set DetElemId=%d manuId=%d",detElemId,manuId));
        return -1;
      }      
    }
    
    for ( Int_t manuChannel = 0; manuChannel < 64; ++manuChannel )
    {
      AliMpPad pad = seg->PadByLocation(AliMpIntPair(manuId,manuChannel),kFALSE);
      
      if (pad.IsValid()) 
      {
        ++nchannels;

        seg->GetNeighbours(pad,tmp,true,true);
        Int_t nofPadNeighbours = tmp.GetEntriesFast();
            
        for ( Int_t i = 0; i < nofPadNeighbours; ++i )
        {
          AliMpPad* pad = static_cast<AliMpPad*>(tmp.At(i));
          Int_t x;
          Bool_t ok = calibParam->PackValues(pad->GetLocation().GetFirst(),pad->GetLocation().GetSecond(),x);
          if (!ok)
          {
            AliError("Could not pack value. Something is seriously wrong. Please check");
            return -1;
          }
          calibParam->SetValueAsInt(manuChannel,i,x);
        }
      }
    }
    }
  
  timer.Print();
  
  return nchannels;
}

//_____________________________________________________________________________
void
AliMUONCDB::WriteLocalTriggerMasks(Int_t startRun, Int_t endRun)
{  
  /// Write local trigger masks to OCDB
  
  AliMUONV1DStore* ltm = new AliMUON1DArray(235);
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
AliMUONCDB::WriteRegionalTriggerMasks(Int_t startRun, Int_t endRun)
{  
  /// Write regional trigger masks to OCDB
  
  AliMUONV1DStore* rtm = new AliMUON1DArray(16);
  Int_t ngenerated = MakeRegionalTriggerMaskStore(*rtm);
  AliInfo(Form("Ngenerated = %d",ngenerated));
  if (ngenerated>0)
  {
    WriteToCDB("MUON/Calib/RegionalTriggerBoardMasks",rtm,startRun,endRun,true);
  }
  delete rtm;
}

//_____________________________________________________________________________
void
AliMUONCDB::WriteGlobalTriggerMasks(Int_t startRun, Int_t endRun)
{  
  /// Write global trigger masks to OCDB
  
  AliMUONVCalibParam* gtm = new AliMUONCalibParamNI(1,16);

  Int_t ngenerated = MakeGlobalTriggerMaskStore(*gtm);
  AliInfo(Form("Ngenerated = %d",ngenerated));
  if (ngenerated>0)
  {
    WriteToCDB("MUON/Calib/GlobalTriggerBoardMasks",gtm,startRun,endRun,true);
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
  
  AliMUONV2DStore* neighbours = new AliMUON2DMap(kTRUE);
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
AliMUONCDB::WritePedestals(Bool_t defaultValues,
                           Int_t startRun, Int_t endRun)
{
  /// generate pedestal values (either 0 if defaultValues=true or random
  /// if defaultValues=false, see makePedestalStore) and
  /// store them into CDB located at cdbpath, with a validity period
  /// ranging from startRun to endRun
  
  AliMUONV2DStore* pedestalStore = new AliMUON2DMap(true);
  Int_t ngenerated = MakePedestalStore(*pedestalStore,defaultValues);
  AliInfo(Form("Ngenerated = %d",ngenerated));
  WriteToCDB("MUON/Calib/Pedestals",pedestalStore,startRun,endRun,defaultValues);
  delete pedestalStore;
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
  
  AliMUONV2DStore* gainStore = new AliMUON2DMap(true);
  Int_t ngenerated = MakeGainStore(*gainStore,defaultValues);
  AliInfo(Form("Ngenerated = %d",ngenerated));  
  WriteToCDB("MUON/Calib/Gains",gainStore,startRun,endRun,defaultValues);
  delete gainStore;
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
  
  AliMUONV1DStore* capaStore = new AliMUON1DMap(16828);
  Int_t ngenerated = MakeCapacitanceStore(*capaStore,defaultValues);
  AliInfo(Form("Ngenerated = %d",ngenerated));
  WriteToCDB("MUON/Calib/Capacitances",capaStore,startRun,endRun,defaultValues);
  delete capaStore;
}

//_____________________________________________________________________________
void
AliMUONCDB::WriteTrigger(Int_t startRun, Int_t endRun)
{
  /// Writes all Trigger related calibration to CDB
  WriteLocalTriggerMasks(startRun,endRun);
  WriteRegionalTriggerMasks(startRun,endRun);
  WriteGlobalTriggerMasks(startRun,endRun);
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
}

