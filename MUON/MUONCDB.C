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

/// By Laurent Aphecetche

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "MUONCDB.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliDCSValue.h"
#include "AliMUON1DArray.h"
#include "AliMUON2DMap.h"
#include "AliMUON2DStoreValidator.h"
#include "AliMUONCalibParam1I.h"
#include "AliMUONCalibParam2F.h"
#include "AliMUONConstants.h"
#include "AliMUONHVNamer.h"
#include "AliMUONObjectPair.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerLut.h"
#include "AliMUONV2DStore.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVDataIterator.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpManuList.h"
#include "AliMpSegmentation.h"
#include "AliMpStationType.h"
#include "AliMpVSegmentation.h"
#include "Riostream.h"
#include "TArrayI.h"
#include "TH1F.h"
#include "TList.h"
#include "TMap.h"
#include "TObjString.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include <map>

#endif

//_____________________________________________________________________________
AliMUONV2DStore* diff(AliMUONV2DStore& store1, AliMUONV2DStore& store2, 
                      const char* /* opt */)
{
  // creates a store which contains store1-store2
  // if opt="abs" the difference is absolute one,
  // if opt="rel" then what is stored is (store1-store2)/store1
  
  AliMUONV2DStore* d = static_cast<AliMUONV2DStore*>(store1.Clone());

  AliMUONVDataIterator* it = d->Iterator();
  
  AliMUONObjectPair* p;
  
  while ( ( p = dynamic_cast<AliMUONObjectPair*>(it->Next() ) ) )
  {
    AliMpIntPair* dm = dynamic_cast<AliMpIntPair*>(p->Key());
    //FIXME: this might happen (if a full manu is missing, for instance)
    //handle it.
    assert(dm!=0);
    AliMUONVCalibParam* param = dynamic_cast<AliMUONVCalibParam*>(p->Value());
    //FIXMENOT: this should *not* happen
    assert(param!=0);

    AliMUONVCalibParam* param2 = dynamic_cast<AliMUONVCalibParam*>(store2.Get(dm->GetFirst(),dm->GetSecond()));
    //FIXME: this might happen. Handle it.
    assert(param2!=0);
    
    for ( Int_t i = 0; i < param->Size(); ++i )
    {
      for ( Int_t j = 0; j < param->Dimension(); ++j )
      {
        param->SetValueAsFloat(i,j,param->ValueAsFloat(i,j) - param2->ValueAsFloat(i,j));
      }      
    }
  }
  return d;
}

//_____________________________________________________________________________
void getBoundaries(const AliMUONV2DStore& store,
                   Float_t& x0min, Float_t& x0max,
                   Float_t& x1min, Float_t& x1max)
{
  x0min=1E30;
  x0max=-1E30;
  x1min=1E30;
  x1max=-1E30;
  
  TList* list = AliMpManuList::ManuList();
  TIter next(list);
  AliMpIntPair* p;
  
  while ( ( p = (AliMpIntPair*)next() ) )
  {
    Int_t detElemId = p->GetFirst();
    Int_t manuId = p->GetSecond();
        
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
          
    AliMUONVCalibParam* value = 
      dynamic_cast<AliMUONVCalibParam*>(store.Get(detElemId,manuId));
    
    if (!value) continue;

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
  }  
  delete list;
}

//_____________________________________________________________________________
void plot(const AliMUONV2DStore& store, const char* name, Int_t nbins)
{
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
  
  TList* list = AliMpManuList::ManuList();
  TIter next(list);
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
        h0->Fill(x);
        if (h1)
        {
          h1->Fill(value->ValueAsFloat(manuChannel,1));
        }
      }
    }
    else
    {
      cout << "Got a null value for DE=" << detElemId << " manuId="
      << manuId << endl;
    }
  }
  
  cout << "Number of channels = " << n << endl;
  for ( Int_t i = 0; i < 7; ++i )
  {
    cout << "Station " << (i+1) << " " << nPerStation[i] << endl;
  }
  
  delete list;
}

//_____________________________________________________________________________
void testReadStore(const AliMUONV2DStore& store, Int_t n)
{
  TList* list = AliMpManuList::ManuList();
  TIter next(list);
  AliMpIntPair* p;
  
  while ( ( p = (AliMpIntPair*)next() ) )
  {
    for ( Int_t i = 0; i < n; ++i )
    {
      store.Get(p->GetFirst(),p->GetSecond());
    }
  }
  delete list;
} 

//_____________________________________________________________________________
Int_t makeHVStore(TMap& aliasMap, Bool_t defaultValues)
{
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
  
  cout << nChannels << " HV channels and " << nSwitch << " switches" << endl;
  
  return nChannels+nSwitch;
}

//_____________________________________________________________________________
Int_t makePedestalStore(AliMUONV2DStore& pedestalStore, Bool_t defaultValues)
{
  TList* list = AliMpManuList::ManuList();
  TIter next(list);
  
  AliMpIntPair* p;
  
  Int_t nchannels(0);
  Int_t nmanus(0);
  
  Bool_t replace = kFALSE;
  
  const Int_t nChannels(64);
  const Float_t kPedestalMeanMean(150);
  const Float_t kPedestalMeanSigma(10);
  const Float_t kPedestalSigmaMean(1.0);
  const Float_t kPedestalSigmaSigma(0.2);
  
  while ( ( p = (AliMpIntPair*)next() ) )
  {
    ++nmanus;
    AliMUONVCalibParam* ped = new AliMUONCalibParam2F(nChannels,AliMUONVCalibParam::InvalidFloatValue());
    
    Int_t detElemId = p->GetFirst();
        
    Int_t manuId = p->GetSecond();
    
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
    
    for ( Int_t manuChannel = 0; manuChannel < nChannels; ++manuChannel )
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
      cout << "Could not set DetElemId=" << detElemId << " manuId="
      << manuId << endl;
    }
  }
  
  delete list;
  cout << nmanus << " Manus and " << nchannels << " channels." << endl;
  return nchannels;
  
}

//_____________________________________________________________________________
Int_t makeGainStore(AliMUONV2DStore& gainStore, Bool_t defaultValues)
{  
  TList* list = AliMpManuList::ManuList();
  TIter next(list);
  
  AliMpIntPair* p;
  
  Int_t nchannels(0);
  Int_t nmanus(0);
  
  Bool_t replace = kFALSE;
  
  const Int_t nChannels(64);
  const Double_t kSaturation(3000);
  const Double_t kGainMean(1.0);
  const Double_t kGainSigma(0.05);
  
  while ( ( p = (AliMpIntPair*)next() ) )
  {
    ++nmanus;
    AliMUONVCalibParam* gain = new AliMUONCalibParam2F(nChannels,AliMUONVCalibParam::InvalidFloatValue());

    Int_t detElemId = p->GetFirst();
    Int_t manuId = p->GetSecond();

    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);

    for ( Int_t manuChannel = 0; manuChannel < nChannels; ++manuChannel )
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
      cout << "Could not set DetElemId=" << detElemId << " manuId="
        << manuId << endl;
    }
  }
  
  delete list;
  cout << nmanus << " Manus and " << nchannels << " channels." << endl;
  return nchannels;
}

//_____________________________________________________________________________
void testMakeStores(Int_t readLoop)
{
  AliMUONV2DStore* pedestalStore = new AliMUON2DMap;
  AliMUONV2DStore* gainStore = new AliMUON2DMap;
  
  TStopwatch timer;
  
  cout << "Creating" << endl;
  
  Bool_t defaultValues = kTRUE;
  
  timer.Start(kTRUE);
  makePedestalStore(*pedestalStore,defaultValues);
  makeGainStore(*gainStore,defaultValues);
  timer.Print();
  
  cout << "Reading..." << endl;
  timer.Start(kTRUE);
  testReadStore(*pedestalStore,readLoop);
  testReadStore(*gainStore,readLoop);
  cout << timer.CpuTime()/readLoop << " CPUs (mean of " << readLoop 
    <<" samples." << endl;
  
  delete pedestalStore;
  delete gainStore;
}

//_____________________________________________________________________________
void generateTrigger(const char* cdbpath)
{
  // 
  // Generate trigger related conditions :
  //
  // - trigger masks for board (locals, regionals, global)
  // - trigger lut
  // - trigger efficiency
  // - trigger switches (to be implemented FIXME)
  //
  const Int_t nlboards = 234;
  AliMUONV1DStore* localBoardMasks = new AliMUON1DArray(nlboards+1);
  
  // Generate fake mask values for 234 localboards and put that into
  // one single container (localBoardMasks)
  for ( Int_t i = 1; i <= nlboards; ++i )
  {
    AliMUONVCalibParam* localBoard = new AliMUONCalibParam1I(8);
    for ( Int_t x = 0; x < 2; ++x )
    {
      for ( Int_t y = 0; y < 4; ++y )
      {
        Int_t index = x*4+y;
        localBoard->SetValueAsInt(index,0,0xFFFF);
      }
    }
    localBoardMasks->Set(i,localBoard,kFALSE);
  }
  
  // Generate values for regional boards
  const Int_t nrboards = 16;
  AliMUONV1DStore* regionalBoardMasks = new AliMUON1DArray(16);
  
  for ( Int_t i = 0; i < nrboards; ++i )
  {
    AliMUONVCalibParam* regionalBoard = new AliMUONCalibParam1I(16);
    for ( Int_t j = 0; j < 16; ++j )
    {
      regionalBoard->SetValueAsInt(j,0,0x3F);
    }
    regionalBoardMasks->Set(i,regionalBoard,kFALSE);
  }
  
  // Generate values for global board
  AliMUONVCalibParam* globalBoardMasks = new AliMUONCalibParam1I(16);
  
  for ( Int_t j = 0; j < 16; ++j )
  {
    globalBoardMasks->SetValueAsInt(j,0,0xFFF);
  }
    
  AliMUONTriggerLut lut;
  lut.ReadFromFile("$(ALICE_ROOT)/MUON/data/lutAptLpt1Hpt1p7.root");
  
  AliMUONTriggerEfficiencyCells cells("$ALICE_ROOT/MUON/data/efficiencyCells.dat");
  
  //--------------------------------------------
  // Store the resulting containers into the CDB
  Int_t ever = 99999999;
  
  AliCDBId id("MUON/Calib/LocalTriggerBoardMasks",0,ever);
  
  AliCDBMetaData md;
  md.SetBeamPeriod(1);
  md.SetAliRootVersion(gROOT->GetVersion());
  md.SetComment("Test with default values");
  md.SetResponsible("Rachid Guernane");
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbpath);
  man->Put(localBoardMasks,id,&md);
  
  id.SetPath("MUON/Calib/RegionalTriggerBoardMasks");
  
  man->Put(regionalBoardMasks,id,(AliCDBMetaData*)md.Clone());
  
  id.SetPath("MUON/Calib/GlobalTriggerBoardMasks");
  
  man->Put(globalBoardMasks,id,(AliCDBMetaData*)md.Clone());

  id.SetPath("MUON/Calib/TriggerLut");
  
  man->Put(&lut,id,(AliCDBMetaData*)md.Clone());
  
  id.SetPath("MUON/Calib/TriggerEfficiency");
  md.SetResponsible("Diego Stocco");
  
  man->Put(&cells,id,(AliCDBMetaData*)md.Clone());
  
  delete localBoardMasks;
  delete regionalBoardMasks;
  delete globalBoardMasks;
}


//_____________________________________________________________________________
void writeToCDB(const char* cdbpath, const char* calibpath, TObject* object, 
                Int_t startRun, Int_t endRun, Bool_t defaultValues)
{
  AliCDBId id(calibpath,startRun,endRun);
  AliCDBMetaData md;
  md.SetBeamPeriod(1);
  md.SetAliRootVersion(gROOT->GetVersion());
  if ( defaultValues )
  {
    md.SetComment("Test with default values");
  }
  else
  {
    md.SetComment("Test with random values");
  }
  md.SetResponsible("Laurent Aphecetche");
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbpath);
  man->Put(object,id,&md);
}

//_____________________________________________________________________________
void writeHV(const char* cdbpath, Bool_t defaultValues,
             Int_t startRun, Int_t endRun)
{
  /// generate HV values (either cste = 1500 V) if defaultValues=true or random
  /// if defaultValues=false, see makeHVStore) and
  /// store them into CDB located at cdbpath, with a validity period
  /// ranging from startRun to endRun
  
  TMap* hvStore = new TMap;
  Int_t ngenerated = makeHVStore(*hvStore,defaultValues);
  cout << "Ngenerated = " << ngenerated << endl;
  
  writeToCDB(cdbpath,"MUON/Calib/HV",hvStore,startRun,endRun,defaultValues);
  
  delete hvStore;
}

//_____________________________________________________________________________
void writePedestals(const char* cdbpath, Bool_t defaultValues,
                    Int_t startRun, Int_t endRun)
{
  /// generate pedestal values (either 0 if defaultValues=true or random
  /// if defaultValues=false, see makePedestalStore) and
  /// store them into CDB located at cdbpath, with a validity period
  /// ranging from startRun to endRun
  
  AliMUONV2DStore* pedestalStore = new AliMUON2DMap;
  Int_t ngenerated = makePedestalStore(*pedestalStore,defaultValues);
  cout << "Ngenerated = " << ngenerated << endl;
  
  writeToCDB(cdbpath,"MUON/Calib/Pedestals",pedestalStore,startRun,endRun,defaultValues);
  delete pedestalStore;
}

//_____________________________________________________________________________
void writeGains(const char* cdbpath, Bool_t defaultValues,
                    Int_t startRun, Int_t endRun)
{
  /// generate gain values (either 1 if defaultValues=true or random
  /// if defaultValues=false, see makePedestalStore) and
  /// store them into CDB located at cdbpath, with a validity period
  /// ranging from startRun to endRun
  
  AliMUONV2DStore* gainStore = new AliMUON2DMap;
  Int_t ngenerated = makeGainStore(*gainStore,defaultValues);
  cout << "Ngenerated = " << ngenerated << endl;
  
  writeToCDB(cdbpath,"MUON/Calib/Gains",gainStore,startRun,endRun,defaultValues);
  delete gainStore;
}

//_____________________________________________________________________________
Bool_t check1I(const AliMUONVCalibParam& calib, Int_t channel)
{
  /// 
  
  return ( calib.ValueAsInt(channel) == 0);
}

//_____________________________________________________________________________
void validate1I(const AliMUONV2DStore& store)
{
  AliMUON2DStoreValidator validator;
  
  TObjArray* a = validator.Validate(store,check1I);
  
  if (a) a->Print();
  
  TList lines;
  
  validator.Report(lines);
  
  lines.Print();
}

//_____________________________________________________________________________
void validate2F(const AliMUONV2DStore& store)
{
  AliMUON2DStoreValidator validator;
  
  TObjArray* a = validator.Validate(store,AliMUONCalibParam2F::InvalidFloatValue());
  
  if (a) a->Print();
  
  TList lines;
  
  validator.Report(lines);
  
  lines.Print();
}

//_____________________________________________________________________________
void dump(const TArrayI& a, const char* what)
{
  cout << what << " " << a.GetSize() << " manus" << endl;
  for ( Int_t i = 0; i < a.GetSize(); ++i ) 
  {
    cout << Form(" %5d ",a[i]);
  }
  cout << endl;
}

//_____________________________________________________________________________
void countManus()
{
  AliMpDEIterator it;
  
  it.First();
  
  Int_t b(0);
  Int_t nb(0);
  
  while (!it.IsDone())
  {
    Int_t detElemId = it.CurrentDEId();
    AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
    if ( stationType != AliMp::kStationTrigger ) 
    {
      const AliMpVSegmentation* seg0 = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::kCath0);
      const AliMpVSegmentation* seg1 = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::kCath1);
      TArrayI a1;
      TArrayI a0;
      seg0->GetAllElectronicCardIDs(a0);
      seg1->GetAllElectronicCardIDs(a1);
      
      cout << Form("DE %5d B %4d NB %4d Total %5d",detElemId,
                   a0.GetSize(),a1.GetSize(),a0.GetSize()+a1.GetSize())            
        << endl;
      
      b += a0.GetSize();
      nb += a1.GetSize();
      
      if ( detElemId == 500 ) 
      {
        dump(a0,"B");
        dump(a1,"NB");
      }
    }
    it.Next();
  }
  
  cout << Form("B %5d NB %5d Total %5d",b,nb,b+nb) << endl;
}

//_____________________________________________________________________________
void count(const AliMUONV2DStore& store)
{
  AliMUONVDataIterator* it = store.Iterator();
  AliMUONObjectPair* pair;
  std::map<int,std::pair<int,int> > demap;
  
  while ( ( pair = static_cast<AliMUONObjectPair*>(it->Next()) ) )
  {
    AliMpIntPair* ip = static_cast<AliMpIntPair*>(pair->First());
    
    Int_t detElemId = ip->GetFirst();
    
    Int_t manuId = ip->GetSecond();
    
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId);
    
    if ( seg->PlaneType() == AliMp::kNonBendingPlane ) 
    {
      demap[detElemId].second++;
    }
    else
    {
      demap[detElemId].first++;
    }    
  }
  
  std::map<int,std::pair<int,int> >::const_iterator mit;
  
  Int_t b(0);
  Int_t nb(0);
  
  for ( mit = demap.begin(); mit != demap.end(); ++mit ) 
  {
    cout << Form("DE %5d B %4d NB %4d Total %5d",mit->first,
                 mit->second.first,mit->second.second,
                 mit->second.first+mit->second.second) << endl;
    b += mit->second.first;
    nb += mit->second.second;    
  }
  
  cout << Form("B %5d NB %5d Total %5d",b,nb,b+nb) << endl;  
}


