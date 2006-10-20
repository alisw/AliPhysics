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
#include "AliMUON1DArray.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParam1I.h"
#include "AliMUONCalibParam2F.h"
#include "AliMUONConstants.h"
#include "AliMUONObjectPair.h"
#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONTriggerLut.h"
#include "AliMUONV2DStore.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVDataIterator.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpManuList.h"
#include "AliMpSegFactory.h"
#include "AliMpStationType.h"
#include "AliMpVSegmentation.h"
#include "Riostream.h"
#include "TH1F.h"
#include "TList.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TSystem.h"
#endif

//_____________________________________________________________________________
Int_t countChannels(AliMpVSegmentation& seg)
{
  Int_t n(0);
  
  for ( Int_t ix = 0; ix < seg.MaxPadIndexX(); ++ix )
  {
    for ( Int_t iy = 0; iy < seg.MaxPadIndexY(); ++iy )
    {
      if ( seg.HasPad(AliMpIntPair(ix,iy)) ) ++n;
    }
  }
  return n;
}

//_____________________________________________________________________________
AliMpSegFactory* segFactory()
{
  static AliMpSegFactory* sf = new AliMpSegFactory();
  return sf;
}

//_____________________________________________________________________________
void countChannels()
{
  AliMpDEIterator it;
  Int_t ntotal(0);
  Int_t ntracker(0);
  Int_t ntrigger(0);
  
  for ( Int_t station = 0; station < AliMUONConstants::NCh(); ++station )
  {
    Int_t n(0);
    it.First(station);
    while (!it.IsDone())
    {
      Int_t de = it.CurrentDE();
      for ( Int_t cathode = 0; cathode < 2; ++cathode )
      {
        AliMpVSegmentation* seg = segFactory()->CreateMpSegmentation(de,cathode);
        n += countChannels(*seg);
      }
      it.Next();
    }
    cout << "Station " << station << " has " << n << " channels" << endl;
    if ( station < AliMUONConstants::NTrackingCh() )
    {
      ntracker += n;
    }
    else
    {
      ntrigger += n;
    }
    ntotal += n;
  }
  cout << "Tracker channels = " << ntracker << endl;
  cout << "Trigger channels = " << ntrigger << endl;
  cout << "Total channels =" << ntotal << endl;
}

//_____________________________________________________________________________
AliMUONV2DStore* read2D(const char* calibType, Int_t runNumber)
{
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(CDBPath);

  AliCDBEntry* entry = man->Get(calibType,runNumber);

  if (entry)
    {
      return (AliMUONV2DStore*)entry->GetObject();
    }
  return 0;
}

//_____________________________________________________________________________
AliMUONV1DStore* read1D(const char* calibType, Int_t runNumber)
{
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(CDBPath);
  
  AliCDBEntry* entry = man->Get(calibType,runNumber);
  
  if (entry)
  {
    return (AliMUONV1DStore*)entry->GetObject();
  }
  return 0;
}

//_____________________________________________________________________________
void checkCDB(const char* calibType)
{
  TString c(calibType);
  Float_t refValue(0);
  
  if ( c == "MUON/Calib/DeadChannels" )
  {
    refValue=5;
  }
   
  AliMUONV2DStore* store = read2D(calibType);
  if (!store) return;
  
  TList* list = AliMpManuList::ManuList();
  TIter next(list);
  AliMpIntPair* p;
  
  while ( ( p = (AliMpIntPair*)next() ) )
  {
    Int_t detElemId = p->GetFirst();
    Int_t manuId = p->GetSecond();
    
    AliMUONVCalibParam* value = 
      dynamic_cast<AliMUONVCalibParam*>(store->Get(detElemId,manuId));
    
    if (value)
    {
      for ( Int_t manuChannel = 0; manuChannel < value->Size(); ++manuChannel )
      {
        Float_t testValue = value->ValueAsFloat(manuChannel,0);
        if ( testValue && testValue != refValue )
        {
          cout << "Got a strange value for DE=" << detElemId << " manuId="
          << manuId << " manuChannel=" << manuChannel << " was expecting "
          << refValue << " and I got " << testValue << endl;
        }
      }
    }
    else
    {
      cout << "Got a null value for DE=" << detElemId << " manuId="
      << manuId << endl;
    }
  }
  
  delete list;
  delete store;
}

//_____________________________________________________________________________
//void testDump(AliMUONV2DStore& store, int n)
//{  
//  AliMUONObjectPair* p;
//  
//  Int_t c(0);
//  
//  for ( Int_t i = 0; i < n; ++i )
//  {
//    AliMUONVDataIterator* it = store.Iterator();
//    
//    while ( ( p = dynamic_cast<AliMUONObjectPair*>(it->Next() ) ) )
//    {
//      AliMpIntPair* dm = dynamic_cast<AliMpIntPair*>(p->Key());
//      if (dm)
//      {
//        Int_t a(dm->GetFirst()+dm->GetSecond());
//        a=2;
//        ++c;
//        AliMUONVCalibParam* param = dynamic_cast<AliMUONVCalibParam*>(p->Value());
//      }
//      delete p;
//    }
//    delete it;
//  }
//    
//  cout << c << endl;
//}

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
        
    AliMpVSegmentation* seg = 
      segFactory()->CreateMpSegmentationByElectronics(detElemId,manuId);
          
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
    Int_t station = AliMpDEmanager::GetChamberId(detElemId);
    
    AliMpVSegmentation* seg = 
      segFactory()->CreateMpSegmentationByElectronics(detElemId,manuId);
    
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
void plotCDB(const char* calibType, Int_t runNumber)
{
  AliMUONV2DStore* store = read2D(calibType,runNumber);
  if (!store) return;

  TString c(calibType);
  c.ReplaceAll("/","_");
  
  plot(*store,c.Data());
  
  delete store;
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
    
    if ( detElemId / 100 == 1 || detElemId / 1000 == 1 ) continue;
    
    if ( detElemId == 501 || detElemId == 903 ) continue;
    
    Int_t manuId = p->GetSecond();
    
    if ( manuId == 4 && detElemId / 700 > 0 ) continue;
    
    AliMpVSegmentation* seg = 
      segFactory()->CreateMpSegmentationByElectronics(detElemId,manuId);
    
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

    AliMpVSegmentation* seg = 
      segFactory()->CreateMpSegmentationByElectronics(detElemId,manuId);

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
Int_t makeDeadStore(AliMUONV2DStore& deadStore, Bool_t defaultValues)
{  
  TList* list = AliMpManuList::ManuList();
  TIter next(list);
  
  AliMpIntPair* p;
  
  Int_t nchannels(0);
  Int_t nmanus(0);
  
  Bool_t replace = kFALSE;
  
  const Int_t nChannels(64);
  const Double_t deadProba = 1.0; // 1%
  
  while ( ( p = (AliMpIntPair*)next() ) )
  {
    ++nmanus;
    AliMUONVCalibParam* dead = new AliMUONCalibParam1I(nChannels,-9999);
    
    Int_t detElemId = p->GetFirst();
    Int_t manuId = p->GetSecond();
    
    AliMpVSegmentation* seg = 
      segFactory()->CreateMpSegmentationByElectronics(detElemId,manuId);
    
    for ( Int_t manuChannel = 0; manuChannel < nChannels; ++manuChannel )
    {
      AliMpPad pad = seg->PadByLocation(AliMpIntPair(manuId,manuChannel),kFALSE);
      if (!pad.IsValid()) continue;
      
      ++nchannels;
            
      if (!defaultValues)
      {
        // probability that this channel is dead ~ 1%
        if ( gRandom->Uniform(100.0) < deadProba ) 
        {
          Int_t reason = 5; // that value could be used to distinguish
                            // why the channel is dead or how it was flagged bad (online,
                            // offline, by chance...). 5 is of course a fake number.
          dead->SetValueAsInt(manuChannel,0,reason);
        }
      }
    }
    Bool_t ok = deadStore.Set(detElemId,manuId,dead,replace);
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
  AliMUONV2DStore* deadStore = new AliMUON2DMap;
  
  TStopwatch timer;
  
  cout << "Creating" << endl;
  
  Bool_t defaultValues = kTRUE;
  
  timer.Start(kTRUE);
  makePedestalStore(*pedestalStore,defaultValues);
  makeGainStore(*gainStore,defaultValues);
  makeDeadStore(*deadStore,defaultValues);
  timer.Print();
  
  cout << "Reading..." << endl;
  timer.Start(kTRUE);
  testReadStore(*pedestalStore,readLoop);
  testReadStore(*gainStore,readLoop);
  testReadStore(*deadStore,readLoop);
  cout << timer.CpuTime()/readLoop << " CPUs (mean of " << readLoop 
    <<" samples." << endl;
  
  delete pedestalStore;
  delete gainStore;
  delete deadStore;
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
void writePedestals(const char* cdbpath, Bool_t defaultValues,
                    Int_t startRun, Int_t endRun)
{
  AliMUONV2DStore* pedestalStore = new AliMUON2DMap;
  Int_t ngenerated = makePedestalStore(*pedestalStore,defaultValues);
  cout << "Ngenerated = " << ngenerated << endl;
  
  writeToCDB(cdbpath,"MUON/Calib/Pedestals",pedestalStore,startRun,endRun,defaultValues);
}

//_____________________________________________________________________________
void writeGains(const char* cdbpath, Bool_t defaultValues,
                    Int_t startRun, Int_t endRun)
{
  AliMUONV2DStore* gainStore = new AliMUON2DMap;
  Int_t ngenerated = makeGainStore(*gainStore,defaultValues);
  cout << "Ngenerated = " << ngenerated << endl;
  
  writeToCDB(cdbpath,"MUON/Calib/Gains",gainStore,startRun,endRun,defaultValues);
}

//_____________________________________________________________________________
void writeDeadChannels(const char* cdbpath, Bool_t defaultValues,
                       Int_t startRun, Int_t endRun)
{
  AliMUONV2DStore* deadStore = new AliMUON2DMap;
  Int_t ngenerated = makeDeadStore(*deadStore,defaultValues);
  cout << "Ngenerated = " << ngenerated << endl;
  
  writeToCDB(cdbpath,"MUON/Calib/DeadChannels",deadStore,startRun,endRun,defaultValues);
}


