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
#include "AliMUON2DMap.h"
#include "AliMUONV2DStore.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONCalibParam1I.h"
#include "AliMUONCalibParam2F.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
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

//static const char* CDBPath = "local://$ALICE_ROOT/MUON/CDB/Random";
static const char* CDBPath = "local://$ALICE_ROOT/MUON/CDB/Default";

//_____________________________________________________________________________
AliMpSegFactory* segFactory()
{
  static AliMpSegFactory* sf = new AliMpSegFactory();
  return sf;
}

//_____________________________________________________________________________
AliMUONV2DStore* readCDB(const char* calibType)
{
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(CDBPath);

  AliCDBEntry* entry = man->Get(calibType,0);

  if (entry)
    {
      return (AliMUONV2DStore*)entry->GetObject();
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
   
  AliMUONV2DStore* store = readCDB(calibType);
  if (!store) return;
  
  TIter next(manuList());
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
  
  delete store;
}


//_____________________________________________________________________________
void plotCDB(const char* calibType)
{
  TString c(calibType);
  
  TH1* h = 0;
  TH1* h2 = 0;
  
  if ( c == "MUON/Calib/Gains" )
  {
    h = new TH1F("gains_mean","mean gain",100,0,1.5);
    h2 = new TH1F("saturation","adc saturation",4096,-0.5,4095.5);
  }
  else if ( c == "MUON/Calib/Pedestals" )
  {
    h = new TH1F("pedestals_mean","pedestals_mean",4096,-0.5,4095.5);
    h2 = new TH1F("pedestals_sigma","pedestals_sigma",100,0,20);
  }
  else if ( c == "MUON/Calib/DeadChannels" )
  {
    h = new TH1F("dead_channels","dead channels per DE",1500,-0.5,1499.5);
  }
  else
  {
    cerr << "Don't know how to deal with " << calibType << endl;
    return;
  }
  
  AliMUONV2DStore* store = readCDB(calibType);
  if (!store) return;

  TIter next(manuList());
  AliMpIntPair* p;
  Int_t n(0);
  Int_t ndead(0);
  Int_t nPerStation[7];
  
  for ( Int_t i = 0; i < 7; ++i ) nPerStation[i]=0;
  
  while ( ( p = (AliMpIntPair*)next() ) )
  {
    Int_t detElemId = p->GetFirst();
    Int_t manuId = p->GetSecond();
    
    Int_t station = detElemId/100 - 1;
    
    AliMUONVCalibParam* value = 
    dynamic_cast<AliMUONVCalibParam*>(store->Get(detElemId,manuId));
    
    if (value)
    {
      for ( Int_t manuChannel = 0; manuChannel < value->Size(); ++manuChannel )
      {
        ++n;
        ++nPerStation[station];
        if (h2)
        {
          h->Fill(value->ValueAsFloat(manuChannel,0));
          h2->Fill(value->ValueAsFloat(manuChannel,1));
        }
        else
        {
          if( value->ValueAsInt(manuChannel) )
          {
            h->Fill(detElemId);
            ++ndead;
          }
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
  
  if (n && c == "MUON/Calib/DeadChannels")
  {
    cout << "Number of dead channels=" << ndead << endl;
  }
  delete store;
}

//_____________________________________________________________________________
void testReadStore(const AliMUONV2DStore& store, Int_t n)
{
  TIter next(manuList());
  AliMpIntPair* p;
  
  while ( ( p = (AliMpIntPair*)next() ) )
  {
    for ( Int_t i = 0; i < n; ++i )
    {
      store.Get(p->GetFirst(),p->GetSecond());
    }
  }
} 

//_____________________________________________________________________________
Int_t makeStores(AliMUONV2DStore& pedestalStore,
                 AliMUONV2DStore& gainStore,
                 AliMUONV2DStore& deadStore,
                 Bool_t defaultValues)
{  
  TIter next(manuList());
  
  AliMpIntPair* p;
  
  Int_t ngenerated(0);
  
  Bool_t replace = kFALSE;
  
  const Int_t nChannels(64);
  
  while ( ( p = (AliMpIntPair*)next() ) )
  {
    AliMUONVCalibParam* ped = new AliMUONCalibParam2F(nChannels);
    AliMUONVCalibParam* gain = new AliMUONCalibParam2F(nChannels);
    AliMUONVCalibParam* dead = new AliMUONCalibParam1I(nChannels);

    for ( Int_t manuChannel = 0; manuChannel < nChannels; ++manuChannel )
    {
      Float_t meanPedestal;
      Float_t sigmaPedestal;
      Float_t meanGain;
      Float_t saturation(3000);
    
      if ( defaultValues ) 
      {
        meanPedestal = 0.0;
        sigmaPedestal = 1.0;
        meanGain = 1.0;
      }
      else
      {
        meanPedestal = -1;
        while ( meanPedestal < 0 )
        {
          meanPedestal = gRandom->Gaus(150,10);
        }
        sigmaPedestal = -1;
        while ( sigmaPedestal < 0 )
        {
          sigmaPedestal = gRandom->Gaus(1,0.2);
        }
        meanGain = -1;
        while ( meanGain < 0 )
        {
          meanGain = gRandom->Gaus(1,0.05);
        }
      }
      ped->SetValueAsFloat(manuChannel,0,meanPedestal);
      ped->SetValueAsFloat(manuChannel,1,sigmaPedestal);
      gain->SetValueAsFloat(manuChannel,0,meanGain);
      gain->SetValueAsFloat(manuChannel,1,saturation);
      
      if (!defaultValues)
      {
        // probability that this channel is dead ~ 1%
        if ( gRandom->Uniform(100.0) < 1.0 ) 
        {
          Int_t reason = 5; // that value could be used to distinguish
          // why the channel is dead or how it was flagged bad (online,
          // offline, by chance...). 5 is of course a fake number.
          dead->SetValueAsInt(manuChannel,0,reason);
        }
      }
    }
    Int_t detElemId = p->GetFirst();
    Int_t manuId = p->GetSecond();
    Bool_t ok1 = pedestalStore.Set(detElemId,manuId,ped,replace);
    Bool_t ok2 = gainStore.Set(detElemId,manuId,gain,replace);
    Bool_t ok3 = deadStore.Set(detElemId,manuId,dead,replace);
    if (!ok1 || !ok2 || !ok3)
    {
      cout << "Could not set DetElemId=" << detElemId << " manuId="
        << manuId << endl;
    }
    else
    {
      ++ngenerated;
    }
  }
  
  return ngenerated;
}

//_____________________________________________________________________________
TList* manuList(Bool_t reset)
{
  static TList* fgmanuList = new TList;
  
  if (reset) 
  {
    fgmanuList->Delete();
    return fgmanuList;
  }
  
  if (!fgmanuList->IsEmpty()) return fgmanuList;
  
  TStopwatch timer;
  
  cout << "Generating manu list. Please wait" << endl;
  
  fgmanuList->SetOwner(kTRUE);
  
  timer.Start();
  
  AliMpDEIterator it;
  
  it.First();
  
  while ( !it.IsDone() )
  {
    Int_t detElemId = it.CurrentDE();
    AliMpStationType stationType = AliMpDEManager::GetStationType(detElemId);
    if ( stationType != kStationTrigger ) 
    {
      for ( Int_t cath = 0; cath <=1 ; ++cath )
      {
        AliMpVSegmentation* seg = segFactory()->CreateMpSegmentation(detElemId,cath);
        
        TArrayI manus;
        
        seg->GetAllElectronicCardIDs(manus);
        
        for ( Int_t im = 0; im < manus.GetSize(); ++im )
        {
          fgmanuList->Add(new AliMpIntPair(detElemId,manus[im]));
        }        
      }
    }
    it.Next();
  }
  
  cout << "Time to make the manu list = ";
  timer.Print();
  
  return fgmanuList;
}

//_____________________________________________________________________________
void testMakeStores(Int_t readLoop)
{
  manuList();
  
  AliMUONV2DStore* pedestalStore = new AliMUON2DMap;
  AliMUONV2DStore* gainStore = new AliMUON2DMap;
  AliMUONV2DStore* deadStore = new AliMUON2DMap;
  
  TStopwatch timer;
  
  cout << "Creating" << endl;
  
  timer.Start(kTRUE);
  makeStores(*pedestalStore,*gainStore,*deadStore,true);
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
void generateCalibrations(const char* cdbpath, Bool_t defaultValues)
{
  //
  //
  //

  AliMUONV2DStore* pedestalStore = new AliMUON2DMap;
  AliMUONV2DStore* gainStore = new AliMUON2DMap;
  AliMUONV2DStore* deadStore = new AliMUON2DMap; 
  
  makeStores(*pedestalStore,*gainStore,*deadStore,defaultValues);
  
  Int_t ever = 99999999;
  
  AliCDBId id("MUON/Calib/Pedestals",0,ever);
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
  man->Put(pedestalStore,id,&md);
  
  id.SetPath("MUON/Calib/Gains");
  
  man->Put(gainStore,id,(AliCDBMetaData*)md.Clone());
  
  id.SetPath("MUON/Calib/DeadChannels");
  
  man->Put(deadStore,id,(AliCDBMetaData*)md.Clone());
  
  delete deadStore;
  delete pedestalStore;
  delete gainStore;
}


