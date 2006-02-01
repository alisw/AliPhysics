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
#include "AliMUON3DMap.h"
#include "AliMUONV3DStore.h"
#include "AliMUONCalibParam.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpPlaneType.h"
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

ClassImp(Triplet)

static const char* CDBPath = "local://$ALICE_ROOT/MUON/CDB/Random";

//_____________________________________________________________________________
AliMpSegFactory* segFactory()
{
  static AliMpSegFactory* sf = new AliMpSegFactory();
  return sf;
}

//_____________________________________________________________________________
AliMUONV3DStore* readCDB(const char* calibType="MUON/Calib/Pedestals")
{
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(CDBPath);

  AliCDBEntry* entry = man->Get(calibType,0);

  if (entry)
    {
      return (AliMUONV3DStore*)entry->GetObject();
    }
  return 0;
}

//_____________________________________________________________________________
void plotCDB(const char* calibType="MUON/Calib/Gains")
{
  TString c(calibType);
  
  TH1* h = 0;
  
  if ( c == "MUON/Calib/Gains" )
  {
    h = new TH1F("gains","gains",100,0,1.5);
  }
  else if ( c == "MUON/Calib/Pedestals" )
  {
    h = new TH1F("pedestals","pedestals",4096,-0.5,4095.5);
  }
  else
  {
    cerr << "Don't know how to deal with " << calibType << endl;
    return;
  }
  
  AliMUONV3DStore* store = readCDB(calibType);
  if (!store) return;

  TIter next(padList());
  Triplet* t;
  
  while ( ( t = (Triplet*)next() ) )
  {
    Int_t detElemId = t->DetElemId();
    Int_t manuId = t->ManuId();
    Int_t manuChannel = t->ManuChannel();
    
    AliMUONCalibParam* value = 
    dynamic_cast<AliMUONCalibParam*>(store->Get(detElemId,manuId,manuChannel));
    if (value)
    {
      h->Fill(value->Mean());
    }
    else
    {
      cout << "Got a null value for DE=" << detElemId << " manuId="
      << manuId << " manuChannel=" << manuChannel << endl;
    }
  }
  
  delete store;
}

//_____________________________________________________________________________
void testReadStore(const AliMUONV3DStore& store, Int_t n)
{
  TIter next(padList());
  Triplet* t;
  
  while ( ( t = (Triplet*)next() ) )
  {
    for ( Int_t i = 0; i < n; ++i )
    {
      store.Get(t->DetElemId(),t->ManuId(),t->ManuChannel());
    }
  }
} 

//_____________________________________________________________________________
Int_t makeStores(AliMUONV3DStore& pedestalStore,
                 AliMUONV3DStore& gainStore,
                 Bool_t defaultValues)
{  
  TIter next(padList());
  
  Triplet* t;
  
  Int_t ngenerated(0);
  
  Bool_t replace = kFALSE;
  
  while ( ( t = (Triplet*)next() ) )
  {
    Float_t meanPedestal;
    Float_t sigmaPedestal;
    Float_t meanGain;
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
    Int_t detElemId = t->DetElemId();
    Int_t manuId = t->ManuId();
    Int_t manuChannel = t->ManuChannel();
    Bool_t ok1 = pedestalStore.Set(detElemId,manuId,manuChannel,
                                   new AliMUONCalibParam(meanPedestal,sigmaPedestal),replace);
    Bool_t ok2 = gainStore.Set(detElemId,manuId,manuChannel,
                               new AliMUONCalibParam(meanGain,0),replace);
    if (!ok1 || !ok2)
    {
      cout << "Could not set DetElemId=" << detElemId << " manuId="
      << manuId << " manuChannel=" << manuChannel << endl;
    }
    else
    {
      ++ngenerated;
    }
  }
  
  return ngenerated;
}

//_____________________________________________________________________________
TList* padList(Bool_t reset)
{
  static TList* fgPadList = new TList;
  
  if (reset) 
  {
    fgPadList->Delete();
    return fgPadList;
  }
  
  if (!fgPadList->IsEmpty()) return fgPadList;
  
  TStopwatch timer;
  
  cout << "Generating pad list. Please wait" << endl;
  
  fgPadList->SetOwner(kTRUE);
  
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
        AliMpPlaneType planeType = AliMpDEManager::GetPlaneType(detElemId,cath);
        
        for ( Int_t ix = 0; ix <= seg->MaxPadIndexX(); ++ix )
        { 
          for ( Int_t iy = 0; iy <= seg->MaxPadIndexY(); ++iy )
          {
            AliMpPad pad = seg->PadByIndices(AliMpIntPair(ix,iy),kFALSE);
            if ( pad.IsValid() )
            {
              Int_t manuId = pad.GetLocation().GetFirst();
              if ( planeType == kNonBendingPlane )
              {
                manuId |= (1<<11);
              }
              Int_t manuChannel = pad.GetLocation().GetSecond();
              fgPadList->Add(new Triplet(detElemId,manuId,manuChannel));
            }
          }
        }
      }
    }
    it.Next();
  }
  
  cout << "Time to make the pad list = ";
  timer.Print();
  
  return fgPadList;
}

//_____________________________________________________________________________
void testMakeStores(Int_t readLoop=10)
{
  padList();
  
  AliMUONV3DStore* pedestalStore = new AliMUON3DMap;
  AliMUONV3DStore* gainStore = new AliMUON3DMap;
  
  TStopwatch timer;
  
  cout << "Creating" << endl;
  
  timer.Start(kTRUE);
  makeStores(*pedestalStore,*gainStore,true);
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
void generateCalibrations(const char* cdbpath = CDBPath,
                          Bool_t defaultValues = kFALSE)
{
  //
  //
  //

  AliMUONV3DStore* pedestalStore = new AliMUON3DMap;
  AliMUONV3DStore* gainStore = new AliMUON3DMap;
 
  makeStores(*pedestalStore,*gainStore,defaultValues);
  
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
  
  delete pedestalStore;
  delete gainStore;
}


