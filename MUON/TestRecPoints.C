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

/* $Id: TestRecPoints.C 23207 2007-12-20 09:59:20Z ivana $ */

/// \ingroup macros
/// \file TestRecPoints.C
/// \brief This macro is to be used to check the trigger and tracking chamber behaviour
/// through the RecPoints analysis
///
/// \author D. Stocco (Torino)

#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT
#include "TString.h"
#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TDirectory.h"

#include "TTree.h"
#include "TClonesArray.h"

// STEER
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"

// tracker
#include "AliGeomManager.h"

// MUON
#include "AliMpCDB.h"
#include "AliMpDDL.h"
#include "AliMpDDLStore.h"
#include "AliMpConstants.h"

// trigger
#include "AliMUONVTriggerStore.h"
#include "AliMUONDigitStoreV1.h"
#include "AliMUONDigitMaker.h"
#include "AliMUONLocalTrigger.h"

// tracker
#include "AliMUONClusterStoreV2.h"
#include "AliMUONClusterFinderMLEM.h"
#include "AliMUONPreClusterFinder.h"
#include "AliMUONSimpleClusterServer.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONDigitStoreV2R.h"
#include "AliMUONVCluster.h"
#include "AliMUONRecoParam.h"

#endif

const Int_t kNtrigChambers = AliMpConstants::NofTriggerChambers();
const Int_t kNcathodes = AliMpConstants::NofCathodes();

const Int_t kNtrackChambers = AliMpConstants::NofTrackingChambers();

enum {kOnlyTracker, kOnlyTrigger, kTrackTrig};

Int_t GetPlane(Int_t ch, Int_t cath){return kNcathodes * ch + cath;}
void ClusterSize(TList&, AliMUONVDigit*, Int_t&, Int_t);

// Main Method
void TestRecPoints(TString baseDir=".", TString outDir=".", Float_t adcCut = 10., Int_t whatToTest=kTrackTrig, Int_t runNumber=0, TString cdbStorage="local://$ALICE_ROOT/OCDB")
{
  const Int_t kNplanes = kNtrigChambers * kNcathodes;
  const Int_t kNslats = 18;
  
  TList histoListTrig;
  
  TH1F* trigStripMult[kNplanes];
  TH1F* trigClusterSize[kNplanes];
  TH1F* trigClusterMult[kNplanes];
  TString histoName, histoTitle;
  
  if(whatToTest!=kOnlyTracker){
    for(Int_t ch=0; ch<kNtrigChambers; ch++){
      for(Int_t cath=0; cath<kNcathodes; cath++){
        Int_t iplane = GetPlane(ch, cath);
        histoName = "stripMult";
        histoName.Append(Form("Ch%iCath%i",ch,cath));
        histoTitle = Form("Strip multiplicity Ch %i Cathode %i", ch, cath);
        trigStripMult[iplane] = new TH1F(histoName.Data(), histoTitle.Data(), kNslats, 0.-0.5, (Float_t)kNslats - 0.5);
        trigStripMult[iplane]->SetXTitle("slat");
        histoListTrig.Add(trigStripMult[iplane]);
        
        histoName = "clusterSize";
        histoName.Append(Form("Ch%iCath%i",ch,cath));
        histoTitle = Form("Cluster size Ch %i Cathode %i", ch, cath);
        trigClusterSize[iplane] = new TH1F(histoName.Data(), histoTitle.Data(), 10, 0.-0.5, (Float_t)10 - 0.5);
        trigClusterSize[iplane]->SetXTitle("cluster size");
        histoListTrig.Add(trigClusterSize[iplane]);
        
        histoName = "clusterMult";
        histoName.Append(Form("Ch%iCath%i",ch,cath));
        histoTitle = Form("Cluster multiplicity Ch %i Cathode %i", ch, cath);
        trigClusterMult[iplane] = new TH1F(histoName.Data(), histoTitle.Data(), kNslats, 0.-0.5, (Float_t)kNslats - 0.5);
        trigClusterMult[iplane]->SetXTitle("slat");
        histoListTrig.Add(trigClusterMult[iplane]);
      }
    }
  }
  
  const Int_t kMaxDetElem = 26;
  TList histoListTrack;
  TH1F* trackClusterSize[kNtrackChambers];
  TH1F* trackClusterMult[kNtrackChambers];
  TH1F* trackClusterYield[kNtrackChambers];
  
  if(whatToTest!=kOnlyTrigger){
    for(Int_t ch=0; ch<kNtrackChambers; ch++){
      histoName = "clusterSize";
      histoName.Append(Form("Ch%i",ch));
      histoTitle = Form("Cluster size Ch %i", ch);
      trackClusterSize[ch] = new TH1F(histoName.Data(), histoTitle.Data(), kMaxDetElem, 0.-0.5, (Float_t)kMaxDetElem - 0.5);
      trackClusterSize[ch]->SetXTitle("cluster size");
      histoListTrack.Add(trackClusterSize[ch]);
      
      histoName = "clusterMult";
      histoName.Append(Form("Ch%i",ch));
      histoTitle = Form("Cluster multiplicity vs. Charge Ch %i", ch);
      trackClusterMult[ch] = new TH1F(histoName.Data(), histoTitle.Data(), 200, 0, 15000.);
      trackClusterMult[ch]->SetXTitle("Charge (ADC)");
      histoListTrack.Add(trackClusterMult[ch]);
      
      histoName = "clusterYield";
      histoName.Append(Form("Ch%i",ch));
      histoTitle = Form("Cluster yield vs. DetElem Ch %i", ch);
      trackClusterYield[ch] = new TH1F(histoName.Data(), histoTitle.Data(), kMaxDetElem, 0.-0.5, (Float_t)kMaxDetElem - 0.5);
      trackClusterYield[ch]->SetXTitle("Detector element #");
      histoListTrack.Add(trackClusterYield[ch]);
    }
  }
  
  
  // Creating Run Loader and opening RecPoints
  TString filename = baseDir.Data();
  filename.Append("/galice.root");
  AliRunLoader * RunLoader = AliRunLoader::Open(filename.Data(),"MUONLoader","UPDATE");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename.Data());
    return;
  }
  
  // Loading MUON subsystem
  AliLoader* MUONLoader = RunLoader->GetDetectorLoader("MUON");
  if(whatToTest!=kOnlyTracker)  MUONLoader->LoadRecPoints("READ");
  if(whatToTest!=kOnlyTrigger)  MUONLoader->LoadDigits("READ");
  
  TTree *treeR = 0x0, *treeD = 0x0;
  AliMUONVTriggerStore* trigStore = 0x0;
  AliMUONLocalTrigger* locTrg = 0x0;
  
  AliMUONVDigitStore* digitStoreTrack = 0x0;
  AliMUONDigitStoreV2R digitStoreTrackCut;
  AliMUONVCluster* cluster = 0x0;
  
  // Load segmentation
  AliCDBManager::Instance()->SetDefaultStorage(cdbStorage.Data());
  AliCDBManager::Instance()->SetRun(runNumber);
  
  AliMpCDB::LoadDDLStore();
  
  AliMUONGeometryTransformer* transformer = 0x0;
  
  AliMUONRecoParam* recoParam = 0x0;
  
  AliMUONClusterStoreV2 clusterStore;
  AliMUONVClusterFinder* clusterFinder = 0x0;
  
  AliMUONSimpleClusterServer* clusterServer = 0x0;
  
  Int_t firstChamber(0);
  Int_t lastChamber(9);
  
  if(whatToTest!=kOnlyTrigger){
    // Import TGeo geometry
    TString geoFilename = baseDir.Data();
    geoFilename.Append("/geometry.root");
    if ( ! AliGeomManager::GetGeometry() ) {
      AliGeomManager::LoadGeometry(geoFilename);
      if (! AliGeomManager::GetGeometry() ) {
        printf(">>> Error : getting geometry from file %s failed\n", geoFilename.Data());
        return;
      }
    }
    transformer = new AliMUONGeometryTransformer();
    // Load geometry data
    transformer->LoadGeometryData();
    // Load reconstruction parameters
    AliCDBPath path("MUON","Calib","RecoParam");
    AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
    if(entry) {
      recoParam = dynamic_cast<AliMUONRecoParam*>(entry->GetObject());
      entry->SetOwner(0);
      AliCDBManager::Instance()->UnloadFromCache(path.GetPath());
    }
    if (!recoParam) {
      printf("Couldn't find RecoParam object in OCDB: create default one");
      recoParam = AliMUONRecoParam::GetLowFluxParam();
    }
    recoParam->Print("FULL");
    clusterFinder = new AliMUONClusterFinderMLEM(kFALSE,new AliMUONPreClusterFinder);
    clusterServer = new AliMUONSimpleClusterServer(clusterFinder,*transformer);
  }
  
  AliMUONDigitStoreV1 digitStore;
  AliMUONVDigit* mDigit;
  
  Int_t clusterSize;
  
  AliMUONDigitMaker digitMaker;
  TList digitsList[kNplanes];
  
  Int_t nevents = RunLoader->GetNumberOfEvents();
  Int_t analysisFrac = nevents/10+1;
  
  printf("\nNumber of events = %i\n\n",nevents);
  
  for(Int_t ievent=0; ievent<nevents; ievent++){
    if(ievent%analysisFrac==0) printf("Analysing event = %i\n",ievent);
    RunLoader->GetEvent(ievent);
    if(whatToTest!=kOnlyTracker){
      digitStore.Clear();
      treeR = MUONLoader->TreeR();
      trigStore = AliMUONVTriggerStore::Create(*treeR);
      
      if ( trigStore == 0x0 ) continue;
      trigStore->Clear();
      trigStore->Connect(*treeR);
      treeR->GetEvent(0);
      
      TIter nextLocal(trigStore->CreateLocalIterator());
      while ( (locTrg = static_cast<AliMUONLocalTrigger*>( nextLocal() )) != NULL )
      {
        TArrayS xyPattern[2];
        locTrg->GetXPattern(xyPattern[0]);
        locTrg->GetYPattern(xyPattern[1]);
        
        Int_t nBoard = locTrg->LoCircuit();
        digitMaker.TriggerDigits(nBoard, xyPattern, digitStore);
      }
      
      TIter next(digitStore.CreateIterator());
      
      while ( ( mDigit = static_cast<AliMUONVDigit*>(next()) ) )
      {
        Int_t detElemId = mDigit->DetElemId();
        Int_t ch = detElemId/100 - 11;
        Int_t cathode = mDigit->Cathode();
        Int_t slat = detElemId%100;
        Int_t iplane = GetPlane(ch, cathode);
        trigStripMult[iplane]->Fill(slat);
        digitsList[iplane].Add(mDigit);
      } // loop on digits
      
      for(Int_t iplane=0; iplane<kNplanes; iplane++){
        while(digitsList[iplane].GetEntries()){
          clusterSize=1;
          mDigit = (AliMUONVDigit*)digitsList[iplane].At(0);
          digitsList[iplane].Remove(mDigit);
          ClusterSize(digitsList[iplane],mDigit,clusterSize,1);
          //if(clusterSize>1) printf("Cluster size = %i\n\n",clusterSize);
          trigClusterSize[iplane]->Fill(clusterSize);
          
          Int_t detElemId = mDigit->DetElemId();
          Int_t slat = detElemId%100;
          trigClusterMult[iplane]->Fill(slat);
        } // loop o sorted digits
      } // loop on planes
    } // trigger part
    
    if(whatToTest!=kOnlyTrigger){
      clusterStore.Clear();
      treeD = MUONLoader->TreeD();
      digitStoreTrack = AliMUONVDigitStore::Create(*treeD);
      
      if ( digitStoreTrack == 0x0 ) continue;
      digitStoreTrack->Clear();
      digitStoreTrackCut.Clear();
      digitStoreTrack->Connect(*treeD);
      treeD->GetEvent(0);
      
      // Cut low charge channels (pedestal subtraction)
      TIter nextDigitTrack(digitStoreTrack->CreateIterator());
      while ( ( mDigit = static_cast<AliMUONVDigit*>(nextDigitTrack()) ) )
      {
        //printf("Digit class = %s",mDigit->ClassName());
        Float_t charge = mDigit->Charge();
        if(charge<adcCut) continue;
        digitStoreTrackCut.Add(*mDigit, AliMUONVDigitStore::kDeny);
      } // loop on digits
      
      TIter nextDigitTrackCut(digitStoreTrackCut.CreateIterator());
      clusterServer->UseDigits(nextDigitTrackCut);
      
      for (Int_t ch = firstChamber; ch <= lastChamber; ++ch ){
        clusterServer->Clusterize(ch, clusterStore, AliMpArea(),recoParam);
      }
      
      TIter next(clusterStore.CreateIterator());
      while ( ( cluster = static_cast<AliMUONVCluster*>(next()) ) )
      {
        Float_t charge = cluster->GetCharge();
        if(charge==0.) continue;
        Int_t ch = cluster->GetChamberId();
        Int_t npads = cluster->GetNDigits();
        Int_t detElemId = cluster->GetDetElemId()%100;
        //printf("charge = %f   pads = %i   detElemId = %i\n",charge,npads,detElemId);
        
        trackClusterSize[ch]->Fill(npads);
        trackClusterMult[ch]->Fill(charge);
        trackClusterYield[ch]->Fill(detElemId);      
      } // loop on clusters
    } // tracker part
  } // loop on events
  
  MUONLoader->UnloadRecPoints();
  MUONLoader->UnloadDigits();
  RunLoader->UnloadAll();
  delete RunLoader;
  TString outFileName = outDir.Data();
  outFileName.Append("/outTestRecPoints.root");
  TFile* outFile = new TFile(outFileName.Data(), "RECREATE");
  TDirectory* dir = 0x0;
  if(whatToTest!=kOnlyTracker){
    outFile->cd();
    dir = outFile->mkdir("Trigger");
    dir->cd();
    histoListTrig.Write();
  }
  if(whatToTest!=kOnlyTrigger){
    outFile->cd();
    dir = outFile->mkdir("Tracker");
    dir->cd();
    histoListTrack.Write();
  }
  outFile->Close();
  printf("\nSee results in %s\n",outFileName.Data());
}

void ClusterSize(TList& list, AliMUONVDigit* refDigit, Int_t& clusterSize, Int_t depthLevel)
{
  AliMUONVDigit* mDigit = 0x0;
  for(Int_t idigit=0; idigit<list.GetEntries(); idigit++){
    mDigit = (AliMUONVDigit*) list.At(idigit);
    if(mDigit->DetElemId() != refDigit->DetElemId()) continue;
    Int_t diffX = TMath::Abs(mDigit->PadX() - refDigit->PadX());
    Int_t diffY = TMath::Abs(mDigit->PadY() - refDigit->PadY());
    if(diffX>1) continue;
    if(diffY>1) continue;
    if(diffX*diffY != 0) continue;
    clusterSize++;
    list.Remove(mDigit);
    //printf("DetElemId = %i   Level = %i   Ref. (%2i,%2i)    Found (%2i,%2i)\n",mDigit->DetElemId(),depthLevel,refDigit->PadX(),refDigit->PadY(),mDigit->PadX(),mDigit->PadY());
    ClusterSize(list, mDigit, clusterSize, depthLevel+1);
    Int_t val = idigit + depthLevel - clusterSize;
    //printf("Level = %i   Current digit = %i\t",depthLevel,idigit);
    idigit = TMath::Max(-1,val);
    //printf("New digit = %i\n",idigit);
  }
}
