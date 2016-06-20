// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel & Bogdan Vulpescu: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <string.h>

#include "AliEveMUONData.h"

#include <AliEveMUONChamberData.h>
#include <AliEveEventManager.h>

#include <AliRawReader.h>
#include <AliRawReaderFile.h>
#include <AliRawReaderDate.h>
#include <AliRawReaderRoot.h>

#include <AliMUONDigitMaker.h>
#include <AliMUONHit.h>
#include <AliMUONVCluster.h>
#include "AliMUONVClusterStore.h"
#include <AliMUONVDigit.h>
#include "AliMUONDigitStoreV1.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrack.h"
#include "AliMUONESDInterface.h"
#include "AliESDMuonTrack.h"
#include "AliESDEvent.h"
#include "TTree.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TFile.h"
#include <iostream>


//______________________________________________________________________________
// AliEveMUONData
//

using std::cout;
using std::endl;
ClassImp(AliEveMUONData)

AliRawReader*            AliEveMUONData::fgRawReader        = 0;

//______________________________________________________________________________
AliEveMUONData::AliEveMUONData() :
  fChambers(14),
  fNTrackList(0)
{
  //
  // Constructor
  //

  for (Int_t i = 0; i < 256; i++) {
    fTrackList[i] = -1;
  }

  CreateAllChambers();

}

//______________________________________________________________________________
AliEveMUONData::~AliEveMUONData()
{
  //
  // Destructor
  //

  DeleteAllChambers();

}

//______________________________________________________________________________
void AliEveMUONData::Reset()
{
  //
  // Reset data
  //

  //DropAllChambers();

  fNTrackList = 0;
  for (Int_t i = 0; i < 256; i++) {
    fTrackList[i] = -1;
  }

}

//______________________________________________________________________________
AliEveMUONData::AliEveMUONData(const AliEveMUONData &mdata) :
  TObject(mdata),
  TEveRefCnt(),
  fChambers(14),
  fNTrackList(0)
{
  //
  // Copy constructor
  //
  memset(fTrackList,0,256*sizeof(Int_t));
}

//______________________________________________________________________________
AliEveMUONData& AliEveMUONData::operator=(const AliEveMUONData &mdata)
{
  //
  // Assignment operator
  //

  if (this != &mdata) {

  }

  return *this;

}

//______________________________________________________________________________
void AliEveMUONData::CreateChamber(Int_t chamber)
{
  //
  // create data for the chamber with id=chamber (0 to 13)
  //

  if (fChambers[chamber] == 0)
    fChambers[chamber] = new AliEveMUONChamberData(chamber);

}

//______________________________________________________________________________
void AliEveMUONData::CreateAllChambers()
{
  //
  // create all 14 chambers data
  //

  for (Int_t c = 0; c < 14; ++c)
    CreateChamber(c);

}

//______________________________________________________________________________
void AliEveMUONData::DropAllChambers()
{
  //
  // release data from all chambers
  //

  for (Int_t c = 0; c < 14; ++c) {

    if (fChambers[c] != 0)
      fChambers[c]->DropData();

  }

}

//______________________________________________________________________________
void AliEveMUONData::DeleteAllChambers()
{
  //
  // delete all chambers data
  //

  for (Int_t c = 0; c < 14; ++c) {

    delete fChambers[c];
    fChambers[c] = 0;

  }

}

//______________________________________________________________________________
void AliEveMUONData::RegisterTrack(Int_t track)
{
  //
  // register (in a list) a track with hits in the chambers
  //

  if (fNTrackList == (256-1)) {
    cout << "Maximum of registered tracks reached..." << endl;
    return;
  }

  Bool_t inList = kFALSE;
  for (Int_t i = 0; i < fNTrackList; i++) {
    if (track == fTrackList[i]) {
      inList = kTRUE;
      break;
    }
  }
  if (!inList) {
    fTrackList[fNTrackList] = track;
    fNTrackList++;
  }

}

//______________________________________________________________________________
void AliEveMUONData::LoadRecPoints(TTree* tree)
{
  //
  // load reconstructed points from the TreeR
  // load local trigger information
  //

  AliMUONVClusterStore *clusterStore = AliMUONVClusterStore::Create(*tree);
  clusterStore->Clear();
  clusterStore->Connect(*tree,kFALSE);

  tree->GetEvent(0);

  AliMUONVCluster *cluster;
  Int_t detElemId;
  Double_t clsX, clsY, clsZ, charge;

  for (Int_t ch = 0; ch < 10; ++ch) {

    if (fChambers[ch] == 0) continue;

    TIter next(clusterStore->CreateChamberIterator(ch,ch));

    while ( ( cluster = static_cast<AliMUONVCluster*>(next()) ) ) {

      detElemId = cluster->GetDetElemId();

      clsX   = cluster->GetX();
      clsY   = cluster->GetY();
      clsZ   = cluster->GetZ();
      charge = cluster->GetCharge();

      fChambers[ch]->RegisterCluster(detElemId,0,clsX,clsY,clsZ,charge);
      fChambers[ch]->RegisterCluster(detElemId,1,clsX,clsY,clsZ,charge);

    }

  }

  delete clusterStore;

}

//______________________________________________________________________________
void AliEveMUONData::LoadRecPointsFromESD(const Char_t *fileName)
{
  //
  // load reconstructed points stored in AliESDs.root
  // load local trigger information
  //

  TFile* esdFile = TFile::Open(fileName);
  if (!esdFile || !esdFile->IsOpen()) {
    cout << "opening ESD file " << fileName << "failed" << endl;
    return;
  }
  TTree* esdTree = (TTree*) esdFile->Get("esdTree");
  if (!esdTree) {
    cout << "no ESD tree found" << endl;
    esdFile->Close();
    return;
  }
  AliESDEvent* esdEvent = new AliESDEvent();
  esdEvent->ReadFromTree(esdTree);

  AliMUONVCluster *cluster;
  AliMUONTrackParam *trackParam;
  AliESDMuonTrack *esdTrack;
  AliMUONTrack muonTrack;
  Int_t detElemId, chamber, nTrackParam;
  Double_t clsX, clsY, clsZ, charge;
  
  if (esdTree->GetEvent(AliEveEventManager::Instance()->GetEventId()) <= 0) {
    cout << "fails to read ESD object for event " << AliEveEventManager::Instance()->GetEventId() << endl;
    return;
  }
    
  Int_t nTracks = Int_t(esdEvent->GetNumberOfMuonTracks());
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    esdTrack = esdEvent->GetMuonTrack(iTrack);
    if (!esdTrack->ContainTrackerData()) continue;
    AliMUONESDInterface::ESDToMUON(*esdTrack,muonTrack);
    nTrackParam = muonTrack.GetTrackParamAtCluster()->GetEntries();
    for(Int_t iCluster = 0; iCluster < nTrackParam; iCluster++) {
      trackParam = (AliMUONTrackParam *) muonTrack.GetTrackParamAtCluster()->At(iCluster);
      cluster = trackParam->GetClusterPtr();
      chamber   = cluster->GetChamberId();
      detElemId = cluster->GetDetElemId();
      charge  = cluster->GetCharge();
      clsX = cluster->GetX();
      clsY = cluster->GetY();
      clsZ = cluster->GetZ();
      
      fChambers[chamber]->RegisterCluster(detElemId,0,clsX,clsY,clsZ,charge);
      fChambers[chamber]->RegisterCluster(detElemId,1,clsX,clsY,clsZ,charge);
      
    }
  }

  delete esdEvent;
  
  esdFile->Close();

}

//______________________________________________________________________________
void AliEveMUONData::LoadHits(TTree* tree)
{
  //
  // load simulation hits from the TreeH
  //

  TClonesArray *hits = 0;
  AliMUONHit  *mhit;
  Int_t cha, detElemId, nhits, ntracks;
  Float_t hitX, hitY, hitZ;

  ntracks = tree->GetEntries();
  tree->SetBranchAddress("MUONHits",&hits);

  for (Int_t it = 0; it < ntracks; it++) {

    tree->GetEvent(it);
    nhits = hits->GetEntriesFast();

    for (Int_t ih = 0; ih < nhits; ih++) {

      mhit = (AliMUONHit*)hits->UncheckedAt(ih);
      hitX = mhit->X();
      hitY = mhit->Y();
      hitZ = mhit->Z();
      detElemId = mhit->DetElemId();
      cha = mhit->Chamber();

      RegisterTrack(mhit->GetTrack());

      fChambers[cha-1]->RegisterHit(detElemId,hitX,hitY,hitZ);

    }
  }

}

//______________________________________________________________________________
void AliEveMUONData::LoadDigits(TTree* tree)
{
  //
  // load digits from the TreeD
  //

  AliMUONVDigitStore *digitStore = AliMUONVDigitStore::Create(*tree);
  digitStore->Clear();
  digitStore->Connect(*tree,0);

  tree->GetEvent(0);

  AliMUONVDigit* digit;
  TIter next(digitStore->CreateIterator());

  Int_t cathode, detElemId, ix, iy, charge, chamber, adc;

  while ( ( digit = static_cast<AliMUONVDigit*>(next() ) ) )
    {
      cathode   = digit->Cathode();
      ix        = digit->PadX();
      iy        = digit->PadY();
      detElemId = digit->DetElemId();
      charge    = (Int_t)digit->Charge();
      adc       = digit->ADC();
      chamber   = detElemId/100 - 1;
      if (chamber > 9) {
	fChambers[chamber]->RegisterDigit(detElemId,cathode,ix,iy,charge);
      } else {
	fChambers[chamber]->RegisterDigit(detElemId,cathode,ix,iy,adc);
      }
    }

  delete digitStore;

}

//______________________________________________________________________________
void AliEveMUONData::LoadRaw(TString fileName)
{
  //
  // load raw data from fileName; tracker and trigger data
  //

  if (fgRawReader == 0) {
    // check extention to choose the rawdata file format
    if (fileName.EndsWith("/")) {
      fgRawReader = new AliRawReaderFile(fileName); // DDL files
    } else if (fileName.EndsWith(".root")) {
      fgRawReader = new AliRawReaderRoot(fileName); // ROOT file
    } else if (!fileName.IsNull()) {
      fgRawReader = new AliRawReaderDate(fileName); // DATE file
    }
  }

  fgRawReader->RewindEvents();
  fgRawReader->Reset();

  Int_t iEvent = 0;
  while (fgRawReader->NextEvent())
  {
    if (iEvent != AliEveEventManager::Instance()->GetEventId())
    {
      iEvent++;
      continue;
    }
    break;
  }

  AliMUONDigitMaker digitMaker;

  digitMaker.SetMakeTriggerDigits(kTRUE);

  AliMUONDigitStoreV1 digitStore;

  digitMaker.Raw2Digits(fgRawReader,&digitStore);

  AliMUONVDigit* digit;
  TIter next(digitStore.CreateIterator());

  Int_t cathode, detElemId, ix, iy, charge, chamber, adc;

  while ( ( digit = static_cast<AliMUONVDigit*>(next() ) ) )
  {
      cathode   = digit->Cathode();
      ix        = digit->PadX();
      iy        = digit->PadY();
      detElemId = digit->DetElemId();
      charge    = (Int_t)digit->Charge();
      adc       = digit->ADC();
      chamber   = detElemId/100 - 1;
      if (chamber > 9) {
	fChambers[chamber]->RegisterDigit(detElemId,cathode,ix,iy,charge);
      } else {
	fChambers[chamber]->RegisterDigit(detElemId,cathode,ix,iy,adc);
      }
  }
}

//______________________________________________________________________________
Int_t AliEveMUONData::GetTrack(Int_t index) const
{
  //
  // return track stack number for "index"-th track with hits in the chambers
  //

  if (index < 256) {
    return fTrackList[index];
  } else {
    return -1;
  }
}

//______________________________________________________________________________
AliEveMUONChamberData* AliEveMUONData::GetChamberData(Int_t chamber)
{
  //
  // return chamber data
  //

  if (chamber < 0 || chamber > 13) return 0;

  //if (fChambers[chamber] == 0) CreateChamber(chamber);

  return fChambers[chamber];
}
