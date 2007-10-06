//
// Sources:
//
// GetTrackerMapping = AliMUONDigitMaker::GetMapping
// GetTriggerMapping = AliMUONDigitMaker::TriggerDigits
// GetTriggerChamber = AliMUONDigitMaker::GetTriggerChamber
// LoadRawTracker    = MUONRawStreamTracker.C
// LoadRawTrigger    = MUONRawStreamTrigger.C
//

#include "MUONData.h"

#include <Alieve/MUONChamberData.h>
#include <Alieve/EventAlieve.h>

#include <AliRawReader.h>
#include <AliRawReaderFile.h>
#include <AliRawReaderDate.h>
#include <AliRawReaderRoot.h>

#include <AliLog.h>

#include <AliMUONDigitMaker.h>
#include <AliMUONHit.h>
#include <AliMUONVCluster.h>
#include "AliMUONVClusterStore.h"
#include <AliMUONVDigit.h>
#include "AliMUONDigitStoreV1.h"
#include "AliMUONVDigitStore.h"
#include "TTree.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TList.h"

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// MUONData
//

ClassImp(MUONData)

AliRawReader*            MUONData::fgRawReader        = 0;

//______________________________________________________________________
MUONData::MUONData() :
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

//______________________________________________________________________
MUONData::~MUONData()
{
  //
  // Destructor
  //

  DeleteAllChambers();

}

//______________________________________________________________________
void MUONData::Reset()
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

//______________________________________________________________________
MUONData::MUONData(const MUONData &mdata) :
  TObject(mdata),
  Reve::ReferenceCount(),
  fChambers(14),
  fNTrackList(0)
{
  //
  // Copy constructor
  //

}

//______________________________________________________________________
MUONData& MUONData::operator=(const MUONData &mdata)
{
  //
  // Assignment operator
  //

  if (this != &mdata) {

  }

  return *this;

}

//______________________________________________________________________
void MUONData::CreateChamber(Int_t chamber)
{
  // 
  // create data for the chamber with id=chamber (0 to 13)
  //

  if (fChambers[chamber] == 0)
    fChambers[chamber] = new MUONChamberData(chamber);

}

//______________________________________________________________________
void MUONData::CreateAllChambers()
{
  //
  // create all 14 chambers data
  //

  for (Int_t c = 0; c < 14; ++c)
    CreateChamber(c);

}

//______________________________________________________________________
void MUONData::DropAllChambers()
{
  // 
  // release data from all chambers 
  //

  for (Int_t c = 0; c < 14; ++c) {

    if (fChambers[c] != 0)
      fChambers[c]->DropData();

  }

}

//______________________________________________________________________
void MUONData::DeleteAllChambers()
{
  //
  // delete all chambers data
  //

  for (Int_t c = 0; c < 14; ++c) {

    delete fChambers[c];
    fChambers[c] = 0;

  }

}

//______________________________________________________________________
void MUONData::RegisterTrack(Int_t track)
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

//______________________________________________________________________
void MUONData::LoadRecPoints(TTree* tree)
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

//______________________________________________________________________
void MUONData::LoadHits(TTree* tree)
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

//______________________________________________________________________
void MUONData::LoadDigits(TTree* tree)
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

//______________________________________________________________________
void MUONData::LoadRaw(TString fileName)
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
    if (iEvent != Alieve::gEvent->GetEventId()) 
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

//______________________________________________________________________
Int_t MUONData::GetTrack(Int_t index)
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

//______________________________________________________________________
MUONChamberData* MUONData::GetChamberData(Int_t chamber)
{
  //
  // return chamber data
  //

  if (chamber < 0 || chamber > 13) return 0;

  //if (fChambers[chamber] == 0) CreateChamber(chamber);

  return fChambers[chamber];

}
