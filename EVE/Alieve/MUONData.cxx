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
#include <AliMUONRawCluster.h>
#include <AliMUONDigit.h>
#include <AliMUONTriggerCrateStore.h>
#include <AliMUONData.h>

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
  Reve::ReferenceCount()
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
void MUONData::LoadDigits(TTree* tree)
{
  // 
  // load digits from the TreeD
  //

  Char_t branchname[30];
  TClonesArray *digits = 0;
  Int_t ndigits;
  AliMUONDigit  *mdig;
  Int_t cathode, detElemId, ix, iy, charge;

  for (Int_t c = 0; c < 14; ++c) {

    if (fChambers[c] == 0) continue;
    sprintf(branchname,"MUONDigits%d",c+1);
    tree->SetBranchAddress(branchname,&digits);
    tree->GetEntry(0);

    ndigits = digits->GetEntriesFast(); 

    for (Int_t id = 0; id < ndigits; id++) {
      mdig  = (AliMUONDigit*)digits->UncheckedAt(id);

      cathode   = mdig->Cathode();
      ix        = mdig->PadX();
      iy        = mdig->PadY();
      detElemId = mdig->DetElemId();      
      charge    = (Int_t)mdig->Signal();

      fChambers[c]->RegisterDigit(detElemId,cathode,ix,iy,charge);
      
    } // end digits loop

  }

}

//______________________________________________________________________
void MUONData::LoadRecPoints(TTree* tree)
{
  //
  // load reconstructed points from the TreeR
  // load local trigger information
  //

  Char_t branchname[30];
  TClonesArray *clusters = 0;
  Int_t nclusters;
  AliMUONRawCluster  *mcls;
  Int_t detElemId;
  Float_t clsX, clsY, clsZ, charge;

  for (Int_t c = 0; c < 10; ++c) {

    if (fChambers[c] == 0) continue;
    sprintf(branchname,"MUONRawClusters%d",c+1);
    tree->SetBranchAddress(branchname,&clusters);
    tree->GetEntry(0);

    nclusters = clusters->GetEntriesFast(); 

    for (Int_t ic = 0; ic < nclusters; ic++) {
      mcls  = (AliMUONRawCluster*)clusters->UncheckedAt(ic);

      detElemId = mcls->GetDetElemId();
      for (Int_t icath = 0; icath < 2; icath++) {
	clsX   = mcls->GetX(icath);
	clsY   = mcls->GetY(icath);
	clsZ   = mcls->GetZ(icath);
	charge = mcls->GetCharge(icath);

	fChambers[c]->RegisterCluster(detElemId,icath,clsX,clsY,clsZ,charge);
      }

    }

  }

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
  while (fgRawReader->NextEvent()) {
    if (iEvent != Alieve::gEvent->GetEventId()) {
      iEvent++;
      continue;
    }
    break;
  }

  AliMUONDigitMaker *digitMaker = new AliMUONDigitMaker();

  AliMUONTriggerCrateStore *crateManager = new AliMUONTriggerCrateStore();
  crateManager->ReadFromFile();

  AliMUONData *muonData = new AliMUONData(0x0,"MUON","MUON");

  digitMaker->SetDisplayFlag();
  digitMaker->SetCrateManager(crateManager);
  digitMaker->SetMUONData(muonData);
  muonData->SetDataContainer("D, GLT");

  digitMaker->Raw2Digits(fgRawReader);

  AliMUONDigit *digit;
  Int_t cathode, detElemId, ix, iy, charge, chamber, ndigits;
  for (chamber = 0; chamber < 14; chamber++) {
    ndigits = (Int_t)muonData->Digits(chamber)->GetEntriesFast();
    for (Int_t id = 0; id < ndigits; id++) {
      digit = static_cast<AliMUONDigit*>(muonData->Digits(chamber)->At(id));
      cathode   = digit->Cathode();
      ix        = digit->PadX();
      iy        = digit->PadY();
      detElemId = digit->DetElemId();      
      charge    = (Int_t)digit->Signal();
      chamber   = detElemId/100 - 1;
      fChambers[chamber]->RegisterDigit(detElemId,cathode,ix,iy,charge);
    }
  }

  delete muonData;
  delete crateManager;
  delete digitMaker;

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
