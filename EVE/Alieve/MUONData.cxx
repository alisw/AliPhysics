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

#include <AliMUONDigit.h>
#include <AliMUONRawStreamTracker.h>
#include <AliMUONRawStreamTrigger.h>
#include <AliMUONDDLTracker.h>
#include <AliMUONBlockHeader.h>
#include <AliMUONDspHeader.h>
#include <AliMUONBusStruct.h>
#include <AliMUONDDLTrigger.h>
#include <AliMUONDarcHeader.h>
#include <AliMUONRegHeader.h>
#include <AliMUONLocalStruct.h>
#include <AliMUONTriggerCrateStore.h>
#include <AliMUONTriggerCrate.h>
#include <AliMUONLocalTriggerBoard.h>
#include <AliMUONTriggerCircuit.h>
#include <mapping/AliMpDDLStore.h>
#include <mapping/AliMpVSegmentation.h>
#include <mapping/AliMpSegmentation.h>
#include <mapping/AliMpPad.h>
#include <mapping/AliMpDEManager.h>

#include "TTree.h"
#include "TString.h"
#include "TClonesArray.h"

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// MUONData
//

ClassImp(MUONData)

AliRawReader*            MUONData::fgRawReader        = 0;
AliMUONRawStreamTracker* MUONData::fgRawStreamTracker = 0;
AliMUONRawStreamTrigger* MUONData::fgRawStreamTrigger = 0;
AliMpDDLStore*           MUONData::fgBusPatchManager  = 0;

//______________________________________________________________________
MUONData::MUONData() :
  fChambers(14)
{
  //
  // Constructor
  //
  
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

      if (c > 9) {
	//printf("cha %d deid %d cath %1d ix %d iy %d q %d \n",c,detElemId,cathode,ix,iy,charge);  
      }

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
    fgRawStreamTracker = new AliMUONRawStreamTracker(fgRawReader);
    fgRawStreamTrigger = new AliMUONRawStreamTrigger(fgRawReader);
    fgBusPatchManager = AliMpDDLStore::Instance();
  }
  
  LoadRawTracker();
  LoadRawTrigger();

}

//______________________________________________________________________
void MUONData::LoadRawTracker()
{
  //
  // load raw data for the tracking chambers
  //

  fgRawReader->RewindEvents();

  AliMUONDigit* digit = new AliMUONDigit();

  Int_t maxEvent = 1000;
  Int_t minDDL = 0, maxDDL = 19;
  Int_t cathode, detElemId, ix, iy, iChamber;

  AliMUONDDLTracker*       ddlTracker = 0x0;
  AliMUONBlockHeader*      blkHeader  = 0x0;
  AliMUONDspHeader*        dspHeader  = 0x0;
  AliMUONBusStruct*        busStruct  = 0x0;

  Int_t iEvent = 0;
  Int_t dataSize, buspatchId;
  
  Event* aevent = Alieve::gEvent;

  while (fgRawReader->NextEvent()) {
    
    if (iEvent != aevent->GetEventId()) {
      iEvent++;
      continue;
    }
    
    if (iEvent == maxEvent)
      break;
    
    // read DDL while < 20 DDL
    while(fgRawStreamTracker->NextDDL()) {
      
      if (fgRawStreamTracker->GetDDL() < minDDL || 
	  fgRawStreamTracker->GetDDL() > maxDDL)
	continue;
      
      //printf("\niDDL %d\n", fgRawStreamTracker->GetDDL());
      
      ddlTracker =  fgRawStreamTracker->GetDDLTracker();
      
      // loop over block structure
      Int_t nBlock = ddlTracker->GetBlkHeaderEntries();
      for(Int_t iBlock = 0; iBlock < nBlock ;iBlock++){
	
	blkHeader = ddlTracker->GetBlkHeaderEntry(iBlock);
	//printf("Block Total length %d\n",blkHeader->GetTotalLength());
	
	// loop over DSP structure
	Int_t nDsp = blkHeader->GetDspHeaderEntries();
	for(Int_t iDsp = 0; iDsp < nDsp ;iDsp++){   //DSP loop
	  
	  dspHeader =  blkHeader->GetDspHeaderEntry(iDsp);
	  //   printf("Dsp length %d even word %d\n",dspHeader->GetTotalLength(), dspHeader->GetEventWord());
	  
	  // loop over BusPatch structure
	  Int_t nBusPatch = dspHeader->GetBusPatchEntries();
	  for(Int_t iBusPatch = 0; iBusPatch < nBusPatch; iBusPatch++) {  
	    
	    busStruct = dspHeader->GetBusPatchEntry(iBusPatch);
	    
	    //printf("busPatchId %d", busStruct->GetBusPatchId());
	    //printf(" BlockId %d", busStruct->GetBlockId());
	    //printf(" DspId %d\n", busStruct->GetDspId());
	    
	    // loop over data
	    dataSize = busStruct->GetLength();
	    buspatchId = busStruct->GetBusPatchId();
	    for (Int_t iData = 0; iData < dataSize; iData++) {
	      
	      Int_t  manuId    = busStruct->GetManuId(iData);
	      Int_t  channelId = busStruct->GetChannelId(iData);
	      Int_t  charge    = busStruct->GetCharge(iData);
	      //printf("manuId: %d, channelId: %d charge: %d\n", manuId, channelId, charge);
	      // set digit charge
	      digit->SetSignal(charge);
	      digit->SetPhysicsSignal(charge);
	      digit->SetADC(charge);
	      // Get Back the hits at pads
	      Int_t error;
	      error = GetTrackerMapping(buspatchId,manuId,channelId,digit); 
	      if (error) {
		printf("Mapping Error\n");
		continue;
	      }

	      cathode   = digit->Cathode();
	      ix        = digit->PadX();
	      iy        = digit->PadY();
	      detElemId = digit->DetElemId();      
	      charge    = (Int_t)digit->Signal();
	      iChamber  = detElemId/100 - 1;

	      fChambers[iChamber]->RegisterDigit(detElemId,cathode,ix,iy,charge);
	      
	    } // iData
	  } // iBusPatch
	} // iDsp
      } // iBlock
    } // NextDDL

    break;

  }  // end event loop
  
  delete digit;

}

//______________________________________________________________________
void MUONData::LoadRawTrigger()
{
  // 
  // load raw data for the trigger chambers
  //

  fgRawReader->RewindEvents();

  Int_t maxEvent = 1000;
  Int_t minDDL = 0, maxDDL = 1;
  Int_t detElemId, iChamber, cathode, charge, ix, iy;

  AliMUONDDLTrigger*       ddlTrigger  = 0x0;
  AliMUONDarcHeader*       darcHeader  = 0x0;
  AliMUONRegHeader*        regHeader   = 0x0;
  AliMUONLocalStruct*      localStruct = 0x0;
  
  // crate manager
  AliMUONTriggerCrateStore* crateManager = new AliMUONTriggerCrateStore();   
  crateManager->ReadFromFile();

  // Loop over events  
  Int_t iEvent = 0;
  TList digitList;
  
  Event* aevent = Alieve::gEvent;
  
  while (fgRawReader->NextEvent()) {
    
    if (iEvent != aevent->GetEventId()) {
      iEvent++;
      continue;
    }
    
    if (iEvent == maxEvent)
      break;
    
    // read DDL while < 2 DDL
    while(fgRawStreamTrigger->NextDDL()) {
      
      if (fgRawStreamTrigger->GetDDL() < minDDL || 
	  fgRawStreamTrigger->GetDDL() > maxDDL)
	continue;
      
      //printf("\niDDL %d\n", fgRawStreamTrigger->GetDDL());
      
      ddlTrigger = fgRawStreamTrigger->GetDDLTrigger();
      darcHeader = ddlTrigger->GetDarcHeader();
      
      //printf("Global output %x\n", (Int_t)darcHeader->GetGlobalOutput());
      
      // loop over regional structures
      Int_t nReg = darcHeader->GetRegHeaderEntries();
      for(Int_t iReg = 0; iReg < nReg ;iReg++){   //REG loop
	
	//printf("RegionalId %d\n", iReg);
	
	regHeader =  darcHeader->GetRegHeaderEntry(iReg);
	//  printf("Reg length %d\n",regHeader->GetHeaderLength());
	
	// crate info
	AliMUONTriggerCrate* crate = crateManager->Crate(fgRawStreamTrigger->GetDDL(), iReg);
	TObjArray *boards = crate->Boards();
	
	// loop over local structures
	Int_t nLocal = regHeader->GetLocalEntries();
	for(Int_t iLocal = 0; iLocal < nLocal; iLocal++) {  
	  
	  localStruct = regHeader->GetLocalEntry(iLocal);
	  
	  // check if trigger 
	  if (localStruct->GetTriggerY() == 0) { // no empty data
	    
	    // local trigger circuit number
	    AliMUONLocalTriggerBoard* localBoard = (AliMUONLocalTriggerBoard*)boards->At(iLocal+1);
	    
	    //printf("LocalId %d\n", localStruct->GetId());
	    /*
	    Int_t iLocCard  = localBoard->GetNumber();
	    Int_t loStripX  = (Int_t)localStruct->GetXPos();
	    Int_t loStripY  = (Int_t)localStruct->GetYPos();
	    Int_t loDev     = (Int_t)localStruct->GetXDev();
	    */
	    //printf("iLocCard: %d, XPos: %d, YPos: %d Dev: %d\n", iLocCard, loStripX, loStripY, loDev);

	    digitList.Clear();
	    if ( GetTriggerMapping(localBoard, localStruct, digitList) ) {
	      for (Int_t iEntry = 0; iEntry < digitList.GetEntries(); iEntry++) {

		AliMUONDigit* digit = (AliMUONDigit*)digitList.At(iEntry);
		cathode   = digit->Cathode();
		ix        = digit->PadX();
		iy        = digit->PadY();
		detElemId = digit->DetElemId();      
		charge    = (Int_t)digit->Signal();
		iChamber  = detElemId/100 - 1;
		
		//printf("cha %d deid %d cath %1d ix %d iy %d q %d \n",iChamber,detElemId,cathode,ix,iy,charge);  

		fChambers[iChamber]->RegisterDigit(detElemId,cathode,ix,iy,charge);

	      }

	    }

	  }
	} // iLocal
      } // iReg
    } // NextDDL

    break;

  }  // end event loop

  delete crateManager;

}

//______________________________________________________________________
Int_t MUONData::GetTrackerMapping(Int_t buspatchId, UShort_t manuId, UChar_t channelId, AliMUONDigit* digit)
{
  //
  // decode digits mapping for the tracking chambers
  //
  
  // getting DE from buspatch
  Int_t detElemId = fgBusPatchManager->GetDEfromBus(buspatchId);
  //AliDebug(3,Form("detElemId: %d busPatchId %d\n", detElemId, buspatchId));

  const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId, manuId);  
  AliMpPad pad = seg->PadByLocation(AliMpIntPair(manuId,channelId),kTRUE);

  if (!pad.IsValid())
  {
    printf("No pad for detElemId: %d, busPatchId %d, manuId: %d, channelId: %d\n",detElemId, buspatchId, manuId, channelId);
    
    return 1;
  } // return error

  // Getting padX, padY and cathode number.
  Int_t padX = pad.GetIndices().GetFirst();
  Int_t padY = pad.GetIndices().GetSecond();
  Int_t iCath = AliMpDEManager::GetCathod(detElemId,seg->PlaneType());

  // storing into digits
  digit->SetPadX(padX);
  digit->SetPadY(padY);
  digit->SetCathode(iCath);
  digit->SetDetElemId(detElemId);
  digit->SetElectronics(manuId,channelId);
  
  //printf("detElemId: %d, busPatchId %d, manuId: %d, channelId: %d, padx: %d pady %d\n",detElemId, buspatchId, manuId, channelId, padX, padY);
  
  return 0;

}

//______________________________________________________________________
Int_t MUONData::GetTriggerMapping(AliMUONLocalTriggerBoard* localBoard, 
				  AliMUONLocalStruct* localStruct,
				  TList& digitList)
{
  //
  // decode digits mapping for the trigger chambers
  //

  Int_t detElemId;
  Int_t nBoard;
  Int_t iCath = -1;
  Int_t iChamber = 0;
  Int_t xyPattern = 0;

  // loop over x1-4 and y1-4
  for (Int_t icase = 0; icase < 8; icase++) {

    // get chamber, cathode and associated trigger response pattern
    GetTriggerChamber(localStruct, xyPattern, iChamber, iCath, icase);
  
    if (!xyPattern) continue;

    // get detElemId
    AliMUONTriggerCircuit triggerCircuit;
    detElemId = triggerCircuit.DetElemId(iChamber, localBoard->GetName());
    nBoard    = localBoard->GetNumber();

    const AliMpVSegmentation* seg 
      = AliMpSegmentation::Instance()
        ->GetMpSegmentation(detElemId, AliMp::GetCathodType(iCath));  

    // loop over the 16 bits of pattern
    for (Int_t ibitxy = 0; ibitxy < 16; ibitxy++) {
    
      if ((xyPattern >> ibitxy) & 0x1) {

	// not quite sure about this
	Int_t offset = 0;
	if (iCath && localBoard->GetSwitch(6)) offset = -8;

	AliMpPad pad = seg->PadByLocation(AliMpIntPair(nBoard,ibitxy+offset),kTRUE);

	AliMUONDigit* digit = new  AliMUONDigit();
	if (!pad.IsValid()) {
	  AliWarning(Form("No pad for detElemId: %d, nboard %d, ibitxy: %d\n",
			  detElemId, nBoard, ibitxy));
	  continue;
	} // 

	Int_t padX = pad.GetIndices().GetFirst();
	Int_t padY = pad.GetIndices().GetSecond();

	// file digit
	digit->SetSignal(1);
	digit->SetPadX(padX);
	digit->SetPadY(padY);
	digit->SetCathode(iCath);
	digit->SetDetElemId(detElemId);
	digit->SetElectronics(nBoard, ibitxy);
	digitList.Add(digit);
	
      }// xyPattern
    }// ibitxy
  }// case

  return 1;

}

//____________________________________________________________________
void MUONData::GetTriggerChamber(AliMUONLocalStruct* localStruct, Int_t& xyPattern, Int_t& iChamber, Int_t& iCath, Int_t icase)
{
  //
  // extract digits pattern
  //  

  // get chamber & cathode number, (chamber starts at 0 !)
  switch(icase) {
  case 0: 
    xyPattern =  localStruct->GetX1();
    iCath = 0;
    iChamber = 10;
    break;
  case 1: 
    xyPattern =  localStruct->GetX2();
    iCath = 0;
    iChamber = 11;
    break;
  case 2: 
    xyPattern =  localStruct->GetX3();
    iCath = 0;
    iChamber = 12;
    break;
  case 3: 
    xyPattern =  localStruct->GetX4();
    iCath = 0;
    iChamber = 13;
    break;
  case 4: 
    xyPattern =  localStruct->GetY1();
    iCath = 1;
    iChamber = 10;
    break;
  case 5: 
    xyPattern =  localStruct->GetY2();
    iCath = 1;
    iChamber = 11;
    break;
  case 6: 
    xyPattern =  localStruct->GetY3();
    iCath = 1;
    iChamber = 12;
    break;
  case 7: 
    xyPattern =  localStruct->GetY4();
    iCath = 1;
    iChamber = 13;
    break;
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
