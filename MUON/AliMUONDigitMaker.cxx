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

/// \class AliMUONDigitMaker
/// MUON Digit maker from rawdata.
///
/// Raw2Digits:
/// Using real mapping  for tracker
/// Indranil Das (Adapted for runloader: Ch. Finck) july 05
///
/// Implemented non-constant buspatch numbers for tracking
/// with correct DDL id.
/// (Ch. Finck, dec 05)
///
/// Add reader for scaler trigger events
/// Use memcpy instead of assignment elt by elt
/// (Ch. Finck, Jan 06)
///
/// Using new interface with AliMUONRawStreamTracker(Trigger)
/// (New interface of AliMUONRawReader class)
/// (further details could be found in Alice-note)
/// (Ch. Finck, March 06)
///
/// Add (S)Digit maker tracker (for free)
/// and for trigger. Create trigger inverse mapping.
/// (Ch. Finck, oct 06) 

#include "AliMUONDigitMaker.h"
#include "AliMUONDigit.h"

#include "AliMUONConstants.h"
#include "AliMUONData.h"

#include "AliMUONRawStreamTracker.h"
#include "AliMUONDDLTracker.h"
#include "AliMUONDspHeader.h"
#include "AliMUONBlockHeader.h"
#include "AliMUONBusStruct.h"

#include "AliMUONRawStreamTrigger.h"
#include "AliMUONDDLTrigger.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONRegHeader.h"
#include "AliMUONLocalStruct.h"

#include "AliMUONTriggerCrateStore.h"
#include "AliMUONTriggerCrate.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONTriggerCircuitNew.h"

#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpPad.h"
#include "AliMpDEManager.h"
#include "AliMpBusPatch.h"

#include "AliRawReader.h"
#include "AliRawDataHeader.h"
#include "AliLog.h"
#include "AliRun.h"

#include <TList.h>


/// \cond CLASSIMP
ClassImp(AliMUONDigitMaker) // Class implementation in ROOT context
/// \endcond

//__________________________________________________________________________
AliMUONDigitMaker::AliMUONDigitMaker(Bool_t flag)
  : TObject(),
    fMUONData(0x0),
    fBusPatchManager(new AliMpBusPatch()),
    fScalerEvent(kFALSE),
    fDigitFlag(flag),
    fRawStreamTracker(new AliMUONRawStreamTracker()),    
    fRawStreamTrigger(new AliMUONRawStreamTrigger()),    
    fDigit(new AliMUONDigit()),
    fLocalTrigger(new AliMUONLocalTrigger()),
    fGlobalTrigger(new AliMUONGlobalTrigger()),
    fCrateManager(0x0),
    fTrackerTimer(),
    fTriggerTimer(),
    fMappingTimer()
{
  /// ctor with AliMUONData as argument
  /// for reconstruction

  AliDebug(1,"");

  // Standard Constructor

  // bus patch 
  fBusPatchManager->ReadBusPatchFile();

 
  fTrackerTimer.Start(kTRUE); fTrackerTimer.Stop();
  fTriggerTimer.Start(kTRUE); fTriggerTimer.Stop();
  fMappingTimer.Start(kTRUE); fMappingTimer.Stop();

}

//__________________________________________________________________________
AliMUONDigitMaker::~AliMUONDigitMaker()
{
  /// clean up
  /// and time processing measure

  delete fRawStreamTracker;
  delete fRawStreamTrigger;

  delete fDigit;
  delete fLocalTrigger;
  delete fGlobalTrigger;

  delete fBusPatchManager;

  AliDebug(1, Form("Execution time for MUON tracker : R:%.2fs C:%.2fs",
               fTrackerTimer.RealTime(),fTrackerTimer.CpuTime()));
  AliDebug(1, Form("   Execution time for MUON tracker (mapping calls part) "
               ": R:%.2fs C:%.2fs",
               fMappingTimer.RealTime(),fMappingTimer.CpuTime()));
  AliDebug(1, Form("Execution time for MUON trigger : R:%.2fs C:%.2fs",
               fTriggerTimer.RealTime(),fTriggerTimer.CpuTime()));

  return;
}

//____________________________________________________________________
Int_t AliMUONDigitMaker::Raw2Digits(AliRawReader* rawReader)
{
  /// Main method to creates digit
  /// for tracker 
  /// and trigger

  // generate digits
  ReadTrackerDDL(rawReader);

  // generate trigger
  ReadTriggerDDL(rawReader);

  return kTRUE;

}

//____________________________________________________________________
Int_t AliMUONDigitMaker::ReadTrackerDDL(AliRawReader* rawReader)
{

  /// reading tracker DDL
  /// filling the TClonesArray in MUONData

  fTrackerTimer.Start(kFALSE);

  // elex info
  Int_t    buspatchId;
  UChar_t  channelId;
  UShort_t manuId;
  Char_t   parity;
  UShort_t charge; 
  Int_t    dataSize;

  Int_t iChamber;

  AliMUONDDLTracker*   ddlTracker = 0x0;
  AliMUONBlockHeader*  blkHeader  = 0x0;
  AliMUONDspHeader*    dspHeader  = 0x0;
  AliMUONBusStruct*    busStruct  = 0x0;


  fRawStreamTracker->SetReader(rawReader);

  while(fRawStreamTracker->NextDDL()) {

    ddlTracker =  fRawStreamTracker->GetDDLTracker();

    Int_t nBlock = ddlTracker->GetBlkHeaderEntries();
    for(Int_t iBlock = 0; iBlock < nBlock ;iBlock++){

      blkHeader = ddlTracker->GetBlkHeaderEntry(iBlock);
 
      Int_t nDsp = blkHeader->GetDspHeaderEntries();

      for(Int_t iDsp = 0; iDsp < nDsp ;iDsp++){   //DSP loop

	dspHeader =  blkHeader->GetDspHeaderEntry(iDsp);

	Int_t nBusPatch = dspHeader->GetBusPatchEntries();

	for(Int_t iBusPatch = 0; iBusPatch < nBusPatch; iBusPatch++) {  

	  busStruct = dspHeader->GetBusPatchEntry(iBusPatch);

	  dataSize   = busStruct->GetLength();
	  buspatchId = busStruct->GetBusPatchId();

	  for (Int_t iData = 0; iData < dataSize; iData++) {

	    // digits info
	    parity    = busStruct->GetParity(iData); // test later for parity
	    manuId    = busStruct->GetManuId(iData);
	    channelId = busStruct->GetChannelId(iData);
	    charge    = busStruct->GetCharge(iData);
	    // set charge
	    fDigit->SetSignal(charge);
	    fDigit->SetPhysicsSignal(charge);
	    fDigit->SetADC(charge);

	    // Get Back the hits at pads
	    Int_t error = GetMapping(buspatchId,manuId,channelId,fDigit); 
	    if (error) {
	      AliWarning("Mapping Error\n");
	      continue;
	    }
	    // debugging 
	    if (AliLog::GetGlobalDebugLevel() == 3) {
	      Int_t padX  = fDigit->PadX();
	      Int_t padY  = fDigit->PadY();
	      Int_t iCath = fDigit->Cathode();  
	      Int_t idDE  = fDigit->DetElemId();

	      AliDebug(1,Form("output  IdDE %d busPatchid %d PadX %d PadY %d iCath %d \n", 
			      idDE, buspatchId, padX, padY, iCath));
		
	      AliDebug(3,Form("idDE %d Padx %d Pady %d, Cath %d, charge %d",
			      idDE, padX, padY, iCath, charge));
	    }

	    // fill digits
	    iChamber = AliMpDEManager::GetChamberId(fDigit->DetElemId());

	    if (fDigitFlag)
	      fMUONData->AddDigit(iChamber, *fDigit);
	    else
	      fMUONData->AddSDigit(iChamber, *fDigit);


	  } // iData
	} // iBusPatch
      } // iDsp
    } // iBlock
  } // NextDDL

  fTrackerTimer.Stop();

  return kTRUE;
}
//____________________________________________________________________
Int_t AliMUONDigitMaker::GetMapping(Int_t busPatchId, UShort_t manuId, 
					 UChar_t channelId, AliMUONDigit* digit )
{
  /// mapping  for tracker

  fMappingTimer.Start(kFALSE);
  
  // getting DE from buspatch
  Int_t detElemId = fBusPatchManager->GetDEfromBus(busPatchId);
  AliDebug(3,Form("detElemId: %d busPatchId %d\n", detElemId, busPatchId));

  const AliMpVSegmentation* seg 
    = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId, manuId);  
  AliMpPad pad = seg->PadByLocation(AliMpIntPair(manuId,channelId),kTRUE);

  if (!pad.IsValid())
  {
    AliWarning(Form("No pad for detElemId: %d, busPatchId %d, manuId: %d, channelId: %d\n",
		  detElemId, busPatchId, manuId, channelId));
    fMappingTimer.Stop();
    return kTRUE;
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
  
  AliDebug(3,Form("detElemId: %d, busPatchId %d, manuId: %d, channelId: %d, padx: %d pady %d\n",
		  detElemId, busPatchId, manuId, channelId, padX, padY));
  StdoutToAliDebug(3,digit->Print(););
  
  fMappingTimer.Stop();
  return kFALSE;
}

//____________________________________________________________________
Int_t AliMUONDigitMaker::ReadTriggerDDL(AliRawReader* rawReader)
{
  /// reading tracker DDL
  /// filling the TClonesArray in MUONData

  AliMUONDDLTrigger*       ddlTrigger      = 0x0;
  AliMUONDarcHeader*       darcHeader      = 0x0;
  AliMUONRegHeader*        regHeader       = 0x0;
  AliMUONLocalStruct*      localStruct     = 0x0;

  Int_t loCircuit;
  TList digitList;


  fTriggerTimer.Start(kFALSE);

  fRawStreamTrigger->SetReader(rawReader);

  while(fRawStreamTrigger->NextDDL()) {

    ddlTrigger = fRawStreamTrigger->GetDDLTrigger();
    darcHeader = ddlTrigger->GetDarcHeader();

    // fill global trigger information in Digit Tree
    if (fDigitFlag) {
      if (darcHeader->GetGlobalFlag()) {
	fGlobalTrigger->SetFromGlobalResponse(darcHeader->GetGlobalOutput());
	fMUONData->AddGlobalTrigger(*fGlobalTrigger);
      }
    }

    Int_t nReg = darcHeader->GetRegHeaderEntries();

    for(Int_t iReg = 0; iReg < nReg ;iReg++){   //reg loop

      // crate info
      if (!fCrateManager) AliFatal("Crate Store not defined");
      AliMUONTriggerCrate* crate = fCrateManager->Crate(fRawStreamTrigger->GetDDL(), iReg);
  
      if (!crate) 
	AliWarning(Form("Missing crate number %d in DDL %d\n", iReg, fRawStreamTrigger->GetDDL()));

      TObjArray *boards  = crate->Boards();


      regHeader =  darcHeader->GetRegHeaderEntry(iReg);

      Int_t nLocal = regHeader->GetLocalEntries();
      for(Int_t iLocal = 0; iLocal < nLocal; iLocal++) {  

	localStruct = regHeader->GetLocalEntry(iLocal);

	// if card has triggered
	if (localStruct->GetTriggerY() == 0) {


	  AliMUONLocalTriggerBoard* localBoard = 
	    (AliMUONLocalTriggerBoard*)boards->At(localStruct->GetId()+1);

	  loCircuit = localBoard->GetNumber();

	  if (fDigitFlag) {
	    // fill local trigger
	    fLocalTrigger->SetLocalStruct(loCircuit, *localStruct);

	    fMUONData->AddLocalTrigger(*fLocalTrigger);

	  } else {
	    // Make SDigit

	    digitList.Clear();
	    
	    if( TriggerDigits(localBoard, localStruct, digitList) ) {

	      for (Int_t iEntry = 0; iEntry < digitList.GetEntries(); iEntry++) {

		AliMUONDigit* digit = (AliMUONDigit*)digitList.At(iEntry);
		
		// filling S container
		Int_t iChamber = AliMpDEManager::GetChamberId(digit->DetElemId());
		fMUONData->AddSDigit(iChamber, *digit);

	      }

	    } // trigger digits
	  } // S flag

	} // if triggerY
      } // iLocal
    } // iReg
  } // NextDDL

  fTriggerTimer.Stop();

  return kTRUE;

}
//____________________________________________________________________
void AliMUONDigitMaker::GetTriggerChamber(AliMUONLocalStruct* localStruct, Int_t& xyPattern, 
					  Int_t& iChamber, Int_t& iCath, Int_t icase)
{
  /// get chamber & cathode number, (chamber starts at 0 !)

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
//____________________________________________________________________
Int_t AliMUONDigitMaker::TriggerDigits(AliMUONLocalTriggerBoard* localBoard, 
				       AliMUONLocalStruct* localStruct,
				       TList& digitList)
{
  /// make (S)Digit for trigger

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
    AliMUONTriggerCircuitNew triggerCircuit;
    detElemId = triggerCircuit.DetElemId(iChamber, localBoard->GetName());
    nBoard    = localBoard->GetNumber();

    const AliMpVSegmentation* seg 
      = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, iCath);  

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
	digit->SetPadX(padX);
	digit->SetPadY(padY);
	digit->SetCathode(iCath);
	digit->SetDetElemId(detElemId);
	digit->SetElectronics(nBoard, ibitxy);
	digitList.Add(digit);
	
      }// xyPattern
    }// ibitxy
  }// case

  return kTRUE;
} 
//____________________________________________________________________
void  AliMUONDigitMaker::GetCrateName(Char_t* name, Int_t iDDL, Int_t iReg) const
{
  /// set crate name from DDL & reg number
  /// method same as in RawWriter, not so nice
  /// should be put in AliMUONTriggerCrateStore

      switch(iReg) {
      case 0:
      case 1:
	sprintf(name,"%d", iReg+1);
	break;
      case 2:
	strcpy(name, "2-3");
	break;
      case 3:
      case 4:
      case 5:
      case 6:
      case 7:
	sprintf(name,"%d", iReg);
	break;
      }

      // crate Right for first DDL
      if (iDDL == 0)
	strcat(name, "R");
      else 
	strcat(name, "L"); 
}
