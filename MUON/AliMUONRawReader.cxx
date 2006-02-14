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

////////////////////////////////////
//
// MUON Raw Data reader in ALICE-MUON
// Class version 3 (further details could be found in Alice-note)
//
// Implemented non-constant buspatch numbers for tracking
// with correct DDL id (first guess)
// (Ch. Finck, dec 2005)
//
//
// Raw2Digits:
// Using real mapping  for tracker
// Indranil Das (Adapted for runloader: Ch. Finck) july 05
// Add reader for scaler trigger events
// Use memcpy instead of assignment elt by elt
// (Ch. Finck, Jan 06)
// 
////////////////////////////////////

#include <fstream>
#include <string>

#include <TClonesArray.h>

#include "AliLoader.h"
#include "AliBitPacking.h" 
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliRun.h"

#include "AliMUON.h"
#include "AliMUONRawReader.h"
#include "AliMUONDigit.h"

#include "AliMUONConstants.h"
#include "AliMUONData.h"

#include "AliMUONSubEventTracker.h"
#include "AliMUONScalerEventTrigger.h"
#include "AliMUONSubEventTrigger.h"
#include "AliMUONDDLTracker.h"
#include "AliMUONDDLTrigger.h"

#include "AliMUONLocalTrigger.h"
#include "AliMUONGlobalTrigger.h"

#include "AliMUONGeometrySegmentation.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryStore.h"
#include "AliMpSegFactory.h"
#include "AliMpPlaneType.h"
#include "AliMpVSegmentation.h"
#include "AliMpHelper.h"
#include "AliMpPad.h"


ClassImp(AliMUONRawReader) // Class implementation in ROOT context
//__________________________________________________________________________
AliMUONRawReader::AliMUONRawReader(AliLoader* loader,  AliMUONData* data)
  : TObject(),
    fScalerEvent(kFALSE)
{
  // Standard Constructor
 
  // initialize loader's
  fLoader = loader;

  // initialize segmentation factory
  fSegFactory = new AliMpSegFactory();

  // initialize container
  fMUONData  = data;

  // ddl pointer
  fDDLTracker = new AliMUONDDLTracker();
  fDDLTrigger = new AliMUONDDLTrigger();

  fBusPatchManager = new AliMpBusPatch();
  fBusPatchManager->ReadBusPatchFile();

}

//__________________________________________________________________________
AliMUONRawReader::AliMUONRawReader()
  : TObject(),
    fMUONData(0),
    fLoader(0),    
    fSegFactory(0),
    fDDLTracker(0),
    fDDLTrigger(0),
    fBusPatchManager(0),
    fScalerEvent(kFALSE)
{
  // Default Constructor
  
}

//_______________________________________________________________________
AliMUONRawReader::AliMUONRawReader (const AliMUONRawReader& rhs)
  : TObject(rhs)
{
// Protected copy constructor

  AliFatal("Not implemented.");
}

//_______________________________________________________________________
AliMUONRawReader & 
AliMUONRawReader::operator=(const AliMUONRawReader& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}

//__________________________________________________________________________
AliMUONRawReader::~AliMUONRawReader(void)
{
  if (fSegFactory) 
    fSegFactory->DeleteSegmentations();
  delete fSegFactory;  

  if (fDDLTracker)
    delete fDDLTracker;
  if (fDDLTrigger)
    delete fDDLTrigger;

  fBusPatchManager->Delete();

  return;
}

//____________________________________________________________________
Int_t AliMUONRawReader::Raw2Digits(AliRawReader* rawReader)
{

  // generate digits
  ReadTrackerDDL(rawReader);

  // generate trigger
  ReadTriggerDDL(rawReader);

  return kTRUE;

}

//____________________________________________________________________
Int_t AliMUONRawReader::ReadTrackerDDL(AliRawReader* rawReader)
{
  // reading tracker DDL
  // filling the TClonesArray in MUONData
  //

  AliMUONSubEventTracker* subEventTracker = new AliMUONSubEventTracker();
  AliMUONDigit* digit = new AliMUONDigit();


  //Read Header Size of DDL,Block,DSP and BusPatch (put k before constant imb'cile)

  Int_t kDDLHeaderSize      = fDDLTracker->GetHeaderSize();
  Int_t kBlockHeaderSize    = fDDLTracker->GetBlkHeaderLength();
  Int_t kDspHeaderSize      = fDDLTracker->GetDspHeaderLength();
  Int_t kBusPatchHeaderSize = subEventTracker->GetHeaderLength();

  Int_t totalDDLSize, totalBlockSize, totalDspSize , totalBusPatchSize, dataSize; 


  Int_t iBusPerDSP[5];//number of bus patches per DSP
  Int_t iDspMax; //number max of DSP per block

  // minimum data size (only header's)
  Int_t blankDDLSize;
  Int_t blankBlockSize;
  Int_t blankDspSize;  

  for(Int_t iDDL = 0; iDDL < 20; iDDL++) { // DDL loop

    AliDebug(3, Form("Chamber %d\n", iDDL/2 +1 ));

    // getting DSP info
    fBusPatchManager->GetDspInfo(iDDL/2, iDspMax, iBusPerDSP);

    //   Each DDL is made with 2 Blocks each of which consists of 5 DSP's at most and each of DSP has at most 5 buspatches.
    //   This information is used to calculate the size of headers (DDL,Block and DSP) which has no interesting data.
    blankDDLSize   = kDDLHeaderSize + 2*kBlockHeaderSize + 2*iDspMax*kDspHeaderSize;
    blankBlockSize = kBlockHeaderSize + iDspMax*kDspHeaderSize;

    for (Int_t i = 0; i < iDspMax; i++) {
      blankDDLSize   += 2*iBusPerDSP[i]*kBusPatchHeaderSize;
      blankBlockSize +=   iBusPerDSP[i]*kBusPatchHeaderSize;
    }

    rawReader->Select(0X9, iDDL, iDDL);  //Select the DDL file to be read  

    rawReader->ReadHeader();

    totalDDLSize = (rawReader->GetDataSize() + sizeof(AliRawDataHeader))/4; // 4 is multiplied to convert byte 2 words

    if(totalDDLSize > blankDDLSize) {      // Compare the DDL header with an empty DDL header size to read the file

      Int_t totalDataWord = rawReader->GetDataSize(); // in bytes
      UInt_t *buffer = new UInt_t[totalDataWord/4];
  
      rawReader->ReadNext((UChar_t*)buffer, totalDataWord); 

      // elex info
      Int_t    buspatchId;
      UChar_t  channelId;
      UShort_t manuId;
      Char_t   parity;
      UShort_t charge; 

      // indexes
      Int_t indexDsp;
      Int_t indexBusPatch;
      Int_t index = 0;

      for(Int_t iBlock = 0; iBlock < 2 ;iBlock++){  // loop over 2 blocks
	totalBlockSize = buffer[index];
	  
	if(totalBlockSize > blankBlockSize) {        // compare block header
	  index += kBlockHeaderSize;

	  for(Int_t iDsp = 0; iDsp < iDspMax ;iDsp++){   //DSP loop

	    totalDspSize = buffer[index];
	    indexDsp = index;

	    blankDspSize =  kDspHeaderSize + iBusPerDSP[iDsp]*kBusPatchHeaderSize; // no data just header

	    if(totalDspSize > blankDspSize) {       // Compare DSP Header
	      index += kDspHeaderSize;
		
	      for(Int_t iBusPatch = 0; iBusPatch < iBusPerDSP[iDsp]; iBusPatch++) {  

		//copy buffer into header structure
		memcpy(subEventTracker->GetBusPatchHeader(), &buffer[index], kBusPatchHeaderSize*4);

		totalBusPatchSize = subEventTracker->GetTotalLength();
		buspatchId        = subEventTracker->GetBusPatchId();
		indexBusPatch     = index;
		

		if(totalBusPatchSize > kBusPatchHeaderSize) {    //Check Buspatch header, not empty events

		  index   += kBusPatchHeaderSize;
		  dataSize = subEventTracker->GetLength();

		  if(dataSize>0) { // check data present

		    //copy buffer into data structure
		    memcpy(subEventTracker->GetData(), &buffer[index], dataSize*4); 
		    index += dataSize;

		    for(Int_t iData = 0; iData < dataSize; iData++) {

		      // digits info
		      parity    = subEventTracker->GetParity(iData); // test later for parity
		      manuId    = subEventTracker->GetManuId(iData);
		      channelId = subEventTracker->GetChannelId(iData);
		      charge    = subEventTracker->GetCharge(iData);
		      // set charge
		      digit->SetSignal(charge);

		      Int_t error = GetMapping(buspatchId,manuId,channelId,digit); // Get Back the hits at pads
		      if (error) continue;

		      // debugging 
		      if (AliLog::GetGlobalDebugLevel() == 3) {
			Int_t padX  = digit->PadX();
			Int_t padY  = digit->PadY();
			Int_t iCath = digit->Cathode();  
			Int_t idDE  = digit->DetElemId();

			AliDebug(1,Form("output  IdDE %d busPatchid %d PadX %d PadY %d iCath %d \n", 
				      idDE, buspatchId, padX, padY, iCath));
		
			AliDebug(3,Form("idDE %d Padx %d Pady %d, Cath %d, charge %d",idDE, padX, padY, iCath, charge));
		      }

		      // fill digits
		      fMUONData->AddDigit(iDDL/2, *digit);

		    } // data loop
		  } // dataSize test
		} // testing buspatch

		index = indexBusPatch + totalBusPatchSize;

	      }  //buspatch loop
		
	    }  // dsp test

	    index = indexDsp + totalDspSize;
	      
	  }  // dsp loop

	}   //block test

	index = totalBlockSize;

      }  //block loop

      delete[] buffer;
    } //loop checking the header size of DDL

    //delete rawReader;
  } // DDL loop


  delete subEventTracker;
  delete digit;

  return kTRUE;
}

//____________________________________________________________________
Int_t AliMUONRawReader::GetMapping(Int_t busPatchId, UShort_t manuId, 
					 UChar_t channelId, AliMUONDigit* digit )
{

 // mapping  for tracker

  // getting DE from buspatch
  Int_t  idDE = fBusPatchManager->GetDEfromBus(busPatchId);
  AliDebug(3,Form("idDE: %d busPatchId %d\n", idDE, busPatchId));

  // segmentation
  Int_t iCath;
  Int_t iCath1 = 0;
  Int_t iCath2 = 1;

  if (idDE < 500) { // should use GetDirection somehow (ChF)
    if ( ((idDE % 100) % 2) != 0 ) {
      iCath1 = 1;
      iCath2 = 0;
    }
  }

  iCath = (manuId > 1000) ? iCath2 : iCath1;

  if (manuId > 1000) manuId -= 1000; // back to normal manuId

  // Could the above logic be simplified ???
  //AliMpVSegmentation* seg = AliMUONSegmentationManager::Segmentation(idDE, plane);
  AliMpVSegmentation* seg = fSegFactory->CreateMpSegmentation(idDE, iCath);  
  AliMpPad pad = seg->PadByLocation(AliMpIntPair(manuId,(Int_t)channelId),kTRUE);

  if(!pad.IsValid()){
    AliWarning(Form("No pad for idDE: %d, busPatchId %d, manuId: %d, channelId: %d\n",
		  idDE, busPatchId, manuId, channelId));
    return kTRUE;
  } // return error

  // Getting padX
  Int_t padX = pad.GetIndices().GetFirst();

 // Getting padY
  Int_t padY = pad.GetIndices().GetSecond();
     
  if (idDE >= 500) { // Since in AliMpSlat pads begin at (0,0) 
    padX++;         // while in AliMUONSt345Seg. they begin at (1,1)
    padY++;
  }
  // storing into digits
  digit->SetPadX(padX);
  digit->SetPadY(padY);
  digit->SetCathode(iCath);
  digit->SetDetElemId(idDE);

  AliDebug(3,Form("idDE: %d, busPatchId %d, manuId: %d, channelId: %d, padx: %d pady %d\n",
		  idDE, busPatchId, manuId, channelId, padX, padY));
  return kFALSE;
}

//____________________________________________________________________
Int_t AliMUONRawReader::ReadTriggerDDL(AliRawReader* rawReader)
{

  // reading DDL for trigger

  AliMUONSubEventTrigger* subEventTrigger = new AliMUONSubEventTrigger();
  AliMUONScalerEventTrigger* scalerEvent = 0x0;

  AliMUONGlobalTrigger* globalTrigger = 0x0;
  AliMUONLocalTrigger* localTrigger = new  AliMUONLocalTrigger();


  //Int_t kDDLHeaderSize = fDDLTrigger->GetHeaderSize();    
  // we dont need this, as size of ddl data is same for triger and no trigger

  Int_t kDDLEnhanceHeaderSize = fDDLTrigger->GetHeaderLength(); 
  Int_t kRegHeaderSize        = subEventTrigger->GetRegHeaderLength() ;

  Int_t loCircuit, loStripX, loDev, loStripY, loLpt, loHpt;
  Char_t loDecision; 

  UShort_t x1Pattern, x2Pattern, x3Pattern, x4Pattern;
  UShort_t y1Pattern, y2Pattern, y3Pattern, y4Pattern;

  // loop over the two ddl's
  for(Int_t iDDL = 0; iDDL < 2; iDDL++) { //DDL loop

    rawReader->Select(0XA,iDDL,iDDL);  //Select the DDL file to be read  

    rawReader->ReadHeader();

    Int_t totalDataWord = rawReader->GetDataSize(); // in bytes
    UInt_t *buffer = new UInt_t[totalDataWord/4];

    rawReader->ReadNext((UChar_t*)buffer, totalDataWord); 
  
    Int_t index = 0;
    fDDLTrigger->SetDDLWord(buffer[index++]);

    if(fDDLTrigger->GetEventType() == 2) {
      fScalerEvent = kTRUE;
      scalerEvent = new  AliMUONScalerEventTrigger();
    } else
      fScalerEvent = kFALSE;

    if(fScalerEvent) {
      // 6 DARC scaler words
      memcpy(scalerEvent->GetDarcScalers(), &buffer[index], scalerEvent->GetDarcScalerLength()*4);
      index += scalerEvent->GetDarcScalerLength();
    }

    // 4 words of global board input + Global board output
    memcpy(fDDLTrigger->GetGlobalInput(), &buffer[index], (kDDLEnhanceHeaderSize-1)*4); 
    index += kDDLEnhanceHeaderSize - 1; // kind tricky cos scaler info in-between Darc header

    if(fScalerEvent) {
      // 10 Global scaler words
      memcpy(scalerEvent->GetGlobalScalers(), &buffer[index], scalerEvent->GetGlobalScalerLength()*4);
      index += scalerEvent->GetGlobalScalerLength();
    }

    // fill global trigger information
    globalTrigger = GetGlobalTriggerPattern(fDDLTrigger->GetGlobalOuput());
    fMUONData->AddGlobalTrigger(*globalTrigger);
    
    // 8 regional boards
    for (Int_t iReg = 0; iReg < 8; iReg++) {           //loop over regeonal card

      memcpy(subEventTrigger->GetRegHeader(), &buffer[index], kRegHeaderSize*4);
      index += kRegHeaderSize;

     // 11 regional scaler word
      if(fScalerEvent) {
	memcpy(scalerEvent->GetRegScalers(), &buffer[index], scalerEvent->GetRegScalerLength()*4);
	index += scalerEvent->GetRegScalerLength();
      }

      // 16 local cards per regional board
      for (Int_t iLoc = 0; iLoc < 16; iLoc++) {         //loop over local card
	  
	Int_t iLocIndex = index;

	// 5 word trigger information
	for(Int_t iData = 0; iData < 5 ;iData++ ){
	  subEventTrigger->SetLocalData(buffer[index++],5*iLoc+iData);   //read local data
	}

	if(buffer[iLocIndex] > 0) {

	  loCircuit = (Int_t)subEventTrigger->GetLocalId(iLoc)+ 16*iReg + 128*iDDL; 
	  loStripX  = (Int_t)subEventTrigger->GetXPos(iLoc);
	  loStripY  = (Int_t)subEventTrigger->GetYPos(iLoc);
	  loDev     = (Int_t)subEventTrigger->GetXDev(iLoc);
	    
	  // fill local trigger
	  localTrigger->SetLoCircuit(loCircuit);
	  localTrigger->SetLoStripX(loStripX );
	  localTrigger->SetLoStripY(loStripY);
	  localTrigger->SetLoDev(loDev);

	  loDecision = subEventTrigger->GetLocalDec(iLoc);
	  loLpt =  loDecision       & 0x3;
	  loHpt = (loDecision >> 2) & 0x3; 
	    
	  // fill local trigger
	  localTrigger->SetLoLpt(loLpt);
	  localTrigger->SetLoHpt(loHpt);

	  //getting pattern from subvent
	  x1Pattern = subEventTrigger->GetX1(iLoc);
	  x2Pattern = subEventTrigger->GetX2(iLoc);
	  x3Pattern = subEventTrigger->GetX3(iLoc);
	  x4Pattern = subEventTrigger->GetX4(iLoc);
	    
	  y1Pattern = subEventTrigger->GetY1(iLoc);
	  y2Pattern = subEventTrigger->GetY2(iLoc);
	  y3Pattern = subEventTrigger->GetY3(iLoc);
	  y4Pattern = subEventTrigger->GetY4(iLoc);

	  // fill local trigger
	  localTrigger->SetX1Pattern(x1Pattern);
	  localTrigger->SetX2Pattern(x2Pattern);
	  localTrigger->SetX3Pattern(x3Pattern);
	  localTrigger->SetX4Pattern(x4Pattern);

	  localTrigger->SetY1Pattern(y1Pattern);
	  localTrigger->SetY2Pattern(y2Pattern);
	  localTrigger->SetY3Pattern(y3Pattern);
	  localTrigger->SetY4Pattern(y4Pattern);
	  fMUONData->AddLocalTrigger(*localTrigger);

	} // if buffer[] > 0

	// 45 regional scaler word
	if(fScalerEvent) {
	  memcpy(scalerEvent->GetLocalScalers(), &buffer[index], scalerEvent->GetLocalScalerLength()*4);
	  index += scalerEvent->GetLocalScalerLength();
	}
	  
      } // local card loop
	
    } // regional card loop
      
    delete [] buffer;
  } // DDL loop

  delete subEventTrigger;
  delete globalTrigger;
  delete localTrigger;

  if(fScalerEvent)
    delete scalerEvent;

  return kTRUE;

}
//____________________________________________________________________
AliMUONGlobalTrigger* AliMUONRawReader::GetGlobalTriggerPattern(Int_t gloTrigPat) const
{
  // global trigger pattern calculation

  Int_t globalSinglePlus[3];  // tot num of single plus
  Int_t globalSingleMinus[3]; // tot num of single minus
  Int_t globalSingleUndef[3]; // tot num of single undefined
  Int_t globalPairUnlike[3];  // tot num of unlike-sign pairs
  Int_t globalPairLike[3];    // tot num of like-sign pairs


  for (Int_t i = 0; i < 3; i++) {

    globalSinglePlus[i]  = (gloTrigPat >> i) & 0x1;
    globalSingleMinus[i] = (gloTrigPat >> (i+3)) & 0x1;
    globalSingleUndef[i] = (gloTrigPat >> (i+6)) & 0x1;
    globalPairUnlike[i]  = (gloTrigPat >> (i+9)) & 0x1;
    globalPairLike[i]    = (gloTrigPat >> (i+12)) & 0x1;
  }

  return (new AliMUONGlobalTrigger(globalSinglePlus, globalSingleMinus,
				   globalSingleUndef, globalPairUnlike, 
				   globalPairLike));  

}

