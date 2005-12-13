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
// MUON Raw Data generator and reader in ALICE-MUON
// This class version 3 (further details could be found in Alice-note)
//
// Implemented non-constant buspatch numbers for tracking
// with correct DDL id (first guess)
// (Ch. Finck, dec 2005)
//
// Digits2Raw:
// Generates raw data for MUON tracker and finally for trigger
// Using real mapping (inverse) for tracker
// For trigger there is no mapping (mapping could be found in AliMUONTriggerCircuit)
// Ch. Finck july 04
//
// Raw2Digits:
// Using real mapping  for tracker
// Indranil Das (Adapted for runloader: Ch. Finck) july 05
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
#include "AliMUONRawData.h"
#include "AliMUONDigit.h"

#include "AliMUONConstants.h"
#include "AliMUONData.h"

#include "AliMUONSubEventTrigger.h"
#include "AliMUONDDLTracker.h"
#include "AliMUONDDLTrigger.h"

#include "AliMUONLocalTrigger.h"
#include "AliMUONGlobalTrigger.h"

#include "AliMUONGeometrySegmentation.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryStore.h"
#include "AliMUONSegmentationManager.h"
#include "AliMpPlaneType.h"
#include "AliMpVSegmentation.h"
#include "AliMpHelper.h"
#include "AliMpPad.h"


ClassImp(AliMUONRawData) // Class implementation in ROOT context
//__________________________________________________________________________
AliMUONRawData::AliMUONRawData(AliLoader* loader)
  : TObject()
{
  // Standard Constructor
 
  // initialize loader's
  fLoader = loader;

  // initialize container
  fMUONData  = new AliMUONData(fLoader,"MUON","MUON");

  // initialize array
  fSubEventArray = new TClonesArray("AliMUONSubEventTracker",1000);


  // ddl pointer
  fDDLTracker = new AliMUONDDLTracker();
  fDDLTrigger = new AliMUONDDLTrigger();

  ReadBusPatchFile();

}

//__________________________________________________________________________
AliMUONRawData::AliMUONRawData()
  : TObject(),
    fMUONData(0),
    fLoader(0),
    fDDLTracker(0),
    fDDLTrigger(0)
{
  // Default Constructor
  fFile[0] = fFile[1] = 0x0;
  
}

//_______________________________________________________________________
AliMUONRawData::AliMUONRawData (const AliMUONRawData& rhs)
  : TObject(rhs)
{
// Protected copy constructor

  AliFatal("Not implemented.");
}

//_______________________________________________________________________
AliMUONRawData & 
AliMUONRawData::operator=(const AliMUONRawData& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}

//__________________________________________________________________________
AliMUONRawData::~AliMUONRawData(void)
{
  if (fMUONData)
    delete fMUONData;
  if (fSubEventArray)
    fSubEventArray->Delete(); //using delete cos allocating memory in copy ctor.

  if (fDDLTracker)
    delete fDDLTracker;
  if (fDDLTrigger)
    delete fDDLTrigger;

  fDetElemIdToBusPatch.Delete();
  fBusPatchToDetElem.Delete();
  fBusPatchToDDL.Delete();

  return;
}
//____________________________________________________________________
Int_t AliMUONRawData::Digits2Raw()
{
 // convert digits of the current event to raw data

  Int_t idDDL;
  Char_t name[20];

  fLoader->LoadDigits("READ");

  fMUONData->SetTreeAddress("D,GLT");


  // tracking chambers

  for (Int_t ich = 0; ich < AliMUONConstants::NTrackingCh(); ich++) {
 
    // open files
    idDDL = ich * 2  + 0x900; // official number for MUON
    sprintf(name, "MUON_%d.ddl",idDDL);
    fFile[0] = fopen(name,"w");

    idDDL = (ich * 2) + 1 + 0x900;
    sprintf(name, "MUON_%d.ddl",idDDL);
    fFile[1] = fopen(name,"w");

    WriteTrackerDDL(ich);
  
    // reset and close
    fclose(fFile[0]);
    fclose(fFile[1]);
    fMUONData->ResetDigits();
  }
 
  // trigger chambers
 
  // open files
  idDDL = 0xA00;// official number for MUTR
  sprintf(name, "MUTR_%d.ddl",idDDL);
  fFile[0] = fopen(name,"w");

  idDDL = 0xA00 + 1;
  sprintf(name, "MUTR_%d.ddl",idDDL);
  fFile[1] = fopen(name,"w");

  WriteTriggerDDL();
  
  // reset and close
  fclose(fFile[0]);
  fclose(fFile[1]);
  fMUONData->ResetTrigger();
  
  fLoader->UnloadDigits();

  return kTRUE;
}
//____________________________________________________________________
Int_t AliMUONRawData::WriteTrackerDDL(Int_t iCh)
{
  // writing DDL for tracker
  // used inverse mapping

  // resets
  TClonesArray* muonDigits = 0;
  fSubEventArray->Clear();

  //
  TArrayI nbInBus;

  nbInBus.Set(5000);

  nbInBus.Reset();

  // DDL header
  AliRawDataHeader header = fDDLTracker->GetHeader();
  Int_t headerSize = fDDLTracker->GetHeaderSize();

  // DDL event one per half chamber
  AliMUONSubEventTracker* subEvent;

  // data format
  Char_t parity = 0x4;
  UShort_t manuId = 0;
  UChar_t channelId = 0;
  UShort_t charge = 0;
  Int_t busPatchId = 0;

  UInt_t word;
  Int_t nEntries = 0;
  Int_t* buffer = 0;
  Int_t index;
  Int_t indexDsp;
  Int_t indexBlk;
  Int_t padX;
  Int_t padY;
  Int_t cathode = 0;
  Int_t detElemId;
  Int_t nDigits;

  const AliMUONDigit* digit;

  AliDebug(3, Form("WriteDDL chamber %d\n", iCh+1));

  // getting digits
  fMUONData->ResetDigits();
  fMUONData->GetDigits();
  muonDigits = fMUONData->Digits(iCh);

  nDigits = muonDigits->GetEntriesFast();
  AliDebug(3,Form("ndigits = %d\n",nDigits));
 
  // loop over digit
  for (Int_t idig = 0; idig < nDigits; idig++) {

    digit = (AliMUONDigit*) muonDigits->UncheckedAt(idig);

    padX = digit->PadX();
    padY = digit->PadY();
    charge = digit->Signal();
    charge &= 0xFFF;
    cathode = digit->Cathode();
    detElemId = digit->DetElemId();

    // inverse mapping
    Int_t error = GetInvMapping(digit, busPatchId, manuId, channelId);
    if (error) continue;

    AliDebug(3,Form("input  IdDE %d busPatchId %d PadX %d PadY %d iCath %d \n", 
		    detElemId, busPatchId, padX, padY, cathode));

    AliDebug(3,Form("busPatchId %d, manuId %d channelId %d\n", busPatchId, manuId, channelId ));

    //packing word
    AliBitPacking::PackWord((UInt_t)parity,word,29,31);
    AliBitPacking::PackWord((UInt_t)manuId,word,18,28);
    AliBitPacking::PackWord((UInt_t)channelId,word,12,17);
    AliBitPacking::PackWord((UInt_t)charge,word,0,11);

    // set sub Event
    subEvent = new AliMUONSubEventTracker();
    subEvent->AddData(word);
    subEvent->SetBusPatchId(busPatchId);
       
    // storing the number of identical buspatches
    nbInBus[busPatchId]++;
    AddData(subEvent);
   
    delete subEvent;
  }

  // sorting by buspatch
  fSubEventArray->Sort();

  // gather datas from same bus patch
  nEntries = fSubEventArray->GetEntriesFast();

  for (Int_t i = 0; i < nEntries; i++) {
    AliMUONSubEventTracker* temp = (AliMUONSubEventTracker*)fSubEventArray->At(i);
    busPatchId = temp->GetBusPatchId();

    // add bus patch header, length and total length managed by subevent class
    temp->SetTriggerWord(0xdeadbeef);
    for (Int_t j = 0; j < nbInBus[busPatchId]-1; j++) {
      AliMUONSubEventTracker* temp1 =  (AliMUONSubEventTracker*)fSubEventArray->At(++i);
      temp->AddData(temp1->GetData(0));
      fSubEventArray->RemoveAt(i) ;
    }
  }
  fSubEventArray->Compress();

  if (AliLog::GetGlobalDebugLevel() == 3) {
    nEntries = fSubEventArray->GetEntriesFast();
    for (Int_t i = 0; i < nEntries; i++) {
      AliMUONSubEventTracker* temp =  (AliMUONSubEventTracker*)fSubEventArray->At(i);
      printf("busPatchid back %d\n",temp->GetBusPatchId());
      for (Int_t j = 0; j < temp->GetLength(); j++) {
	printf("manuId back %d, ",temp->GetManuId(j));
	printf("channelId back %d, ",temp->GetChannelId(j));
	printf("charge back %d\n",temp->GetCharge(j));
      }
    }
    printf("\n");
  }
  
  // getting info for the number of buspatches
  Int_t iBusPatch;
  Int_t length;
  Int_t iBusPerDSP[5];//number of bus patches per DSP
  Int_t iDspMax; //number max of DSP per block
 
  Int_t iFile = 0;
  GetDspInfo(iCh, iDspMax, iBusPerDSP);

  TArrayI* vec = GetBusfromDE((iCh+1)*100);

  Int_t iBus0AtCh = vec->At(0); //get first bus patch id for a given ich
	
  AliDebug(3,Form("iBus0AtCh %d", iBus0AtCh));

  iBusPatch = iBus0AtCh - 1; // starting point for each chamber

  // nEntries = fSubEventArray->GetEntriesFast();
  AliMUONSubEventTracker* temp = 0x0;

  // open DDL file, on per 1/2 chamber
  for (Int_t iDDL = 0; iDDL < 2; iDDL++) {
    

    // filling buffer
    buffer = new Int_t [(2048+24)*50]; // 24 words in average for one buspatch and 2048 manu info at most

    indexBlk = 0;
    indexDsp = 0;
    index = 0;

    // two blocks A and B per DDL
    for (Int_t iBlock = 0; iBlock < 2; iBlock++) {

      // block header
      length = fDDLTracker->GetBlkHeaderLength();
      memcpy(&buffer[index],fDDLTracker->GetBlkHeader(),length*4);
      indexBlk = index;
      index += length; 

      // 5 DSP's max per block
      for (Int_t iDsp = 0; iDsp < iDspMax; iDsp++) {

	// DSP header
	length = fDDLTracker->GetDspHeaderLength();
	memcpy(&buffer[index],fDDLTracker->GetDspHeader(),length*4);
	indexDsp = index;
	index += length; 

	// 5 buspatches max per DSP
	for (Int_t i = 0; i < iBusPerDSP[iDsp]; i++) {

	  iBusPatch ++;
	  if ((fBusPatchToDDL(iBusPatch) % 2) == 1) // comparing to DDL file
	    iFile = 0;
	  else
	    iFile = 1;

	  AliDebug(3,Form("iCh %d iDDL %d iBlock %d iDsp %d busPatchId %d", iCh, iDDL, iBlock, iDsp, iBusPatch));

	  nEntries = fSubEventArray->GetEntriesFast();

	  for (Int_t iEntries = 0; iEntries < nEntries; iEntries++) { // method "bourrique"...
	    temp = (AliMUONSubEventTracker*)fSubEventArray->At(iEntries);
	    busPatchId = temp->GetBusPatchId();
	    if (busPatchId == iBusPatch) break;
	    busPatchId = -1;
	    AliDebug(3,Form("busPatchId %d", temp->GetBusPatchId()));
	  } 
	 
	  // check if buspatchid has digit
	  if (busPatchId != -1) {
	    // add bus patch structure
	    length = temp->GetHeaderLength();
	    memcpy(&buffer[index],temp->GetAddress(),length*4);
	    index += length;
	    for (Int_t j = 0; j < temp->GetLength(); j++) {
	      buffer[index++] =  temp->GetData(j);
	      AliDebug(3,Form("busPatchId %d, manuId %d channelId %d\n", temp->GetBusPatchId(), 
			      temp->GetManuId(j), temp->GetChannelId(j) ));
	    }
	    //	      fSubEventArray->RemoveAt(iEntries);
	    //	      fSubEventArray->Compress();
	  } else {
	    // writting anyhow buspatch structure (empty ones)
	    buffer[index++] = 4; // total length
	    buffer[index++] = 0; // raw data length
	    buffer[index++] = iBusPatch; // bus patch
	    buffer[index++] = 0xdeadbeef; // trigger word
	  }
	} // bus patch
	buffer[indexDsp] = index - indexDsp; // dsp length
	buffer[indexDsp+1] = index - indexDsp - fDDLTracker->GetDspHeaderLength();
	if ((index - indexDsp) % 2 == 0)
	  buffer[indexDsp+7] = 0;
	else
	  buffer[indexDsp+7] = 1;
      } // dsp
      buffer[indexBlk] = index - indexBlk; // block length
      buffer[indexBlk+1] = index - indexBlk - fDDLTracker->GetBlkHeaderLength();
    }
    
    //writting onto disk
    // write DDL 1 & 2
    header.fSize = (index + headerSize) * 4;// total length in bytes
    fwrite((char*)(&header),headerSize*4,1,fFile[iFile]);
    fwrite(buffer,sizeof(int),index,fFile[iFile]);
   
    delete[] buffer;
  }

  return kTRUE;
}
//____________________________________________________________________
Int_t AliMUONRawData::WriteTriggerDDL()
{

 // DDL event one per half chamber
  AliMUONSubEventTrigger* subEvent = 0x0;


  // stored local id number 
  TArrayI isFired(256);
  isFired.Reset();


 // DDL header
  AliRawDataHeader header = fDDLTrigger->GetHeader();
  Int_t headerSize = fDDLTrigger->GetHeaderSize();
  Int_t length;
  TClonesArray* localTrigger;
  TClonesArray* globalTrigger;
  AliMUONGlobalTrigger* gloTrg;
  AliMUONLocalTrigger* locTrg = 0x0;

  // getting information from trigger
  fMUONData->GetTriggerD();

  // global trigger for trigger pattern
  globalTrigger = fMUONData->GlobalTrigger(); 
  gloTrg = (AliMUONGlobalTrigger*)globalTrigger->UncheckedAt(0);
  Int_t gloTrigPat = GetGlobalTriggerPattern(gloTrg);

  // local trigger 
  localTrigger = fMUONData->LocalTrigger();    

  UInt_t word;
  Int_t* buffer = 0;
  Int_t index;
  Int_t iEntries = 0;
  Int_t iLocCard, locCard;
  Char_t locDec, trigY, posY, posX,regOut;
  Int_t devX;
  Int_t version = 1; // software version
  Int_t eventType = 1; // trigger type: 1 for physics ?
  Int_t serialNb = 0xF; // serial nb of card: all bits on for the moment
  Int_t globalFlag = 1; // set to 2 if global info present in DDL else set to 1

  Int_t nEntries = (Int_t) (localTrigger->GetEntries());// 234 local cards
  // stored the local card id that's fired
  for (Int_t i = 0; i <  nEntries; i++) {
    locTrg = (AliMUONLocalTrigger*)localTrigger->At(i);
    isFired[locTrg->LoCircuit()] = 1; // storing local boards with informations
  }

  if (!nEntries)
    AliError("No Trigger information available");

  buffer = new Int_t [672]; // [16(local)*5 words + 3 words]*8(reg) + 8 words = 672

  // open DDL file, on per 1/2 chamber
  for (Int_t iDDL = 0; iDDL < 2; iDDL++) {
    
    index = 0; 

    // DDL enhanced header
    word = 0;
    AliBitPacking::PackWord((UInt_t)iDDL+1,word,28,31); //see AliMUONDDLTrigger.h for details
    AliBitPacking::PackWord((UInt_t)serialNb,word,24,27);
    AliBitPacking::PackWord((UInt_t)version,word,16,23);
    AliBitPacking::PackWord((UInt_t)eventType,word,12,15);

    if (iDDL == 0) // suppose global info in DDL one
      globalFlag = 2;
    else 
      globalFlag = 1;
    AliBitPacking::PackWord((UInt_t)globalFlag,word,8,11);
    fDDLTrigger->SetDDLWord(word);

    if (iDDL == 0)
      fDDLTrigger->SetGlobalOutput(gloTrigPat);// no global input for the moment....
    else 
      fDDLTrigger->SetGlobalOutput(0);
    length = fDDLTrigger->GetHeaderLength(); 
    memcpy(&buffer[index],fDDLTrigger->GetEnhancedHeader(),length*4);
    index += length; 

    // 8 regional cards per DDL
    for (Int_t iReg = 0; iReg < 8; iReg++) {

      subEvent = new AliMUONSubEventTrigger();

      // Regional card header
      word = 0;
      regOut  = 0;
      AliBitPacking::PackWord((UInt_t)serialNb,word,24,28); //see  AliMUONSubEventTrigger.h for details
      AliBitPacking::PackWord((UInt_t)version,word,16,23);
      AliBitPacking::PackWord((UInt_t)iReg,word,12,15);
      AliBitPacking::PackWord((UInt_t)regOut,word,0,7); // whenever regional output will be implemented

      subEvent->SetRegWord(word);
      memcpy(&buffer[index++],subEvent->GetAddress(),4);

      buffer[index++] = 0;// 2 words of regional input
      buffer[index++] = 0;

      // 16 local card per regional board
      for (Int_t iLoc = 0; iLoc < 16; iLoc++) {

	iLocCard = iLoc + iReg*16 + iDDL*128;

	if (isFired[iLocCard]) {
	  locTrg = (AliMUONLocalTrigger*)localTrigger->At(iEntries);
	  locCard = locTrg->LoCircuit();
	  locDec  = locTrg->GetLoDecision();
	  trigY = 0;
	  posY = locTrg->LoStripY();
	  posX = locTrg->LoStripX();
	  devX = locTrg->LoDev();
	  AliDebug(4,Form("loctrg %d, posX %d, posY %d, devX %d\n", 
			  locTrg->LoCircuit(),locTrg->LoStripX(),locTrg->LoStripY(),locTrg->LoDev()));
	} else { //no trigger (see PRR chpt 3.4)
	  locCard = -1;
	  locDec = 0;
	  trigY = 1;
	  posY = 15;
	  posX = 0;
	  devX = 0x8;
	}

	//packing word
	word = 0;
	AliBitPacking::PackWord((UInt_t)(iLocCard % 16),word,19,22); //card id number in crate
	AliBitPacking::PackWord((UInt_t)locDec,word,15,18);
	AliBitPacking::PackWord((UInt_t)trigY,word,14,14);
	AliBitPacking::PackWord((UInt_t)posY,word,10,13);
	AliBitPacking::PackWord((UInt_t)devX,word,5,9);
 	AliBitPacking::PackWord((UInt_t)posX,word,0,4);

	if (locCard == iLocCard) {
	  // add local cards structure
	  buffer[index++] = (locTrg->GetX1Pattern() | (locTrg->GetX2Pattern() << 16));
	  buffer[index++] = (locTrg->GetX3Pattern() | (locTrg->GetX4Pattern() << 16));
	  buffer[index++] = (locTrg->GetY1Pattern() | (locTrg->GetY2Pattern() << 16));
	  buffer[index++] = (locTrg->GetY3Pattern() | (locTrg->GetY4Pattern() << 16));
	  buffer[index++] = (Int_t)word; // data word
	  if (iEntries < nEntries-1)
	    iEntries++;
	} else {
	  buffer[index++] = 0; // 4 words for x1, x2, y1, y2
	  buffer[index++] = 0; 
	  buffer[index++] = 0; 
	  buffer[index++] = 0; 
	  buffer[index++] = (Int_t)word; // data word

	}
      } // local card 

      delete subEvent;	

    } // Regional card
    
    buffer[index++] = fDDLTrigger->GetEoD(); // End of DDL word
    buffer[index++] = fDDLTrigger->GetEoD(); // End of DDL word for 64 bits transfer purpose

    // writting onto disk
    // write DDL 1
    header.fSize = (index + headerSize) * 4;// total length in bytes
    fwrite((char*)(&header),headerSize*4,1,fFile[iDDL]);
    fwrite(buffer,sizeof(int),index,fFile[iDDL]);
  
  }
  delete[] buffer;

  return kTRUE;
}

//____________________________________________________________________
Int_t AliMUONRawData::GetInvMapping(const AliMUONDigit* digit,
				     Int_t &busPatchId, UShort_t &manuId, UChar_t &channelId)
{

  // Inverse mapping for tracker

  // information from digits
  Int_t iCath = digit->Cathode();
  Int_t idDE  = digit->DetElemId();
  Int_t padX  = digit->PadX();
  Int_t padY  = digit->PadY();

  if (idDE >= 500) { // Since in AliMpSlat pads begin at (0,0) 
    padX--;         // while in AliMUONSt345Seg. they begin at (1,1)
    padY--;
  }

  // segmentation
  AliMpPlaneType plane;
  AliMpPlaneType plane1 = kBendingPlane;
  AliMpPlaneType plane2 = kNonBendingPlane;

  if (idDE < 500) { // should use GetDirection somehow (ChF)
    if ( ((idDE % 100) % 2) != 0 ) {
      plane1 = kNonBendingPlane;
      plane2 = kBendingPlane;
    }
  }
  // station 345 bending == cath0 for the moment
   plane = (iCath == 0) ? plane1 : plane2;

  AliMpVSegmentation* seg = AliMUONSegmentationManager::Segmentation(idDE, plane);
  AliMpPad pad = seg->PadByIndices(AliMpIntPair(padX,padY),kTRUE);

  if(!pad.IsValid()) {
     AliWarning(Form("No elec. for idDE: %d, padx: %d pady %d, charge: %d\n",
		  idDE, digit->PadX(), digit->PadY(), digit->Signal()));
    return kTRUE;
  }

  // Getting Manu id
  manuId = pad.GetLocation().GetFirst();
  manuId &= 0x7FF; // 11 bits 

  // Getting channel id
  channelId =  pad.GetLocation().GetSecond();
  channelId &= 0x3F; // 6 bits

  // Getting buspatch id
  TArrayI* vec = GetBusfromDE(idDE);
  Int_t pos;

  if (idDE < 500) { // station 1 & 2
    // set 32 manus for one bus patch ? (ChF)
    pos = manuId/32;
  } else {
    // offset of 100 in manuId for following bus patch
    pos = manuId/100;
  }

 //  if (pos >(int_t) vec.size())
//     AliWarning("pos greater than size\n");
  busPatchId = vec->At(pos);

  if (plane ==  kNonBendingPlane) // for Non-Bending manuid+= 1000;
    manuId += 1000; // tmp solution til one finds something better  (ChF)
  
  AliDebug(3,Form("idDE: %d, busPatchId %d, manuId: %d, channelId:%d\n",
		  idDE, busPatchId, manuId, channelId));

  AliDebug(3,Form("idDE: %d, busPatchId %d, manuId: %d, channelId: %d, padx: %d pady: %d, charge: %d\n",
		  idDE, busPatchId, manuId, channelId, digit->PadX(), digit->PadY(), digit->Signal()));

  return kFALSE; // no error
}

//____________________________________________________________________
Int_t AliMUONRawData::GetGlobalTriggerPattern(const AliMUONGlobalTrigger* gloTrg) const
{
  // global trigger pattern calculation

  Int_t gloTrigPat = 0;

  if (gloTrg->SinglePlusLpt())  gloTrigPat|= 0x1;
  if (gloTrg->SinglePlusHpt())  gloTrigPat|= 0x2;
  if (gloTrg->SinglePlusApt())  gloTrigPat|= 0x4;
 
  if (gloTrg->SingleMinusLpt()) gloTrigPat|= 0x8;
  if (gloTrg->SingleMinusHpt()) gloTrigPat|= 0x10;
  if (gloTrg->SingleMinusApt()) gloTrigPat|= 0x20;
 
  if (gloTrg->SingleUndefLpt()) gloTrigPat|= 0x40;
  if (gloTrg->SingleUndefHpt()) gloTrigPat|= 0x80;
  if (gloTrg->SingleUndefApt()) gloTrigPat|= 0x100;
 
  if (gloTrg->PairUnlikeLpt())  gloTrigPat|= 0x200;
  if (gloTrg->PairUnlikeHpt())  gloTrigPat|= 0x400;
  if (gloTrg->PairUnlikeApt())  gloTrigPat|= 0x800;

  if (gloTrg->PairLikeLpt())    gloTrigPat|= 0x1000;
  if (gloTrg->PairLikeHpt())    gloTrigPat|= 0x2000;
  if (gloTrg->PairLikeApt())    gloTrigPat|= 0x4000;

  return gloTrigPat;
}

//____________________________________________________________________
Int_t AliMUONRawData::Raw2Digits(AliRawReader* rawReader)
{

  // generate digits
  ReadTrackerDDL(rawReader);

  // generate trigger
  ReadTriggerDDL(rawReader);

  return kTRUE;

}

//____________________________________________________________________
Int_t AliMUONRawData::ReadTrackerDDL(AliRawReader* rawReader)
{
  // reading tracker DDL
  // filling the TClonesArray in MUONData
  //

  AliMUONSubEventTracker* subEventTracker = new AliMUONSubEventTracker();
  AliMUONDigit* digit = new AliMUONDigit();


  //Read Header Size of DDL,Block,DSP and BusPatch.

  Int_t ddlHeaderSize      = fDDLTracker->GetHeaderSize();
  Int_t blockHeaderSize    = fDDLTracker->GetBlkHeaderLength();
  Int_t dspHeaderSize      = fDDLTracker->GetDspHeaderLength();
  Int_t buspatchHeaderSize = subEventTracker->GetHeaderLength();

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
    GetDspInfo(iDDL/2, iDspMax, iBusPerDSP);

    //   Each DDL is made with 2 Blocks each of which consists of 5 DSP's at most and each of DSP has at most 5 buspatches.
    //   This information is used to calculate the size of headers (DDL,Block and DSP) which has no interesting data.
    blankDDLSize   = ddlHeaderSize + 2*blockHeaderSize + 2*iDspMax*dspHeaderSize;
    blankBlockSize = blockHeaderSize + iDspMax*dspHeaderSize;

    for (Int_t i = 0; i < iDspMax; i++) {
      blankDDLSize   += 2*iBusPerDSP[i]*buspatchHeaderSize;
      blankBlockSize +=   iBusPerDSP[i]*buspatchHeaderSize;
    }

    rawReader->Select(0X9, iDDL, iDDL);  //Select the DDL file to be read  

    rawReader->ReadHeader();

    totalDDLSize = (rawReader->GetDataSize() + sizeof(AliRawDataHeader))/4; // 4 is multiplied to convert byte 2 words

    if(totalDDLSize > blankDDLSize) {      // Compare the DDL header with an empty DDL header size to read the file

      Int_t totalDataWord = rawReader->GetDataSize()/4 ;
      UInt_t *buffer = new UInt_t[totalDataWord];
      for(Int_t i = 0; i < totalDataWord; i++) { 
	UInt_t& temp = buffer[i]; 
	rawReader->ReadNextInt(temp);      // takes the whole result into buffer variable for future analysis
      }

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
	  index += blockHeaderSize;

	  for(Int_t iDsp = 0; iDsp < iDspMax ;iDsp++){   //DSP loop

	    totalDspSize = buffer[index];
	    indexDsp = index;

	    blankDspSize =  dspHeaderSize + iBusPerDSP[iDsp]*buspatchHeaderSize; // no data just header

	    if(totalDspSize > blankDspSize) {       // Compare DSP Header
	      index += dspHeaderSize;
		
	      for(Int_t iBusPatch = 0; iBusPatch < iBusPerDSP[iDsp]; iBusPatch++) {  

		totalBusPatchSize = buffer[index];
		buspatchId        = buffer[index+2];
		indexBusPatch     = index;

		if(totalBusPatchSize > buspatchHeaderSize) {    //Check Buspatch header

		  index   += buspatchHeaderSize;
		  dataSize = totalBusPatchSize - buspatchHeaderSize;

		  if(dataSize>0) { // check data present

		    for(Int_t iData = 0; iData < dataSize; iData++) {

		      subEventTracker->SetData(buffer[index++],iData);   //Set to extract data
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
Int_t AliMUONRawData::GetMapping(Int_t busPatchId, UShort_t manuId, 
					 UChar_t channelId, AliMUONDigit* digit )
{

 // mapping  for tracker

  // getting DE from buspatch
  Int_t  idDE = GetDEfromBus(busPatchId);
  AliDebug(3,Form("idDE: %d busPatchId %d\n", idDE, busPatchId));

  // segmentation
  Int_t iCath;
  Int_t iCath1 = 0;
  Int_t iCath2 = 1;

  AliMpPlaneType plane;

  if (manuId > 1000) { // again tmp solution (ChF) (+1000 for Non-Bending plane
    plane = kNonBendingPlane;
  } else {
    plane = kBendingPlane;
  }

  if (idDE < 500) { // should use GetDirection somehow (ChF)
    if ( ((idDE % 100) % 2) != 0 ) {
      iCath1 = 1;
      iCath2 = 0;
    }
  }

  iCath = (manuId > 1000) ? iCath2 : iCath1;

  if (manuId > 1000) manuId -= 1000; // back to normal manuId

  AliMpVSegmentation* seg = AliMUONSegmentationManager::Segmentation(idDE, plane);
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
Int_t AliMUONRawData::ReadTriggerDDL(AliRawReader* rawReader)
{

  // reading DDL for trigger

  AliMUONSubEventTrigger* subEventTrigger = new AliMUONSubEventTrigger();
  AliMUONGlobalTrigger* globalTrigger = 0x0;
  AliMUONLocalTrigger* localTrigger = new  AliMUONLocalTrigger();


  //Int_t ddlHeaderSize = fDDLTrigger->GetHeaderSize();    
  // we dont need this, as size of ddl data is same for triger and no trigger

  Int_t ddlEnhanceHeaderSize = fDDLTrigger->GetHeaderLength(); 
  Int_t regHeaderLength      = subEventTrigger->GetRegHeaderLength() ;

  Int_t loCircuit, loStripX, loDev, loStripY, loLpt, loHpt;
  Char_t loDecision; 

  UShort_t x1Pattern, x2Pattern, x3Pattern, x4Pattern;
  UShort_t y1Pattern, y2Pattern, y3Pattern, y4Pattern;


  // loop over the two ddl's
  for(Int_t iDDL = 0; iDDL < 2; iDDL++) { //DDL loop

    rawReader->Select(0XA,iDDL,iDDL);  //Select the DDL file to be read  

    rawReader->ReadHeader();

    Int_t totalDataWord = rawReader->GetDataSize()/4 ;
    UInt_t *buffer = new UInt_t[totalDataWord];
    for(Int_t i=0;i<totalDataWord;i++){
      UInt_t& temp = buffer[i]; 
      rawReader->ReadNextInt(temp);      // takes the whole result into buffer variable for future analysis
    }

    // rawReader->ReadNext((UChar_t*)buffer, totalDataWord);     // method is protected ????
  
    Int_t index = 0;

    // fill DDL header informations
    memcpy(fDDLTrigger->GetEnhancedHeader(), &buffer[index], ddlEnhanceHeaderSize*4); 

    // fill global trigger information
    globalTrigger = GetGlobalTriggerPattern(fDDLTrigger->GetGlobalOuput());
    fMUONData->AddGlobalTrigger(*globalTrigger);

    index += ddlEnhanceHeaderSize;

    // 8 regional boards
    for (Int_t iReg = 0; iReg < 8; iReg++) {           //loop over regeonal card


      subEventTrigger->SetRegWord(buffer[index]);      //read regional data 

      index += regHeaderLength;

      // 16 local cards per regional board
      for (Int_t iLoc = 0; iLoc < 16; iLoc++) {         //loop over local card
	  
	Int_t iLocIndex = index;

	// 5 word trigger information
	for(Int_t iData = 0; iData < 5 ;iData++ ){
	  subEventTrigger->SetLocalData(buffer[index++],5*iLoc+iData);   //read local data
	}

	if(buffer[iLocIndex] > 0) {

	  loCircuit = (Int_t)subEventTrigger->GetLocalId(iLoc)+ 16*iReg + 128*iDDL; 
	  loStripX =  (Int_t)subEventTrigger->GetXPos(iLoc);
	  loStripY = (Int_t)subEventTrigger->GetYPos(iLoc);
	  loDev = (Int_t)subEventTrigger->GetXDev(iLoc);
	    
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

	}
	  
      } // local card loop
	
    } // regional card loop
      
    delete [] buffer;
  } // DDL loop

  delete subEventTrigger;
  delete globalTrigger;
  delete localTrigger;

  return kTRUE;

}
//____________________________________________________________________
AliMUONGlobalTrigger* AliMUONRawData::GetGlobalTriggerPattern(Int_t gloTrigPat) const
{
  // global trigger pattern calculation

  Int_t globalSinglePlus[3];  // tot num of single plus
  Int_t globalSingleMinus[3]; // tot num of single minus
  Int_t globalSingleUndef[3]; // tot num of single undefined
  Int_t globalPairUnlike[3];  // tot num of unlike-sign pairs
  Int_t globalPairLike[3];    // tot num of like-sign pairs


  for (Int_t i = 0; i < 3; i++) {
    globalSinglePlus[i]  = gloTrigPat & (0x1 << i);
    globalSingleMinus[i] = gloTrigPat & (0x1 << i+3);
    globalSingleUndef[i] = gloTrigPat & (0x1 << i+6);
    globalPairUnlike[i]  = gloTrigPat & (0x1 << i+9);
    globalPairLike[i]    = gloTrigPat & (0x1 << i+12);
  }

  return (new AliMUONGlobalTrigger(globalSinglePlus, globalSingleMinus,
				   globalSingleUndef, globalPairUnlike, 
				   globalPairLike));  

}
//____________________________________________________________________
Int_t AliMUONRawData::GetDEfromBus(Int_t busPatchId)
{
  // getting DE id from bus patch
  Long_t it = fBusPatchToDetElem.GetValue(busPatchId);

 if ( it ) 
   return (Int_t)it;
 else 
   return -1;
}

//____________________________________________________________________
TArrayI*  AliMUONRawData::GetBusfromDE(Int_t idDE)
{
  // getting bus patch from DE id 

  return (TArrayI*)fDetElemIdToBusPatch.GetValue(idDE);
}
//____________________________________________________________________
Int_t AliMUONRawData::GetDDLfromBus(Int_t busPatchId)
{
  // getting DE id from bus patch
  Long_t it = fBusPatchToDDL.GetValue(busPatchId);

 if ( it ) 
   return (Int_t)it;
 else 
   return -1;
}

//____________________________________________________________________
void AliMUONRawData::GetDspInfo(Int_t iCh, Int_t& iDspMax, Int_t* iBusPerDSP)
{
  // calculates the number of DSP & buspatch per block

  Int_t iBusPerBlk = fMaxBusPerCh[iCh]/4; //per half chamber; per block

  iDspMax =  iBusPerBlk/5; //number max of DSP per block
  if (iBusPerBlk % 5 != 0)
    iDspMax += 1;
  
  for (Int_t i = 0; i < iDspMax; i++) {
    if ((iBusPerBlk -= 5) > 0) 
      iBusPerDSP[i] = 5;
    else 
      iBusPerDSP[i] = iBusPerBlk + 5;
  }
  
}
//____________________________________________________________________
void AliMUONRawData::ReadBusPatchFile()
{

  // idDE <> buspatch map
  
  // reading file
   TString dirPath = gSystem->Getenv("ALICE_ROOT");
   dirPath += "/MUON/mapping/data/"; 

   TString infile = dirPath + "DetElemIdToBusPatch.dat";

   ifstream in(infile, ios::in);
   if (!in) AliError("DetElemIdToBusPatch.dat not found.");
       
   char line[80];

   Int_t iChprev = 1;
   Int_t maxBusPatch = 0;

   while ( in.getline(line,80) ) {

      if ( line[0] == '#' ) continue;

      TString tmp(AliMpHelper::Normalize(line));

      Int_t blankPos  = tmp.First(' ');
      Int_t blankPos1 = tmp.Last(' ');

      TString sDE(tmp(0, blankPos));

      Int_t idDE = atoi(sDE.Data());
      
      if (idDE/100 != iChprev) {
	fMaxBusPerCh[iChprev-1] = maxBusPatch-iChprev*100+1;
	iChprev = idDE/100;
      }

      TString sDDL(tmp(blankPos1 + 1, tmp.Length()-blankPos1));

      Int_t iDDL = atoi(sDDL.Data());

      TString busPatch(tmp(blankPos + 1,blankPos1-blankPos-1));
      AliDebug(3,Form("idDE %d buspatch %s iDDL %d\n", idDE, busPatch.Data(), iDDL));

      TArrayI busPatchList;
      // decoding range of buspatch
      AliMpHelper::DecodeName(busPatch,';',busPatchList);
      
      // filling buspatch -> idDE
      for (Int_t i = 0; i < busPatchList.GetSize(); i++) {
	fBusPatchToDetElem.Add((Long_t)busPatchList[i],(Long_t)idDE);
	fBusPatchToDDL.Add((Long_t)busPatchList[i],(Long_t)iDDL);
	maxBusPatch = busPatchList[i];
      }
   
      // filling idDE -> buspatch list (vector)
      fDetElemIdToBusPatch.Add((Long_t)idDE, (Long_t)(new TArrayI(busPatchList))); 

    }
   
   fMaxBusPerCh[iChprev-1] = maxBusPatch-iChprev*100+1;

  in.close();

}
