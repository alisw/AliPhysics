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
// MUON Raw Data generator in ALICE-MUON
//
// This class version 1 (further details could be found in Alice-note coming soon right in our direction)
// Generates raw data for MUON tracker and finally for trigger
// * a simple mapping is used (see below)
// * the bus patch id is calculated with an absolute number 0 - 999
// * one DDL per 1/2 chamber is created for both cathode.
// For trigger there is no mapping (mapping could be found in AliMUONTriggerCircuit)
// don't need for digits2raw but needed in raw2Reco
// the position are given per local card
////////////////////////////////////

#include <TClonesArray.h>
#include "AliMUONRawData.h"
#include "AliMUONDigit.h"

#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONConstants.h"
#include "AliMUONData.h"
#include "AliLoader.h"
#include "AliBitPacking.h" 

#include "AliMUONSubEventTrigger.h"
#include "AliMUONDDLTracker.h"
#include "AliMUONDDLTrigger.h"

#include "AliMUONLocalTrigger.h"
#include "AliMUONGlobalTrigger.h"

#include "AliMUONGeometrySegmentation.h"
#include "AliMUONSegmentManuIndex.h"
#include "AliLog.h"
#include "AliRun.h"

const Int_t AliMUONRawData::fgkDefaultPrintLevel = 0;

ClassImp(AliMUONRawData) // Class implementation in ROOT context

//__________________________________________________________________________
AliMUONRawData::AliMUONRawData(AliLoader* loader)
  : TObject(),
    fDebug(0)
{
  // Standard Constructor
 
  fPrintLevel = fgkDefaultPrintLevel;

  // initialize loader's
  fLoader = loader;

  // initialize container
  fMUONData  = new AliMUONData(fLoader,"MUON","MUON");

  // initialize array
  fSubEventArray[0] = new TClonesArray("AliMUONSubEventTracker",1000);
  fSubEventArray[1] = new TClonesArray("AliMUONSubEventTracker",1000);

  // initialize array not used for the moment
  fSubEventTrigArray[0] = new TClonesArray("AliMUONSubEventTrigger",1000);
  fSubEventTrigArray[1] = new TClonesArray("AliMUONSubEventTrigger",1000);

  // ddl pointer
  fDDLTracker = new AliMUONDDLTracker();
  fDDLTrigger = new  AliMUONDDLTrigger();
}

//__________________________________________________________________________
AliMUONRawData::AliMUONRawData()
  : TObject(),
    fMUONData(0),
    fPrintLevel(fgkDefaultPrintLevel),
    fDebug(0),
    fLoader(0)
{
  // Default Constructor
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
  if (fSubEventArray[0])
    fSubEventArray[0]->Delete(); //using delete cos allocating memory in copy ctor.
  if (fSubEventArray[1])
    fSubEventArray[1]->Delete();

  if (fSubEventTrigArray[0])
    fSubEventTrigArray[0]->Delete();
  if (fSubEventTrigArray[1])
    fSubEventTrigArray[1]->Delete();

  if (fDDLTracker)
    delete fDDLTracker;
  if (fDDLTrigger)
    delete fDDLTrigger;

  return;
}
//____________________________________________________________________
Int_t AliMUONRawData::WriteRawData()
{
 // convert digits of the current event to raw data

  Int_t idDDL;
  Char_t name[20];

  fLoader->LoadDigits("READ");

  fMUONData->SetTreeAddress("D,GLT");


  // tracking chambers

  for (Int_t ich = 0; ich < AliMUONConstants::NTrackingCh(); ich++) {
 
    // open files
    idDDL = ich * 2  + 0x900;
    sprintf(name, "MUON_%d.ddl",idDDL);
    fFile1 = fopen(name,"w");

    idDDL = (ich * 2) + 1 + 0x900;
    sprintf(name, "MUON_%d.ddl",idDDL);
    fFile2 = fopen(name,"w");

    WriteTrackerDDL(ich);
  
    // reset and close
    fclose(fFile1);
    fclose(fFile2);
    fMUONData->ResetDigits();
  }
 
  // trigger chambers
 
  // open files
  idDDL = 0xA00;
  sprintf(name, "MUTR_%d.ddl",idDDL);
  fFile1 = fopen(name,"w");

  idDDL = 0xA00 + 1;
  sprintf(name, "MUTR_%d.ddl",idDDL);
  fFile2 = fopen(name,"w");

  WriteTriggerDDL();
  
  // reset and close
  fclose(fFile1);
  fclose(fFile2);
  fMUONData->ResetTrigger();
  
  fLoader->UnloadDigits();

  return kTRUE;
}
//____________________________________________________________________
Int_t AliMUONRawData::WriteTrackerDDL(Int_t iCh)
{
  // resets
  TClonesArray* muonDigits = 0;
  fSubEventArray[0]->Clear();
  fSubEventArray[1]->Clear();

  //
  TArrayI nbInBus[2];

  nbInBus[0].Set(5000);
  nbInBus[1].Set(5000);

  nbInBus[0].Reset();
  nbInBus[1].Reset();

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

  AliMUON *pMUON;
  AliMUONChamber* iChamber = 0x0;
  AliMUONGeometrySegmentation* segmentation2[2];

  segmentation2[0]=iChamber->SegmentationModel2(1); // cathode 0
  segmentation2[1]=iChamber->SegmentationModel2(2); // cathode 1

   pMUON = (AliMUON*) gAlice->GetModule("MUON");
   iChamber =  &(pMUON->Chamber(iCh));


   AliDebug(1, Form("WriteDDL chamber %d\n", iCh+1));

   for (Int_t iCath = 0; iCath < 2; iCath++) {

    fMUONData->ResetDigits();
    fMUONData->GetCathode(iCath);
    muonDigits = fMUONData->Digits(iCh);

    nDigits = muonDigits->GetEntriesFast();
    if (fPrintLevel == 2)
      printf("ndigits = %d\n",nDigits);

    // open DDL file, on per 1/2 chamber
 
    for (Int_t idig = 0; idig < nDigits; idig++) {

      digit = (AliMUONDigit*) muonDigits->UncheckedAt(idig);

      padX = digit->PadX();
      padY = digit->PadY();
      charge = digit->Signal();
      charge &= 0xFFF;
      cathode = digit->Cathode();
      detElemId = digit->DetElemId();

      // mapping
      if (detElemId == 0) {
	AliWarning("\ndetElemId = 0, old segmentation !\n");
	GetDummyMapping(iCh, iCath, digit, busPatchId, manuId, channelId);
      } else {
      // mapping (not far from real one)
	AliMUONSegmentManuIndex* connect = segmentation2[iCath]->GetMpConnection(detElemId, padX, padY);
	if (connect != 0x0) {
	  busPatchId = connect->GetBusPatchId(); 
	  manuId     = connect->GetManuId();
	  channelId  = connect->GetManuChannelId();
	  AliDebug(3,Form("busPatchId %d, manuId: %d, channelId: %d\n", busPatchId, manuId, channelId));
	} else {
	  busPatchId = 0; 
	  manuId     = 0;
	  channelId  = 0;
	}

      }
      //packing word
      AliBitPacking::PackWord((UInt_t)parity,word,29,31);
      AliBitPacking::PackWord((UInt_t)manuId,word,18,28);
      AliBitPacking::PackWord((UInt_t)channelId,word,12,17);
      AliBitPacking::PackWord((UInt_t)charge,word,0,11);

      // set sub Event
      subEvent = new AliMUONSubEventTracker();
      subEvent->AddData(word);
      subEvent->SetBusPatchId(busPatchId);
      if (digit->PadX() > 0) {
	nbInBus[0][busPatchId]++;
	AddData1(subEvent);
      } else {
	nbInBus[1][busPatchId]++;
	AddData2(subEvent);
      }
      delete subEvent;
    }
  }
  fSubEventArray[0]->Sort();
  fSubEventArray[1]->Sort();

  // gather datas from same bus patch
   for (Int_t iDDL = 0; iDDL < 2; iDDL++) {
    nEntries = fSubEventArray[iDDL]->GetEntriesFast();

    for (Int_t i = 0; i < nEntries; i++) {
      AliMUONSubEventTracker* temp = (AliMUONSubEventTracker*)fSubEventArray[iDDL]->At(i);
      busPatchId = temp->GetBusPatchId();

      // add bus patch header, length and total length managed by subevent class
      temp->SetTriggerWord(0xdeadbeef);
      for (Int_t j = 0; j < nbInBus[iDDL][busPatchId]-1; j++) {
	AliMUONSubEventTracker* temp1 =  (AliMUONSubEventTracker*)fSubEventArray[iDDL]->At(++i);
	temp->AddData(temp1->GetData(0));
	fSubEventArray[iDDL]->RemoveAt(i) ;
      }
    }
    fSubEventArray[iDDL]->Compress();

    if (fPrintLevel == 3) {
      nEntries = fSubEventArray[iDDL]->GetEntriesFast();
      for (Int_t i = 0; i < nEntries; i++) {
	AliMUONSubEventTracker* temp =  (AliMUONSubEventTracker*)fSubEventArray[iDDL]->At(i);
	printf("busPatchid back %d\n",temp->GetBusPatchId());
	for (Int_t j = 0; j < temp->GetLength(); j++) {
	  printf("manuId back %d, ",temp->GetManuId(j));
	  printf("channelId back %d, ",temp->GetChannelId(j));
	  printf("charge back %d\n",temp->GetCharge(j));
	}
      }
      printf("\n");
    }
  
  }
  
  Int_t iBusPatch;
  Int_t iEntries;
  Int_t length;

  for (Int_t iDDL = 0; iDDL < 2; iDDL++) {
 

    // filling buffer
    nEntries = fSubEventArray[iDDL]->GetEntriesFast();
    buffer = new Int_t [(2048+24)*50]; // 24 words in average for one buspatch and 2048 manu info at most

    indexBlk = 0;
    indexDsp = 0;
    index = 0;
    iBusPatch = 0;
    iEntries = 0;

    for (Int_t iBlock = 0; iBlock < 2; iBlock++) {

      // block header
      //    fDDLTracker->SetTotalBlkLength(0xFFFFFFFF);
      length = fDDLTracker->GetBlkHeaderLength();
      memcpy(&buffer[index],fDDLTracker->GetBlkHeader(),length*4);
      indexBlk = index;
      index += length; 

      for (Int_t iDsp = 0; iDsp < 5; iDsp++) {

	// DSP header
	//	fDDLTracker->SetTotalDspLength(0xEEEEEEEE);
	length = fDDLTracker->GetDspHeaderLength();
	memcpy(&buffer[index],fDDLTracker->GetDspHeader(),length*4);
	indexDsp = index;
	index += length; 

	for (Int_t i = 0; i < 5; i++) {

	  iBusPatch = i + iBlock*25 + iDsp*5 + 50*(2*iCh + iDDL);

	  AliMUONSubEventTracker* temp = (AliMUONSubEventTracker*)fSubEventArray[iDDL]->At(iEntries);
	  if (nEntries > 0) 
	    busPatchId = temp->GetBusPatchId();
	   else
	    busPatchId = -1;

	  if (busPatchId == iBusPatch) {
	    // add bus patch structure
	    length = temp->GetHeaderLength();
	    memcpy(&buffer[index],temp->GetAddress(),length*4);
	    index += length;
	    for (Int_t j = 0; j < temp->GetLength(); j++) 
	      buffer[index++] =  temp->GetData(j);
	    if (iEntries < nEntries-1)
	      iEntries++;
	  } else {
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
    if (iDDL == 0) {
      // write DDL 1
      header.fSize = (index + headerSize) * 4;// total length in bytes
      fwrite((char*)(&header),headerSize*4,1,fFile1);
      fwrite(buffer,sizeof(int),index,fFile1);
    } 
    if (iDDL == 1) {
      // write DDL 2
      header.fSize = (index + headerSize) * 4;// total length in bytes
      fwrite((char*)(&header),headerSize*4,1,fFile2);
      fwrite(buffer,sizeof(int),index,fFile2);
    }
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
  Char_t locDec, trigY, posY, devX, posX,regOut;
  Int_t version = 1; // software version
  Int_t eventType = 1; // trigger type: 1 for physics ?
  Int_t serialNb = 0xF; // serial nb of card: all bits on for the moment
  Int_t globalFlag = 1; // set to 2 if global info present in DDL else set to 1

  Int_t nEntries = (Int_t) (localTrigger->GetEntries());// 234 local cards
  // stored the local card id that's fired
  for (Int_t i = 0; i <  nEntries; i++) {
    locTrg = (AliMUONLocalTrigger*)localTrigger->At(i);
    isFired[locTrg->LoCircuit()] = 1;
  }

  if (!nEntries)
    AliError("No Trigger information available");

  buffer = new Int_t [672]; // [16(local)*5 words + 3 words]*8(reg) + 8 words = 672

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

      for (Int_t iLoc = 0; iLoc < 16; iLoc++) {

	iLocCard = iLoc + iReg*16 + iDDL*128;

	if (isFired[iLocCard]) {
	  locTrg = (AliMUONLocalTrigger*)localTrigger->At(iEntries);
	  locCard = locTrg->LoCircuit();
	  locDec = locTrg->GetLoDecision();
	  trigY = 0;
	  posY = locTrg->LoStripY();
	  posX = locTrg->LoStripX();
	  devX = locTrg->LoDev();
	  if (fPrintLevel == 4) 
	    printf("loctrg %d, posX %d, posY %d, devX %d\n", 
		   locTrg-> LoCircuit(),locTrg->LoStripX(),locTrg->LoStripY(),locTrg->LoDev());
	} else { //no trigger (see PRR chpt 3.4)
	  locCard = -1;
	  locDec = 0;
	  trigY = 1;
	  posY = 15;
	  posX = 0;
	  devX = 0x8000;
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

    
    if (iDDL == 0) {
      // write DDL 1
      header.fSize = (index + headerSize) * 4;// total length in bytes
      fwrite((char*)(&header),headerSize*4,1,fFile1);
      fwrite(buffer,sizeof(int),index,fFile1);
    } 
    if (iDDL == 1) {
      // write DDL 2
      header.fSize = (index + headerSize) * 4;// total length in bytes
      fwrite((char*)(&header),headerSize*4,1,fFile2);
      fwrite(buffer,sizeof(int),index,fFile2);
    }
  }
  delete[] buffer;

  return kTRUE;
}
//____________________________________________________________________
void AliMUONRawData::GetDummyMapping(Int_t iCh, Int_t iCath, const AliMUONDigit* digit,
				     Int_t &busPatchId, UShort_t &manuId, UChar_t &channelId)
{
// Dummy mapping for tracker

  Int_t offsetX = 0; // offet row
  Int_t offsetY = 0; // offset columns
  Int_t offsetCath = 0; //offset from one cathod to the other
  Int_t maxChannel = 0; // maximum nb of channel in 1/2 chamber
  Int_t id;
      switch (iCh+1) {
      case 1:
      case 2:
      case 3:
      case 4:
	offsetX = 512;
	offsetY = 256;
	offsetCath = 65536;
	maxChannel = (offsetY * offsetX + 2* offsetY + offsetCath);
	break;
      case 5:
      case 6:
      case 7:
      case 8:
      case 9:
      case 10:
	offsetX = 1024;
	offsetY = 0;
	offsetCath = 65536;
	maxChannel = (256 * offsetX + offsetX + offsetCath);
	break;
      }
      // dummy mapping
      // manu Id directly from a matrix 8*8, same segmentation for B and NB
      // 50 buspatches for 1/2 chamber

      id =  (TMath::Abs(digit->PadX()) * offsetX + digit->PadY() + offsetY +
	     offsetCath * iCath);
      Int_t chPerBus = maxChannel/50;
      busPatchId = id/chPerBus; // start at zero 
      if (digit->PadX() > 0)
	busPatchId += 50*iCh*2;
      else 
	busPatchId += 50*(2*iCh+1);
      // 64 manu cards for one buspatch
      manuId = (id % chPerBus)/64; //start at zero 
      manuId &= 0x7FF; // 11 bits 

      // channel id
      channelId = (id % chPerBus) % 64; //start at zero 
      channelId &= 0x3F; // 6 bits

     
      AliDebug(2,Form("id: %d, busPatchId %d, manuId: %d, channelId: %d, maxchannel: %d, chPerBus %d\n",
		      id, busPatchId, manuId, channelId, maxChannel, chPerBus));

      AliDebug(2,Form("id: %d, busPatchId %d, manuId: %d, channelId: %d, padx: %d pady %d, charge %d\n",
		      id, busPatchId, manuId, channelId, digit->PadX(), digit->PadY(), digit->Signal()));

}

//____________________________________________________________________
Int_t AliMUONRawData::GetGlobalTriggerPattern(const AliMUONGlobalTrigger* gloTrg)
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
