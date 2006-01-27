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
// MUON Raw Data generaton in ALICE-MUON
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
// Use memcpy instead of assignment elt by elt
// Introducing variable DSP numbers, real manu numbers per buspatch for st12
// Implemented scaler event for Trigger
// Ch. Finck , Jan. 06
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
#include "AliMUONRawWriter.h"
#include "AliMUONDigit.h"

#include "AliMUONConstants.h"
#include "AliMUONData.h"

#include "AliMUONSubEventTrigger.h"
#include "AliMUONScalerEventTrigger.h"
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

ClassImp(AliMUONRawWriter) // Class implementation in ROOT context

Int_t AliMUONRawWriter::fgManuPerBusSwp1B[12]  = {1, 27, 53, 79, 105, 131, 157, 183, 201, 214, 224, 232};
Int_t AliMUONRawWriter::fgManuPerBusSwp1NB[12] = {1, 27, 53, 79, 105, 131, 157, 183, 201, 214, 225, 233};

Int_t AliMUONRawWriter::fgManuPerBusSwp2B[12]  = {1, 27, 53, 79, 105, 131, 157, 183, 201, 214, 226, 246};
Int_t AliMUONRawWriter::fgManuPerBusSwp2NB[12] = {1, 27, 53, 79, 105, 131, 157, 183, 201, 214, 227, 245};


//__________________________________________________________________________
AliMUONRawWriter::AliMUONRawWriter(AliLoader* loader,  AliMUONData* data)
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

  // initialize array
  fSubEventArray = new TClonesArray("AliMUONSubEventTracker",1000);


  // ddl pointer
  fDDLTracker = new AliMUONDDLTracker();
  fDDLTrigger = new AliMUONDDLTrigger();

  fBusPatchManager = new AliMpBusPatch();
  fBusPatchManager->ReadBusPatchFile();
}

//__________________________________________________________________________
AliMUONRawWriter::AliMUONRawWriter()
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
  fFile[0] = fFile[1] = 0x0;
  
}

//_______________________________________________________________________
AliMUONRawWriter::AliMUONRawWriter (const AliMUONRawWriter& rhs)
  : TObject(rhs)
{
// Protected copy constructor

  AliFatal("Not implemented.");
}

//_______________________________________________________________________
AliMUONRawWriter & 
AliMUONRawWriter::operator=(const AliMUONRawWriter& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}

//__________________________________________________________________________
AliMUONRawWriter::~AliMUONRawWriter(void)
{
  if (fSegFactory) 
    fSegFactory->DeleteSegmentations();
  delete fSegFactory;  

  if (fSubEventArray)
    fSubEventArray->Delete(); //using delete cos allocating memory in copy ctor.

  if (fDDLTracker)
    delete fDDLTracker;
  if (fDDLTrigger)
    delete fDDLTrigger;

  fBusPatchManager->Delete();

  return;
}
//____________________________________________________________________
Int_t AliMUONRawWriter::Digits2Raw()
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
Int_t AliMUONRawWriter::WriteTrackerDDL(Int_t iCh)
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
  fBusPatchManager->GetDspInfo(iCh, iDspMax, iBusPerDSP);

  TArrayI* vec = fBusPatchManager->GetBusfromDE((iCh+1)*100);

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
	  if ((fBusPatchManager->GetDDLfromBus(iBusPatch) % 2) == 1) // comparing to DDL file
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
	    memcpy(&buffer[index],temp->GetBusPatchHeader(),length*4);
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
Int_t AliMUONRawWriter::GetInvMapping(const AliMUONDigit* digit,
				     Int_t &busPatchId, UShort_t &manuId, UChar_t &channelId)
{

  // Inverse mapping for tracker

  Int_t* ptr = 0;

  // information from digits
  Int_t iCath = digit->Cathode();
  Int_t idDE  = digit->DetElemId();
  Int_t padX  = digit->PadX();
  Int_t padY  = digit->PadY();
  Int_t iCh   = idDE/100;

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

   if (idDE < 500) {
     if (iCh == 1 || iCh == 2)
       if (plane == kBendingPlane)
	 ptr = &fgManuPerBusSwp1B[0];
       else 
	 ptr = &fgManuPerBusSwp1NB[0];
     else 
       if (plane == kBendingPlane)
	 ptr = &fgManuPerBusSwp2B[0];
       else 
	 ptr = &fgManuPerBusSwp2NB[0];
   }

  //AliMpVSegmentation* seg = AliMUONSegmentationManager::Segmentation(idDE, plane);
  AliMpVSegmentation* seg = fSegFactory->CreateMpSegmentation(idDE, iCath);
  AliMpPad pad = seg->PadByIndices(AliMpIntPair(padX,padY),kTRUE);

  if (!pad.IsValid()) {
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
  TArrayI* vec = fBusPatchManager->GetBusfromDE(idDE);
  Int_t pos = 0;

  if (idDE < 500) { // station 1 & 2
    for (Int_t i = 0; i < 12; i++)
      if (manuId >= *(ptr + pos++))
	  break;
 //    while(*(ptr + pos) <= manuId)
//       pos++;
//     pos--;
  } else {
    // offset of 100 in manuId for following bus patch
    pos = manuId/100;
  }

  if (pos >(Int_t) vec->GetSize())
    AliWarning(Form("pos greater %d than size %d manuId %d idDE %d \n", 
		    pos, (Int_t)vec->GetSize(), manuId, idDE));
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
Int_t AliMUONRawWriter::WriteTriggerDDL()
{

 // DDL event one per half chamber
  AliMUONSubEventTrigger*    subEvent    = 0x0;
  AliMUONScalerEventTrigger* scalerEvent = 0x0;


  // stored local id number 
  TArrayI isFired(256);
  isFired.Reset();


 // DDL header
  AliRawDataHeader header = fDDLTrigger->GetHeader();
  Int_t headerSize = fDDLTrigger->GetHeaderSize();

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

  if(fScalerEvent)
    eventType = 2; //set to generate scaler events

  Int_t nEntries = (Int_t) (localTrigger->GetEntries());// 234 local cards
  // stored the local card id that's fired
  for (Int_t i = 0; i <  nEntries; i++) {
    locTrg = (AliMUONLocalTrigger*)localTrigger->At(i);
    isFired[locTrg->LoCircuit()] = 1; // storing local boards with informations
  }

  if (!nEntries)
    AliError("No Trigger information available");

  if(fScalerEvent)
    // [16(local)*50 words + 14 words]*8(reg) + 6 + 10 + 6 words scaler event 6534 words
    buffer = new Int_t [6534];
  else
    // [16(local)*5 words + 3 words]*8(reg) + 8 words = 672 
    buffer = new Int_t [672];

  if(fScalerEvent) {
    scalerEvent = new  AliMUONScalerEventTrigger();
    scalerEvent->SetNumbers(); // set some numbers for scalers
  }

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
    buffer[index++]= word;

    if (iDDL == 0)
      fDDLTrigger->SetGlobalOutput(gloTrigPat);// no global input for the moment....
    else 
      fDDLTrigger->SetGlobalOutput(0);

    if (fScalerEvent) {
      // 6 DARC scaler words
      memcpy(&buffer[index], scalerEvent->GetDarcScalers(),scalerEvent->GetDarcScalerLength()*4);
      index += scalerEvent->GetDarcScalerLength();
    }

    // 4 words of global board input + Global board output
    memcpy(&buffer[index], fDDLTrigger->GetGlobalInput(), (fDDLTrigger->GetHeaderLength()-1)*4); 
    index += fDDLTrigger->GetHeaderLength() - 1; // kind tricky cos scaler info in-between Darc header

    if (fScalerEvent) {
      // 10 Global scaler words
      memcpy(scalerEvent->GetGlobalScalers(), &buffer[index], scalerEvent->GetGlobalScalerLength()*4);
      index += scalerEvent->GetGlobalScalerLength();
    }

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
      memcpy(&buffer[index],subEvent->GetRegHeader(),subEvent->GetRegHeaderLength()*4);
      index += subEvent->GetRegHeaderLength();

      // 11 regional scaler word
      if (fScalerEvent) {
	memcpy(&buffer[index], scalerEvent->GetRegScalers(), scalerEvent->GetRegScalerLength()*4);
	index += scalerEvent->GetRegScalerLength();
      }

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
	// 45 regional scaler word
	if (fScalerEvent) {
	  memcpy(&buffer[index], scalerEvent->GetLocalScalers(), scalerEvent->GetLocalScalerLength()*4);
	  index += scalerEvent->GetLocalScalerLength();
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

  if (fScalerEvent)
    delete scalerEvent;

  return kTRUE;
}

//____________________________________________________________________
Int_t AliMUONRawWriter::GetGlobalTriggerPattern(const AliMUONGlobalTrigger* gloTrg) const
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
