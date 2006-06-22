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

/* $Id$ */

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
// Ch. Finck, July 04
// Use memcpy instead of assignment elt by elt
// Introducing variable DSP numbers, real manu numbers per buspatch for st12
// Implemented scaler event for Trigger
// Ch. Finck, Jan. 06
// 
////////////////////////////////////

#include "AliMUONRawWriter.h"

#include "AliBitPacking.h" 
#include "AliRawReader.h"
#include "AliDAQ.h"
#include "AliLog.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"

#include "AliMUONDarcHeader.h"
#include "AliMUONRegHeader.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONDspHeader.h"
#include "AliMUONBlockHeader.h"

#include "AliMUONData.h"
#include "AliMUONDigit.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONLocalTrigger.h"

#include "AliMpBusPatch.h"
#include "AliMpDEManager.h"
#include "AliMpPad.h"
#include "AliMpPlaneType.h"
#include "AliMpSegFactory.h"
#include "AliMpStationType.h"
#include "AliMpVSegmentation.h"

#include "TClonesArray.h"

ClassImp(AliMUONRawWriter) // Class implementation in ROOT context

Int_t AliMUONRawWriter::fgManuPerBusSwp1B[12]  = {1, 27, 53, 79, 105, 131, 157, 183, 201, 214, 224, 232};
Int_t AliMUONRawWriter::fgManuPerBusSwp1NB[12] = {1, 27, 53, 79, 105, 131, 157, 183, 201, 214, 225, 233};

Int_t AliMUONRawWriter::fgManuPerBusSwp2B[12]  = {1, 27, 53, 79, 105, 131, 157, 183, 201, 214, 226, 246};
Int_t AliMUONRawWriter::fgManuPerBusSwp2NB[12] = {1, 27, 53, 79, 105, 131, 157, 183, 201, 214, 227, 245};


//__________________________________________________________________________
AliMUONRawWriter::AliMUONRawWriter(AliMUONData* data)
: TObject(),
  fScalerEvent(kFALSE)
{
  //
  // Standard Constructor
  //
  AliDebug(1,"Standard ctor");
      
  // initialize container
  fMUONData  = data;

  // initialize array
  fBusArray = new TClonesArray("AliMUONBusStruct",1000);
  fBusArray->SetOwner(kTRUE);

  // ddl tracker pointers
  fBlockHeader     = new AliMUONBlockHeader();
  fDspHeader       = new AliMUONDspHeader();
  fBusStruct       = new AliMUONBusStruct();

  // setting data key to default value (only for writting)
  fBlockHeader->SetDataKey(fBlockHeader->GetDefaultDataKey());
  fDspHeader->SetDataKey(fDspHeader->GetDefaultDataKey());
  fBusStruct->SetDataKey(fBusStruct->GetDefaultDataKey());

  // ddl trigger pointers
  fDarcHeader      = new AliMUONDarcHeader();
  fRegHeader       = new AliMUONRegHeader();
  fLocalStruct     = new AliMUONLocalStruct();

  // bus patch & Seg managers
  fBusPatchManager = new AliMpBusPatch();
  fBusPatchManager->ReadBusPatchFile();

  fSegFactory = new AliMpSegFactory();

  // timers
  fTrackerTimer.Start(kTRUE); fTrackerTimer.Stop();
  fTriggerTimer.Start(kTRUE); fTriggerTimer.Stop();
  fMappingTimer.Start(kTRUE); fMappingTimer.Stop();
  
}

//__________________________________________________________________________
AliMUONRawWriter::AliMUONRawWriter()
  : TObject(),
    fMUONData(0),
    fBlockHeader(0),
    fDspHeader(0),
    fBusStruct(0),
    fDarcHeader(0),
    fRegHeader(0),
    fLocalStruct(0),
    fBusPatchManager(0),
    fScalerEvent(kFALSE),
    fSegFactory(0x0)
{
  //
  // Default Constructor
  //
  AliDebug(1,"Default ctor");   
  fFile[0] = fFile[1] = 0x0;  
  fTrackerTimer.Start(kTRUE); fTrackerTimer.Stop();
  fTriggerTimer.Start(kTRUE); fTriggerTimer.Stop();
  fMappingTimer.Start(kTRUE); fMappingTimer.Stop();
}

//_______________________________________________________________________
AliMUONRawWriter::AliMUONRawWriter (const AliMUONRawWriter& rhs)
  : TObject(rhs)
{
  //
  // Protected copy constructor
  //
  AliFatal("Not implemented.");
}

//_______________________________________________________________________
AliMUONRawWriter & 
AliMUONRawWriter::operator=(const AliMUONRawWriter& rhs)
{
  //
  // Protected assignement operator
  //
  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}

//__________________________________________________________________________
AliMUONRawWriter::~AliMUONRawWriter(void)
{
  //
  // Destructor
  //
  AliDebug(1,"dtor");
  
  delete fBusArray;
  
  delete fBlockHeader;
  delete fDspHeader;
  delete fBusStruct;
  delete fDarcHeader;
  delete fRegHeader;
  delete fLocalStruct;

  delete fBusPatchManager;
  
  delete fSegFactory;
  
  AliInfo(Form("Execution time for MUON tracker : R:%.2fs C:%.2fs",
               fTrackerTimer.RealTime(),fTrackerTimer.CpuTime()));
  AliInfo(Form("   Execution time for MUON tracker (mapping calls part) "
               ": R:%.2fs C:%.2fs",
               fMappingTimer.RealTime(),fMappingTimer.CpuTime()));
  AliInfo(Form("Execution time for MUON trigger : R:%.2fs C:%.2fs",
               fTriggerTimer.RealTime(),fTriggerTimer.CpuTime()));
}

//____________________________________________________________________
Int_t AliMUONRawWriter::Digits2Raw()
{
  //
  // convert digits of the current event to raw data
  //
  Int_t idDDL;
  Char_t name[20];

  fMUONData->GetLoader()->LoadDigits("READ");

  fMUONData->SetTreeAddress("D,GLT");

  fMUONData->ResetDigits();
  fMUONData->ResetTrigger();
  
  // This will get both tracker and trigger digits.
  fMUONData->GetDigits();
  
  // tracking chambers

  for (Int_t ich = 0; ich < AliMUONConstants::NTrackingCh(); ich++) 
  {
    // open files
    idDDL = ich * 2  + AliDAQ::DdlIDOffset("MUONTRK");
    sprintf(name, "MUON_%d.ddl",idDDL);
    fFile[0] = fopen(name,"w");

    idDDL = (ich * 2) + 1 + AliDAQ::DdlIDOffset("MUONTRK");
    sprintf(name, "MUON_%d.ddl",idDDL);
    fFile[1] = fopen(name,"w");
    
    WriteTrackerDDL(ich);
  
    // reset and close
    fclose(fFile[0]);
    fclose(fFile[1]);
  }
 
  // trigger chambers
 
  // open files
  idDDL = AliDAQ::DdlIDOffset("MUONTRG");
  sprintf(name, "MUTR_%d.ddl",idDDL);
  fFile[0] = fopen(name,"w");

  idDDL = AliDAQ::DdlIDOffset("MUONTRG") + 1;
  sprintf(name, "MUTR_%d.ddl",idDDL);
  fFile[1] = fopen(name,"w");

  WriteTriggerDDL();
  
  // reset and close
  fclose(fFile[0]);
  fclose(fFile[1]);

  fMUONData->ResetDigits();
  fMUONData->ResetTrigger();  
  fMUONData->GetLoader()->UnloadDigits();

  return kTRUE;
}

//____________________________________________________________________
Int_t AliMUONRawWriter::WriteTrackerDDL(Int_t iCh)
{
  // writing DDL for tracker
  // used inverse mapping
  //
  fTrackerTimer.Start(kFALSE);
  

  static const Int_t kMAXADC = (1<<12)-1; // We code the charge on a 12 bits ADC.

  // resets
  TClonesArray* muonDigits = 0;

  fBusArray->Delete();


  //
  TArrayI nbInBus;

  nbInBus.Set(5000);

  nbInBus.Reset();

  // DDL header
  Int_t headerSize = sizeof(fHeader)/4;

  // DDL event one per half chamber

  // raw data
  Char_t parity = 0x4;
  UShort_t manuId = 0;
  UChar_t channelId = 0;
  UShort_t charge = 0;
  Int_t busPatchId = 0;
  UInt_t word;


  // Dsp length
  Int_t totalDspLength;
  Int_t dspLength;

  // block length
  Int_t totalBlkLength;
  Int_t blkLength; 
  
  // total DDL length
  Int_t totalDDLLength;

  // indexes
  Int_t index;
  Int_t indexDsp;
  Int_t indexBlk;

  // digits
  Int_t nEntries = 0;
  Int_t* buffer = 0;
  Int_t padX;
  Int_t padY;
  Int_t cathode = 0;
  Int_t detElemId;
  Int_t nDigits;

  const AliMUONDigit* digit;

  AliDebug(3, Form("WriteDDL chamber %d\n", iCh+1));

  muonDigits = fMUONData->Digits(iCh);

  nDigits = muonDigits->GetEntriesFast();
  AliDebug(3,Form("ndigits = %d\n",nDigits));
 
  // loop over digit
  for (Int_t idig = 0; idig < nDigits; idig++) {

    digit = (AliMUONDigit*) muonDigits->UncheckedAt(idig);

    padX = digit->PadX();
    padY = digit->PadY();
    charge = digit->ADC();
    if ( charge > kMAXADC )
    {
      // This is most probably an error in the digitizer (which should insure
      // the adc is below kMAXADC), so make it a (non-fatal) error indeed.
      AliError(Form("adc value %d above %x. Setting to %x",
                      charge,kMAXADC,kMAXADC));
      charge = kMAXADC;
    }
    cathode = digit->Cathode();
    detElemId = digit->DetElemId();

    // inverse mapping
    busPatchId = GetBusPatch(*digit);
    if (busPatchId<0) continue;

    if ( digit->ManuId() > 0x7FF || digit->ManuId() < 0 ||
         digit->ManuChannel() > 0x3F || digit->ManuChannel() < 0 )
    {
      StdoutToAliError(digit->Print(););
      AliFatal("ManuId,ManuChannel are invalid for the digit above.");
    }
    
    manuId = ( digit->ManuId() & 0x7FF ); // 11 bits
    channelId = ( digit->ManuChannel() & 0x3F ); // 6 bits

    AliDebug(3,Form("input  IdDE %d busPatchId %d PadX %d PadY %d iCath %d \n", 
		    detElemId, busPatchId, padX, padY, cathode));

    AliDebug(3,Form("busPatchId %d, manuId %d channelId %d\n", busPatchId, manuId, channelId ));

    //packing word
    word = 0;
    AliBitPacking::PackWord((UInt_t)manuId,word,18,28);
    AliBitPacking::PackWord((UInt_t)channelId,word,12,17);
    AliBitPacking::PackWord((UInt_t)charge,word,0,11);

    // parity word
    parity = word & 0x1;
    for (Int_t i = 1; i <= 30; i++) 
      parity ^=  ((word >> i) & 0x1);
    AliBitPacking::PackWord((UInt_t)parity,word,31,31);

    // set sub Event
    fBusStruct->SetLength(0);
    fBusStruct->AddData(word);
    fBusStruct->SetBusPatchId(busPatchId);
       
    // storing the number of identical buspatches
    nbInBus[busPatchId]++;
    AddData(*fBusStruct);
   
  }

  // sorting by buspatch
  fBusArray->Sort();

  // gather datas from same bus patch
  nEntries = fBusArray->GetEntriesFast();

  for (Int_t i = 0; i < nEntries; i++) {
    AliMUONBusStruct* temp = (AliMUONBusStruct*)fBusArray->At(i);
    busPatchId = temp->GetBusPatchId();

    // add bus patch header, length and total length managed by subevent class
    for (Int_t j = 0; j < nbInBus[busPatchId]-1; j++) {
      AliMUONBusStruct* temp1 =  (AliMUONBusStruct*)fBusArray->At(++i);
      temp->AddData(temp1->GetData(0));
      fBusArray->RemoveAt(i) ;
    }
  }
  fBusArray->Compress();

  if (AliLog::GetGlobalDebugLevel() == 3) {
    nEntries = fBusArray->GetEntriesFast();
    for (Int_t i = 0; i < nEntries; i++) {
      AliMUONBusStruct* temp =  (AliMUONBusStruct*)fBusArray->At(i);
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

  // nEntries = fBusArray->GetEntriesFast();

  AliMUONBusStruct* busStructPtr = 0x0;

  // open DDL file, on per 1/2 chamber
  for (Int_t iDDL = 0; iDDL < 2; iDDL++) {
    
    totalDDLLength = 0;

    // filling buffer
    buffer = new Int_t [(2048+24)*50]; // 24 words at most for one buspatch and 2048 manu info at most

    indexBlk = 0;
    indexDsp = 0;
    index = 0;

    // two blocks A and B per DDL
    for (Int_t iBlock = 0; iBlock < 2; iBlock++) {

      // block header
      length = fBlockHeader->GetHeaderLength();
      memcpy(&buffer[index],fBlockHeader->GetHeader(),length*4);
      indexBlk = index;
      index += length; 

      // 5 DSP's max per block
      for (Int_t iDsp = 0; iDsp < iDspMax; iDsp++) {

	// DSP header
	length = fDspHeader->GetHeaderLength();
	memcpy(&buffer[index],fDspHeader->GetHeader(),length*4);
	indexDsp = index;
	index += length; 

	// 5 buspatches max per DSP
	for (Int_t i = 0; i < iBusPerDSP[iDsp]; i++) {

	  iBusPatch++;
	  if ((fBusPatchManager->GetDDLfromBus(iBusPatch) % 2) == 0) // comparing to DDL file
	    iFile = 1;
	  else
	    iFile = 0;

	  AliDebug(3,Form("iCh %d iDDL %d iBlock %d iDsp %d busPatchId %d", iCh, iDDL, iBlock, iDsp, iBusPatch));

	  nEntries = fBusArray->GetEntriesFast();
	  busPatchId = -1;

	  for (Int_t iEntries = 0; iEntries < nEntries; iEntries++) { // method "bourrique"...
	    busStructPtr = (AliMUONBusStruct*)fBusArray->At(iEntries);
	    busPatchId = busStructPtr->GetBusPatchId();
	    if (busPatchId == iBusPatch) break;
	    busPatchId = -1;
	    AliDebug(3,Form("busPatchId %d", busStructPtr->GetBusPatchId()));
	  } 
	 
	  // check if buspatchid has digit
	  if (busPatchId != -1) {
	    // add bus patch structure header
	    length = busStructPtr->GetHeaderLength();
	    memcpy(&buffer[index],busStructPtr->GetHeader(),length*4);
	    index += length;

	    // add bus patch data
	    for (Int_t j = 0; j < busStructPtr->GetLength(); j++) {
	      buffer[index++] =  busStructPtr->GetData(j);
	      AliDebug(3,Form("busPatchId %d, manuId %d channelId %d\n", 
			      busStructPtr->GetBusPatchId(), 
			      busStructPtr->GetManuId(j), busStructPtr->GetChannelId(j) ));
	    }
	    //	      fBusArray->RemoveAt(iEntries);
	    //	      fBusArray->Compress();
	  } else {
	    // writting anyhow buspatch structure (empty ones)
	    buffer[index++] = busStructPtr->GetDefaultDataKey(); // fill it also for empty data size
	    buffer[index++] = busStructPtr->GetHeaderLength(); // header length
	    buffer[index++] = 0; // raw data length
	    buffer[index++] = iBusPatch; // bus patch
	  }
	} // bus patch

	// check if totalLength even
	// set padding word in case
	// Add one word 0xBEEFFACE at the end of DSP structure
	totalDspLength  = index - indexDsp;
	if ((totalDspLength % 2) == 1) { 
	  buffer[indexDsp + fDspHeader->GetHeaderLength() - 2] = 1;
	  buffer[index++] = fDspHeader->GetDefaultPaddingWord();
	  totalDspLength++;
	}

	dspLength          = totalDspLength - fDspHeader->GetHeaderLength();

	buffer[indexDsp+1] = totalDspLength; // dsp total length
	buffer[indexDsp+2] = dspLength; // data length  
	   
      } // dsp

      totalBlkLength  = index - indexBlk;
      blkLength       = totalBlkLength - fBlockHeader->GetHeaderLength();
      totalDDLLength += totalBlkLength;

      buffer[indexBlk+1] = totalBlkLength; // total block length
      buffer[indexBlk+2] = blkLength;

    } // block
    
    //writting onto disk
    // write DDL 1 & 2
    fHeader.fSize = (totalDDLLength + headerSize) * 4;// total length in bytes
    fwrite((char*)(&fHeader),headerSize*4,1,fFile[iFile]);
    fwrite(buffer,sizeof(int),index,fFile[iFile]);
   
    delete[] buffer;
  }

  fTrackerTimer.Stop();
  return kTRUE;
}

//____________________________________________________________________
Int_t AliMUONRawWriter::GetBusPatch(const AliMUONDigit& digit)
{
  //
  // Determine the BusPatch this digit belongs to.
  //
  fMappingTimer.Start(kFALSE);
  
  Int_t* ptr = 0;

  // information from digits
  Int_t detElemId  = digit.DetElemId();

  AliMpVSegmentation* seg = 
    fSegFactory->CreateMpSegmentationByElectronics(detElemId, digit.ManuId());
  
  AliMpPlaneType plane = seg->PlaneType();

  AliMpStationType stationType = AliMpDEManager::GetStationType(detElemId);

  if ( stationType == kStation1 || stationType == kStation2 )
  {
    if (plane == kBendingPlane) 
    {
      ptr = &fgManuPerBusSwp1B[0];
    }
    else 
    {
      ptr = &fgManuPerBusSwp1NB[0];
    }
  }
  else
  {
    if (plane == kBendingPlane)
    {
      ptr = &fgManuPerBusSwp2B[0];
    }
    else
    {
      ptr = &fgManuPerBusSwp2NB[0];
    }
  }

  // Getting buspatch id
  TArrayI* vec = fBusPatchManager->GetBusfromDE(detElemId);
  Int_t pos = 0;

  Int_t m = ( digit.ManuId() & 0x3FF ); // remove bit 10
                                //FIXME : how can we remove that condition
  // on the 10-th bit ? All the rest need not any knowledge about it,
  // can't we find a way to get manu<->buspatch transparent to this too ?
  
  if ( stationType == kStation1 || stationType == kStation2 )
  {
    for (Int_t i = 0; i < 12; i++)
    {
      if (m >= *(ptr + pos++)) break;
    }
  }
  else 
  {
    // offset of 100 in manuId for following bus patch
    pos = m/100;
  }

  if (pos >(Int_t) vec->GetSize())
  {
    AliError(Form("pos greater %d than size %d manuId %d detElemId %d \n", 
		    pos, (Int_t)vec->GetSize(), digit.ManuId(), detElemId));
    AliError(Form("Chamber %s Plane %s manuId %d m %d",
                    StationTypeName(stationType).Data(),
                    PlaneTypeName(plane).Data(),
                    digit.ManuId(),
                    m));
    return -1;
  }
  
  Int_t busPatchId = vec->At(pos);

  fMappingTimer.Stop();
  
  return busPatchId;
}

//____________________________________________________________________
Int_t AliMUONRawWriter::WriteTriggerDDL()
{
  //
  // Write trigger DDL
  //
  fTriggerTimer.Start(kFALSE);
  
 // DDL event one per half chamber

  // stored local id number 
  TArrayI isFired(256);
  isFired.Reset();


 // DDL header size
  Int_t headerSize = sizeof(AliRawDataHeader)/4;

  TClonesArray* localTrigger;
  TClonesArray* globalTrigger;
  AliMUONGlobalTrigger* gloTrg;
  AliMUONLocalTrigger* locTrg = 0x0;

  // global trigger for trigger pattern
  globalTrigger = fMUONData->GlobalTrigger(); 
  gloTrg = (AliMUONGlobalTrigger*)globalTrigger->UncheckedAt(0);
  Int_t gloTrigPat = gloTrg->GetGlobalPattern();

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
    AliInfo("No Trigger information available");

  if(fScalerEvent)
    // [16(local)*51 words + 16 words]*8(reg) + 6 + 12 + 6 words scaler event 6672 words
    buffer = new Int_t [6680];
  else
    // [16(local)*6 words + 5 words]*8(reg) + 10 words = 818 
    buffer = new Int_t [818];


  // open DDL file, on per 1/2 chamber
  for (Int_t iDDL = 0; iDDL < 2; iDDL++) {
    
    index = 0; 

    word = 0;
    // set darc status word
    AliBitPacking::PackWord((UInt_t)iDDL+1,word,28,31); //see AliMUONDDLTrigger.h for details
    AliBitPacking::PackWord((UInt_t)serialNb,word,24,27);
    AliBitPacking::PackWord((UInt_t)version,word,16,23);
    AliBitPacking::PackWord((UInt_t)eventType,word,12,15);

    if (iDDL == 0) // suppose global info in DDL one
      globalFlag = 2;
    else 
      globalFlag = 1;

    AliBitPacking::PackWord((UInt_t)globalFlag,word,8,11);
    fDarcHeader->SetWord(word);

    memcpy(&buffer[index], fDarcHeader->GetHeader(), (fDarcHeader->GetDarcHeaderLength())*4); 
    index += fDarcHeader->GetDarcHeaderLength();

    if (iDDL == 0)
     fDarcHeader->SetGlobalOutput(gloTrigPat);// no global input for the moment....
    else 
     fDarcHeader->SetGlobalOutput(0);

    if (fScalerEvent) {
      // 6 DARC scaler words
      memcpy(&buffer[index], fDarcHeader->GetDarcScalers(),fDarcHeader->GetDarcScalerLength()*4);
      index += fDarcHeader->GetDarcScalerLength();
    }
    // end of darc word
    buffer[index++] = fDarcHeader->GetEndOfDarc();

    // 4 words of global board input + Global board output
    memcpy(&buffer[index], fDarcHeader->GetGlobalInput(), (fDarcHeader->GetGlobalHeaderLength())*4); 
    index += fDarcHeader->GetGlobalHeaderLength(); 

    if (fScalerEvent) {
      // 10 Global scaler words
      memcpy(fDarcHeader->GetGlobalScalers(), &buffer[index], fDarcHeader->GetGlobalScalerLength()*4);
      index += fDarcHeader->GetGlobalScalerLength();
    }

    // end of global word
    buffer[index++] = fDarcHeader->GetEndOfGlobal();

    // 8 regional cards per DDL
    for (Int_t iReg = 0; iReg < 8; iReg++) {

      // Regional card header
      word = 0;

      // set darc status word
      fRegHeader->SetDarcWord(word);

      regOut  = 0;
      AliBitPacking::PackWord((UInt_t)serialNb,word,24,28); //see  AliMUONLocalStruct.h for details
      AliBitPacking::PackWord((UInt_t)version,word,16,23);
      AliBitPacking::PackWord((UInt_t)iReg,word,12,15);
      AliBitPacking::PackWord((UInt_t)regOut,word,0,7); // whenever regional output will be implemented

      fRegHeader->SetWord(word);
      memcpy(&buffer[index],fRegHeader->GetHeader(),fRegHeader->GetHeaderLength()*4);
      index += fRegHeader->GetHeaderLength();

      // 11 regional scaler word
      if (fScalerEvent) {
	memcpy(&buffer[index], fRegHeader->GetScalers(), fRegHeader->GetScalerLength()*4);
	index += fRegHeader->GetScalerLength();
      }

      // end of regional word
      buffer[index++] = fRegHeader->GetEndOfReg();

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
	  locCard = -1; // not possible on 4 bits
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
	  memcpy(&buffer[index], fLocalStruct->GetScalers(), fLocalStruct->GetScalerLength()*4);
	  index += fLocalStruct->GetScalerLength();
	}

	// end of local structure words
 	buffer[index++] = fLocalStruct->GetEndOfLocal();

      } // local card 

    } // Regional card
    

    // writting onto disk
    // write DDL 1
    fHeader.fSize = (index + headerSize) * 4;// total length in bytes
    fwrite((char*)(&fHeader),headerSize*4,1,fFile[iDDL]);
    fwrite(buffer,sizeof(int),index,fFile[iDDL]);
  
  }
  delete[] buffer;

  fTriggerTimer.Stop();
  
  return kTRUE;
}
//____________________________________________________________________
void AliMUONRawWriter::SetScalersNumbers()
{
  // set numbers for scaler events for trigger headers
  // since this is provided by the experiment
  // put dummy numbers to check the monitoring

  fDarcHeader->SetScalersNumbers();
  fRegHeader->SetScalersNumbers();
  fLocalStruct->SetScalersNumbers();
 
  fScalerEvent = kTRUE;
}
