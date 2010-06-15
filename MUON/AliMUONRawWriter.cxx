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

//-----------------------------------------------------------------------------
/// \class AliMUONRawWriter
/// MUON Raw Data generaton in ALICE-MUON
/// Raw data structure could be found in Alice-note.
///
/// Implemented non-constant buspatch numbers for tracking
/// with correct DDL id (first guess)
/// (Ch. Finck, dec 2005)
///
/// Digits2Raw:
/// Generates raw data for MUON tracker and finally for trigger
/// Using real mapping (inverse) for tracker
/// For trigger there is no mapping (mapping could be found in AliMUONTriggerCircuit)
/// Ch. Finck, July 04
/// Use memcpy instead of assignment elt by elt
/// Introducing variable DSP numbers, real manu numbers per buspatch for st12
/// Implemented scaler event for Trigger
/// Ch. Finck, Jan. 06
/// Using bus itr in DDL instead of simple incrementation
/// treat correctly the DDL & buspatch for station 3.
/// Using informations from AliMUONTriggerCrateStore for 
/// empty slots and non-notified cards in trigger crates.
/// Ch. Finck, August 06.
/// Using AliMpDDLStore::GetBusPatchId.
///
/// \author Ch. Finck, Feb. 07.
//-----------------------------------------------------------------------------


#include "AliMUONRawWriter.h"

#include "AliMUONBlockHeader.h"
#include "AliMUONBusStruct.h"
#include "AliMUONConstants.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONDspHeader.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONRegionalTrigger.h"
#include "AliMUONRegHeader.h"

#include "AliMUONVTriggerStore.h"
#include "AliCodeTimer.h"

#include "AliMpDDLStore.h"
#include "AliMpDDL.h"
#include "AliMpRegionalTrigger.h"
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"
#include "AliMpDetElement.h"
#include "AliMpDEManager.h"
#include "AliMpExMap.h"
#include "AliMpConstants.h"
#include "AliMpPlaneType.h"
#include "AliMpSegmentation.h"
#include "AliMpStationType.h"
#include "AliMpVSegmentation.h"

#include "AliRawReader.h"
#include "AliRawDataHeaderSim.h"
#include "AliBitPacking.h" 
#include "AliDAQ.h"
#include "AliLog.h"

#include "TObjArray.h"
#include "TStopwatch.h"
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONRawWriter) // Class implementation in ROOT context
/// \endcond

//__________________________________________________________________________
AliMUONRawWriter::AliMUONRawWriter()
  : TObject(),
    fBlockHeader(new AliMUONBlockHeader()),
    fDspHeader(new AliMUONDspHeader()),
    fDarcHeader(new AliMUONDarcHeader()),
    fRegHeader(new AliMUONRegHeader()),
    fLocalStruct(new AliMUONLocalStruct()),
    fDDLStore(AliMpDDLStore::Instance()),
    fScalerEvent(kFALSE),
    fHeader(0x0),
    fBufferSize((((43*AliMpConstants::ManuNofChannels() + 4)*5 + 10)*5 + 8)*2),
    fBuffer(new Int_t [fBufferSize])
{
  /// Standard Constructor

  AliDebug(1,"Standard ctor");

  // setting data key to default value (only for writting)
  fBlockHeader->SetDataKey(fBlockHeader->GetDefaultDataKey());
  fDspHeader->SetDataKey(fDspHeader->GetDefaultDataKey());

}

//__________________________________________________________________________
AliMUONRawWriter::~AliMUONRawWriter(void)
{
  /// Destructor

  AliDebug(1,"dtor");
  
  delete fBlockHeader;
  delete fDspHeader;
  delete fDarcHeader;
  delete fRegHeader;
  delete fLocalStruct;
  delete[] fBuffer;
}

//____________________________________________________________________
void  AliMUONRawWriter::LocalWordPacking(UInt_t& word, UInt_t locId, UInt_t locDec, 
					 UInt_t trigY, UInt_t posY, UInt_t posX, 
					 UInt_t sdevX, UInt_t devX)
{
/// pack local trigger word

    AliBitPacking::PackWord(locId,word,19,22); //card id number in crate
    AliBitPacking::PackWord(locDec,word,15,18);
    AliBitPacking::PackWord(trigY,word,14,14);
    AliBitPacking::PackWord(posY,word,10,13);
    AliBitPacking::PackWord(sdevX,word,9,9);
    AliBitPacking::PackWord(devX,word,5,8);
    AliBitPacking::PackWord(posX,word,0,4);

}

//____________________________________________________________________
Int_t AliMUONRawWriter::Digits2Raw(const AliMUONVDigitStore* digitStore,
                                   const AliMUONVTriggerStore* triggerStore)
{
  /// convert digits of the current event to raw data

  AliCodeTimerAuto("",0)
  
  Int_t idDDL;
  Char_t name[255];

  // tracking chambers
  
  if ( digitStore ) 
  {
    AliCodeTimerAuto("for Tracker",1)

    AliMpExMap busPatchMap;

    Int_t nDDLs = AliDAQ::NumberOfDdls("MUONTRK");
    
    Int_t nofBusPatches(0);
    
    for (Int_t iDDL = 0; iDDL < nDDLs; ++iDDL ) 
    {
      AliMpDDL* ddl = fDDLStore->GetDDL(iDDL);
      nofBusPatches += ddl->GetNofBusPatches();
    }
    
    busPatchMap.SetSize(nofBusPatches);
    
    Digits2BusPatchMap(*digitStore,busPatchMap);

    for (Int_t iDDL = 0; iDDL < nDDLs; ++iDDL ) 
    {
      WriteTrackerDDL(busPatchMap,iDDL);
    }
    AliDebug(1,"Tracker written");
  }
 
  if ( triggerStore )
  {
    AliCodeTimerAuto("for Trigger",1)

    // trigger chambers
    
    AliFstream* file[2];
    
    // open files
    idDDL = 0;// MUTR
    strcpy(name,AliDAQ::DdlFileName("MUONTRG",idDDL));
    file[0] = new AliFstream(name);
    
    idDDL = 1;// MUTR
    strcpy(name,AliDAQ::DdlFileName("MUONTRG",idDDL));
    file[1] = new AliFstream(name);
      
    WriteTriggerDDL(*triggerStore,file);
      
    // reset and close
    delete file[0];
    delete file[1];
      
    AliDebug(1,"Trigger written");
  }
  
  return kTRUE;
}

//______________________________________________________________________________
void 
AliMUONRawWriter::Digits2BusPatchMap(const AliMUONVDigitStore& digitStore,
                                     AliMpExMap& busPatchMap)
{
  /// Create bus patch structures corresponding to digits in the store
  
  AliCodeTimerAuto("",0)
  
  static const Int_t kMAXADC = (1<<12)-1; // We code the charge on a 12 bits ADC.
    
  // DDL event one per half chamber
  
  // raw data
  Char_t parity = 0x4;
  UShort_t manuId = 0;
  UChar_t channelId = 0;
  UShort_t charge = 0;
  Int_t busPatchId = 0;
  Int_t currentBusPatchId = -1;
  UInt_t word;
  
  AliMUONBusStruct* busStruct(0x0);
  
  TIter next(digitStore.CreateTrackerIterator());
  AliMUONVDigit* digit;
  
  while ( ( digit = static_cast<AliMUONVDigit*>(next()) ) )
  {
    charge = digit->ADC();
    if ( charge > kMAXADC )
    {
      // This is most probably an error in the digitizer (which should insure
      // the adc is below kMAXADC), so make it a (non-fatal) error indeed.
      AliError(Form("adc value %d above 0x%x for DE %d . Setting to 0x%x. Digit is:",
                    charge,kMAXADC,digit->DetElemId(),kMAXADC));
      StdoutToAliError(digit->Print());
      charge = kMAXADC;
    }
    
    // inverse mapping
    busPatchId = GetBusPatch(*digit);

    if (busPatchId<0) continue;
    
    if ( digit->ManuId() > 0x7FF ||
         digit->ManuChannel() > 0x3F )
    {
      StdoutToAliError(digit->Print(););
      AliFatal("ManuId,ManuChannel are invalid for the digit above.");
    }
    
    manuId = ( digit->ManuId() & 0x7FF ); // 11 bits
    channelId = ( digit->ManuChannel() & 0x3F ); // 6 bits
    
    //packing word
    word = 0;
    AliBitPacking::PackWord((UInt_t)manuId,word,18,28);
    AliBitPacking::PackWord((UInt_t)channelId,word,12,17);
    AliBitPacking::PackWord((UInt_t)charge,word,0,11);
    
    // parity word
    parity = word & 0x1;
    for (Int_t i = 1; i <= 30; ++i) 
    {
      parity ^=  ((word >> i) & 0x1);
    }
    AliBitPacking::PackWord((UInt_t)parity,word,31,31);

    if ( currentBusPatchId != busPatchId ) 
    {
      busStruct = 
        static_cast<AliMUONBusStruct*>(busPatchMap.GetValue(busPatchId));
      currentBusPatchId = busPatchId;
    }
    
    if (!busStruct)
    {
      busStruct = new AliMUONBusStruct;
      busStruct->SetDataKey(busStruct->GetDefaultDataKey());
      busStruct->SetBusPatchId(busPatchId);
      busStruct->SetLength(0);
      busPatchMap.Add(busPatchId,busStruct);
    }
    
    // set sub Event
    busStruct->AddData(word);
  }
}

//______________________________________________________________________________
void
AliMUONRawWriter::WriteTrackerDDL(AliMpExMap& busPatchMap, Int_t iDDL)
{
  /// Write DDL file for one tracker DDL
  
  // buffer size (max'ed out)
  // (((43 manus max per bus patch *64 channels + 4 bus patch words) * 5 bus patch 
  //   + 10 dsp words)*5 dsps + 8 block words)*2 blocks 
 
  AliCodeTimerAuto("",0)

  if (fHeader == 0x0) {
    AliError("Raw data header must be set");
    return;
  }
  memset(fBuffer,0,fBufferSize*sizeof(Int_t));
  
  AliMpDDL* ddl = fDDLStore->GetDDL(iDDL);
  Int_t iDspMax = ddl->GetMaxDsp();
  Int_t iBusPerDSP[5]; //number of bus patches per DSP
  ddl->GetBusPerDsp(iBusPerDSP);
  Int_t busIter = 0;
  
  Int_t totalDDLLength = 0;
  
  Int_t index = 0;
  
  // two blocks A and B per DDL
  for (Int_t iBlock = 0; iBlock < 2; ++iBlock) 
  {
    // block header
    Int_t length = fBlockHeader->GetHeaderLength();
    memcpy(&fBuffer[index],fBlockHeader->GetHeader(),length*4);
    Int_t indexBlk = index;
    index += length; 
    
    // 5 DSP's max per block
    for (Int_t iDsp = 0; iDsp < iDspMax; ++iDsp) 
    {
      // DSP header
      Int_t dspHeaderLength = fDspHeader->GetHeaderLength();
      memcpy(&fBuffer[index],fDspHeader->GetHeader(),dspHeaderLength*4);
      Int_t indexDsp = index;
      index += dspHeaderLength; 
      
      // 5 buspatches max per DSP
      for (Int_t i = 0; i < iBusPerDSP[iDsp]; ++i) 
      {
        Int_t iBusPatch = ddl->GetBusPatchId(busIter++);
        
        // iteration over bus patch in DDL
        if (iBusPatch == -1) 
        {
          AliWarning(Form("Error in bus itr in DDL %d\n", iDDL));
          continue;
        }
        
        AliMUONBusStruct* busStructPtr = static_cast<AliMUONBusStruct*>(busPatchMap.GetValue(iBusPatch));
        
        // check if buspatchid has digit
        if (busStructPtr) 
        {
          // add bus patch structure header
          Int_t busHeaderLength = busStructPtr->GetHeaderLength();
          memcpy(&fBuffer[index],busStructPtr->GetHeader(),busHeaderLength*4);
          index += busHeaderLength;
          
          // add bus patch data
          Int_t busLength = busStructPtr->GetLength();
          memcpy(&fBuffer[index],busStructPtr->GetData(),busLength*4);
          index += busLength;
        } 
        else 
        {
          // writting anyhow buspatch structure (empty ones)
          fBuffer[index++] = busStructPtr->GetDefaultDataKey(); // fill it also for empty data size
          fBuffer[index++] = busStructPtr->GetHeaderLength(); // header length
          fBuffer[index++] = 0; // raw data length
          fBuffer[index++] = iBusPatch; // bus patch
        }
      } // bus patch
      
      // check if totalLength even
      // set padding word in case
      // Add one word 0xBEEFFACE at the end of DSP structure
      Int_t totalDspLength  = index - indexDsp;
      if ((totalDspLength % 2) == 1) 
      { 
        fBuffer[indexDsp + fDspHeader->GetHeaderLength() - 2] = 1;
        fBuffer[index++] = fDspHeader->GetDefaultPaddingWord();
        totalDspLength++;
      }
      
      Int_t dspLength     = totalDspLength - fDspHeader->GetHeaderLength();
      
      fBuffer[indexDsp+1] = totalDspLength; // dsp total length
      fBuffer[indexDsp+2] = dspLength; // data length  
      
    } // dsp
    
    Int_t totalBlkLength  = index - indexBlk;
    Int_t blkLength       = totalBlkLength - fBlockHeader->GetHeaderLength();
    totalDDLLength       += totalBlkLength;
    
    fBuffer[indexBlk+1] = totalBlkLength; // total block length
    fBuffer[indexBlk+2] = blkLength;
        
  } // block
  
    // add twice the end of CRT structure data key
    // hope it's good placed (ChF)
    fBuffer[index++] = fBlockHeader->GetDdlDataKey();
    fBuffer[index++] = fBlockHeader->GetDdlDataKey();
    totalDDLLength  += 2;
  
  // writting onto disk
  // total length in bytes
  // DDL header

  Int_t headerSize = sizeof(AliRawDataHeader)/4;
  
  fHeader->fSize = (totalDDLLength + headerSize) * 4;
  
  AliFstream* file = new AliFstream(AliDAQ::DdlFileName("MUONTRK",iDDL));
  
  file->WriteBuffer((char*)fHeader,headerSize*4);
  file->WriteBuffer((char*)fBuffer,sizeof(int)*index);
  delete file;
}

//______________________________________________________________________________
Int_t AliMUONRawWriter::GetBusPatch(const AliMUONVDigit& digit) const
{
  /// Determine the BusPatch this digit belongs to.

    return fDDLStore->GetBusPatchId(digit.DetElemId(),digit.ManuId());
}

//______________________________________________________________________________
Int_t AliMUONRawWriter::WriteTriggerDDL(const AliMUONVTriggerStore& triggerStore, AliFstream* file[2])
{
  /// Write trigger DDL
  
  AliCodeTimerAuto("",0)

  if (fHeader == 0x0) {
    AliError("Raw data header must be set");
    return 0;
  }

 // DDL event one per half chamber

 // DDL header size
  Int_t headerSize = sizeof(AliRawDataHeader)/4;

  // global trigger for trigger pattern
  AliMUONGlobalTrigger* gloTrg = triggerStore.Global();
  if (!gloTrg) 
  {
    return 0;
  }
  
  Int_t gloTrigResp = gloTrg->GetGlobalResponse();
  UInt_t *gloTrigInput = gloTrg->GetGlobalInput();

  UInt_t word;
  Int_t* buffer = 0;
  Int_t index;
  Int_t locCard;
  UChar_t locDec, trigY, posY, posX, regOut;
  UInt_t regInpLpt;
  UInt_t regInpHpt;

  UInt_t devX;
  UChar_t sdevX;
  UInt_t version = 1; // software version
  UInt_t eventPhys = 1; // trigger type: 1 for physics, 0 for software
  UInt_t serialNb = 0xF; // serial nb of card: all bits on for the moment
  Int_t globalFlag = 0; // set to 1 if global info present in DDL else set to 0

  // size of headers
  static const Int_t kDarcHeaderLength   = fDarcHeader->GetDarcHeaderLength();
  static const Int_t kGlobalHeaderLength = fDarcHeader->GetGlobalHeaderLength();
  static const Int_t kDarcScalerLength   = fDarcHeader->GetDarcScalerLength();
  static const Int_t kGlobalScalerLength = fDarcHeader->GetGlobalScalerLength();
  static const Int_t kRegHeaderLength    = fRegHeader->GetHeaderLength();
  static const Int_t kRegScalerLength    = fRegHeader->GetScalerLength();
  static const Int_t kLocHeaderLength    = fLocalStruct->GetLength();
  static const Int_t kLocScalerLength    = fLocalStruct->GetScalerLength();

  // [16(local)*6 words + 6 words]*8(reg) + 8 words = 824 
  static const Int_t kBufferSize = (16 * (kLocHeaderLength+1) +  (kRegHeaderLength+1))* 8 
      +  kDarcHeaderLength + kGlobalHeaderLength + 2;

  // [16(local)*51 words + 16 words]*8(reg) + 8 + 10 + 8 words scaler event 6682 words
  static const Int_t kScalerBufferSize = (16 * (kLocHeaderLength +  kLocScalerLength +1) +  
					 (kRegHeaderLength + kRegScalerLength +1))* 8 +
                                         (kDarcHeaderLength + kDarcScalerLength + 
					  kGlobalHeaderLength + kGlobalScalerLength + 2);
  if(fScalerEvent) {
    eventPhys = 0; //set to generate scaler events
    fHeader->fWord2 |= (0x1 << 14); // set L1SwC bit on
  }
  if(fScalerEvent)
    buffer = new Int_t [kScalerBufferSize];
  else
    buffer = new Int_t [kBufferSize];

  // reset crate

  // open DDL file, on per 1/2 chamber
  for ( Int_t iDDL = 0; iDDL < 2; ++iDDL ) 
  {
    index = 0; 

    if (iDDL == 0) // suppose global info in DDL one
      globalFlag = 1;
    else 
      globalFlag = 0;

    word = 0;
    // set darc status word
    // see AliMUONDarcHeader.h for details
    AliBitPacking::PackWord((UInt_t)eventPhys,word,30,30);
    AliBitPacking::PackWord((UInt_t)serialNb,word,20,23);
    AliBitPacking::PackWord((UInt_t)globalFlag,word,10,10);
    AliBitPacking::PackWord((UInt_t)version,word,12,19);
    fDarcHeader->SetWord(word);

    memcpy(&buffer[index], fDarcHeader->GetHeader(), (kDarcHeaderLength)*4); 
    index += kDarcHeaderLength;

    // no global input for the moment....
    if (iDDL == 0) {
      fDarcHeader->SetGlobalOutput(gloTrigResp);
      for (Int_t ii = 0; ii < 4; ii++) {
	fDarcHeader->SetGlobalInput(gloTrigInput[ii],ii);
      }
    } else {
      fDarcHeader->SetGlobalOutput(0);
    }

    if (fScalerEvent) {
      // 6 DARC scaler words
      memcpy(&buffer[index], fDarcHeader->GetDarcScalers(),kDarcScalerLength*4);
      index += kDarcScalerLength;
    }
    // end of darc word
    buffer[index++] = fDarcHeader->GetEndOfDarc();

    // 4 words of global board input + Global board output
    memcpy(&buffer[index], fDarcHeader->GetGlobalInput(), (kGlobalHeaderLength)*4); 
    index += kGlobalHeaderLength; 

    if (fScalerEvent) {
      // 10 Global scaler words
      memcpy(&buffer[index], fDarcHeader->GetGlobalScalers(), kGlobalScalerLength*4);
      index += kGlobalScalerLength;
    }

    // end of global word
    buffer[index++] = fDarcHeader->GetEndOfGlobal();
    const AliMpRegionalTrigger* reg = AliMpDDLStore::Instance()->GetRegionalTrigger(); 

    Int_t nCrate = reg->GetNofTriggerCrates()/2;
    // 8 regional cards per DDL
    for (Int_t iReg = 0; iReg < nCrate; ++iReg) {

        // crate info
      AliMpTriggerCrate* crate = AliMpDDLStore::Instance()->GetTriggerCrate(iDDL, iReg);

      if (!crate) 
	AliWarning(Form("Missing crate number %d in DDL %d\n", iReg, iDDL));

      // regional info tree, make sure that no reg card missing
      AliMUONRegionalTrigger* regTrg  = triggerStore.FindRegional(crate->GetId());
      if (!regTrg) 
        AliError(Form("Missing regional board %d in trigger Store\n", crate->GetId()));
    
      // Regional card header
      word = 0;

      // set darc status word
      fRegHeader->SetDarcWord(word);

      regOut    = regTrg->GetOutput();
      regInpLpt = regTrg->GetLocalOutput(0);
      regInpHpt = regTrg->GetLocalOutput(1);

      // fill darc word, not darc status for the moment (empty)
      //see  AliMUONRegHeader.h for details
      AliBitPacking::PackWord((UInt_t)eventPhys,word,31,31); 
      AliBitPacking::PackWord((UInt_t)serialNb,word,20,25); 
      AliBitPacking::PackWord((UInt_t)version,word,8,15);
      AliBitPacking::PackWord((UInt_t)crate->GetId(),word,16,19);
      AliBitPacking::PackWord((UInt_t)regOut,word,0,7); 
      fRegHeader->SetWord(word);


      // fill header later, need local response
      Int_t indexReg = index;
      index += kRegHeaderLength;

      // 11 regional scaler word
      if (fScalerEvent) {
	memcpy(&buffer[index], fRegHeader->GetScalers(), kRegScalerLength*4);
	index += kRegScalerLength;
      }

      // end of regional word
      buffer[index++] = fRegHeader->GetEndOfReg();
      
      // 16 local card per regional board
      //      UShort_t localMask = 0x0;
      
      Int_t nLocalBoard = AliMpConstants::LocalBoardNofChannels();

      for (Int_t iLoc = 0; iLoc < nLocalBoard; iLoc++) {
	  
	// slot zero for Regional card
	Int_t localBoardId = crate->GetLocalBoardId(iLoc);

	if (localBoardId) { // if not empty slot
	  AliMpLocalBoard* localBoard = AliMpDDLStore::Instance()->GetLocalBoard(localBoardId);

	  if (localBoard->IsNotified()) {// if notified board 
	    AliMUONLocalTrigger* locTrg = triggerStore.FindLocal(localBoardId);

	    locCard = locTrg->LoCircuit();
	    locDec  = locTrg->GetLoDecision();
	    trigY   = locTrg->LoTrigY();
	    posY    = locTrg->LoStripY();
	    posX    = locTrg->LoStripX();
	    devX    = locTrg->LoDev();
	    sdevX   = locTrg->LoSdev();
		  
	    AliDebug(4,Form("loctrg %d, posX %d, posY %d, devX %d\n", 
			    locTrg->LoCircuit(),locTrg->LoStripX(),locTrg->LoStripY(),locTrg->LoDev()));  
	    //packing word
	    word = 0;
	    LocalWordPacking(word, (UInt_t)iLoc, (UInt_t)locDec, (UInt_t)trigY, (UInt_t)posY, 
			     (UInt_t)posX, (UInt_t)sdevX, (UInt_t)devX);

	    buffer[index++] = (locTrg->GetX1Pattern() | (locTrg->GetX2Pattern() << 16));
	    buffer[index++] = (locTrg->GetX3Pattern() | (locTrg->GetX4Pattern() << 16));
	    buffer[index++] = (locTrg->GetY1Pattern() | (locTrg->GetY2Pattern() << 16));
	    buffer[index++] = (locTrg->GetY3Pattern() | (locTrg->GetY4Pattern() << 16));
	    buffer[index++] = (Int_t)word; // data word
		      
		
	  }
	  // fill copy card X-Y inputs from the notified cards 
	  if (localBoard->GetInputXfrom() && localBoard->GetInputYfrom()) 
	  {
	    // not triggered
	    locDec = 0; trigY = 1; posY = 15; 	 
	    posX   = 0; devX  = 0; sdevX = 1;
	    LocalWordPacking(word, (UInt_t)iLoc, (UInt_t)locDec, (UInt_t)trigY, (UInt_t)posY, 
			     (UInt_t)posX, (UInt_t)sdevX, (UInt_t)devX);

	    Int_t localFromId = localBoard->GetInputXfrom();
	    AliMUONLocalTrigger* locTrgfrom  = triggerStore.FindLocal(localFromId);

	    buffer[index++] = 0; // copy only X3-4 & Y1-4
	    buffer[index++] = (locTrgfrom->GetX3Pattern() | (locTrgfrom->GetX4Pattern() << 16));
	    buffer[index++] = (locTrgfrom->GetY1Pattern() | (locTrgfrom->GetY2Pattern() << 16));
	    buffer[index++] = (locTrgfrom->GetY3Pattern() | (locTrgfrom->GetY4Pattern() << 16));
	    buffer[index++] = word;
	  }

	} else { 
	  // fill with 10CDEAD word for empty slots
	  for (Int_t i = 0; i < fLocalStruct->GetLength(); i++)
	      buffer[index++] = fLocalStruct->GetDisableWord(); 
	}// condition localBoard
	  
	// 45 regional scaler word
	if (fScalerEvent) {
	  memcpy(&buffer[index], fLocalStruct->GetScalers(), kLocScalerLength*4);
	  index += kLocScalerLength;
	}
	  
	// end of local structure words
	buffer[index++] = fLocalStruct->GetEndOfLocal();
	  
      } // local card 
      // fill regional header with local output
      fRegHeader->SetInput(regInpLpt, 0);
      fRegHeader->SetInput(regInpHpt, 1);
      memcpy(&buffer[indexReg],fRegHeader->GetHeader(),kRegHeaderLength*4);
      
    } // Regional card
    

    // writting onto disk
    // write DDL's
    fHeader->fSize = (index + headerSize) * 4;// total length in bytes
    file[iDDL]->WriteBuffer((char*)fHeader,headerSize*4);
    file[iDDL]->WriteBuffer((char*)buffer,sizeof(int)*index);
  
  }
  delete[] buffer;

  return kTRUE;
}

//____________________________________________________________________
void AliMUONRawWriter::SetScalersNumbers()
{
  /// set numbers for scaler events for trigger headers
  /// since this is provided by the experiment
  /// put dummy numbers to check the monitoring

  fDarcHeader->SetScalersNumbers();
  fRegHeader->SetScalersNumbers();
  fLocalStruct->SetScalersNumbers();
 
  fScalerEvent = kTRUE;
}

