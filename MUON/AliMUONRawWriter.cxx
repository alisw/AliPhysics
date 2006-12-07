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

/// \class AliMUONRawWriter
/// MUON Raw Data generaton in ALICE-MUON
/// This class version 3 (further details could be found in Alice-note)
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

#include "AliMUONRawWriter.h"

#include "AliMUONBlockHeader.h"
#include "AliMUONBusStruct.h"
#include "AliMUONConstants.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONData.h"
#include "AliMUONDigit.h"
#include "AliMUONDspHeader.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONRegionalTrigger.h"
#include "AliMUONRegHeader.h"
#include "AliMUONTriggerCrate.h"
#include "AliMUONTriggerCrateStore.h"

#include "AliMpBusPatch.h"
#include "AliMpDEManager.h"
#include "AliMpExMap.h"
#include "AliMpConstants.h"
#include "AliMpPlaneType.h"
#include "AliMpSegmentation.h"
#include "AliMpStationType.h"
#include "AliMpVSegmentation.h"

#include "AliRawReader.h"
#include "AliBitPacking.h" 
#include "AliDAQ.h"
#include "AliLog.h"

#include "TList.h"
#include "TObjArray.h"
#include "TStopwatch.h"

/// \cond CLASSIMP
ClassImp(AliMUONRawWriter) // Class implementation in ROOT context
/// \endcond

Int_t AliMUONRawWriter::fgManuPerBusSwp1B[12]  = {1, 27, 53, 79, 105, 131, 157, 183, 201, 214, 224, 232};
Int_t AliMUONRawWriter::fgManuPerBusSwp1NB[12] = {1, 27, 53, 79, 105, 131, 157, 183, 201, 214, 225, 233};

Int_t AliMUONRawWriter::fgManuPerBusSwp2B[12]  = {1, 27, 53, 79, 105, 131, 157, 183, 201, 214, 226, 246};
Int_t AliMUONRawWriter::fgManuPerBusSwp2NB[12] = {1, 27, 53, 79, 105, 131, 157, 183, 201, 214, 227, 245};

namespace 
{
  enum ETimer { kWriteTracker, kWriteTrigger, kDigitLoop, kGetBusPatch, kTest, kLast };
}

//__________________________________________________________________________
AliMUONRawWriter::AliMUONRawWriter(AliMUONData* data)
  : TObject(),
    fMUONData(data),
    fBlockHeader(new AliMUONBlockHeader()),
    fDspHeader(new AliMUONDspHeader()),
    fDarcHeader(new AliMUONDarcHeader()),
    fRegHeader(new AliMUONRegHeader()),
    fLocalStruct(new AliMUONLocalStruct()),
    fBusPatchManager(new AliMpBusPatch()),
    fCrateManager(new AliMUONTriggerCrateStore()),
    fScalerEvent(kFALSE),
    fHeader(),
    fTimers(new TStopwatch[kLast])

{
  /// Standard Constructor

  AliDebug(1,"Standard ctor");
  fFile[0] = fFile[1] = 0x0;  
  fFile[2] = fFile[3] = 0x0;  

  // setting data key to default value (only for writting)
  fBlockHeader->SetDataKey(fBlockHeader->GetDefaultDataKey());
  fDspHeader->SetDataKey(fDspHeader->GetDefaultDataKey());

  // bus patch managers
  fBusPatchManager->ReadBusPatchFile();

  // Crate manager
  fCrateManager->ReadFromFile();

  // timers
  for ( Int_t i = 0; i < kLast; ++i )
  {
    fTimers[i].Start(kTRUE); 
    fTimers[i].Stop();
  }
  
}

//__________________________________________________________________________
AliMUONRawWriter::AliMUONRawWriter()
  : TObject(),
    fMUONData(0),
    fBlockHeader(0),
    fDspHeader(0),
    fDarcHeader(0),
    fRegHeader(0),
    fLocalStruct(0),
    fBusPatchManager(0),
    fCrateManager(0x0),
    fScalerEvent(kFALSE),
    fHeader(),
    fTimers(0)
{
  /// Default Constructor

  AliDebug(1,"Default ctor");   
  fFile[0] = fFile[1] = 0x0;  
  fFile[2] = fFile[3] = 0x0;  

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

  delete fBusPatchManager;
  delete fCrateManager;

  for ( Int_t i = 0; i < kLast; ++i )
  {
    AliDebug(1, Form("Execution time (timer %d) : R:%7.2fs C:%7.2fs",i,
                 fTimers[i].RealTime(),fTimers[i].CpuTime()));
  }
  
  delete[] fTimers;
}

//______________________________________________________________________________
//void
//AliMUONRawWriter::CheckDigits()
//{
//  std::map<int,std::map<int,int> > m;
//  
//  for (Int_t iSt = 0; iSt < AliMUONConstants::NTrackingCh()/2; ++iSt) 
//  {
//    for (Int_t iCh = iSt*2; iCh <= iSt*2 + 1; ++iCh) 
//    {      
//      TClonesArray* muonDigits = fMUONData->Digits(iCh);
//      for (Int_t idig = 0; idig < muonDigits->GetEntriesFast(); idig++) 
//      {        
//        AliMUONDigit* digit = (AliMUONDigit*) muonDigits->UncheckedAt(idig);
//        Int_t busPatchId = GetBusPatch(*digit);
//        m[busPatchId][digit->ManuId()]++;
//      }
//    } 
//  }
//  
//  std::map<int,std::map<int,int> >::const_iterator it;
//  
//  Int_t nManuMax(0);
//  
//  for ( it = m.begin(); it != m.end(); ++it )
//  {
//    AliDebug(1,Form("BusPatch %3d has %3d manus",it->first,it->second.size()));
//    nManuMax = std::max((Int_t)it->second.size(),nManuMax);
//    std::map<int,int>::const_iterator it2;
//    for ( it2 = it->second.begin(); it2 != it->second.end(); ++it2 )
//    {
//      AliDebug(1,Form("        BusPatch %3d Manu %4d Nch %3d",it->first,it2->first,it2->second));
//    }
//  }
//  AliDebug(1,Form("Max manus per busPatch : %3d",nManuMax));
//}

//____________________________________________________________________
Int_t AliMUONRawWriter::Digits2Raw()
{
  /// convert digits of the current event to raw data

  Int_t idDDL;
  Char_t name[255];

  fMUONData->GetLoader()->LoadDigits("READ");

  fMUONData->SetTreeAddress("D,GLT");

  fMUONData->ResetDigits();
  fMUONData->ResetTrigger();
  
  // This will get both tracker and trigger digits.
  fMUONData->GetDigits();
  
//  CheckDigits();

  // tracking chambers
  
  for (Int_t iSt = 0; iSt < AliMUONConstants::NTrackingCh()/2; ++iSt) {

    // open files for one station
    // cos station 3, 1/4 of DE's from 2 chambers has same DDL number 
    idDDL = iSt * 4;
    strcpy(name,AliDAQ::DdlFileName("MUONTRK",idDDL));
    fFile[0] = fopen(name,"w");

    idDDL = (iSt * 4) + 1;
    strcpy(name,AliDAQ::DdlFileName("MUONTRK",idDDL));
    fFile[1] = fopen(name,"w");

    idDDL =  (iSt * 4) + 2;;
    strcpy(name,AliDAQ::DdlFileName("MUONTRK",idDDL));
    fFile[2] = fopen(name,"w");

    idDDL =  (iSt * 4) + 3;
    strcpy(name,AliDAQ::DdlFileName("MUONTRK",idDDL));
    fFile[3] = fopen(name,"w");

    WriteTrackerDDL(iSt);
  
    // reset and close when station has been processed
    fclose(fFile[0]);
    fclose(fFile[1]);
    fclose(fFile[2]);
    fclose(fFile[3]);
     
  }
 
  AliDebug(1,"Tracker written");
  
  // trigger chambers
 
  // open files
  idDDL = 0;// MUTR
  strcpy(name,AliDAQ::DdlFileName("MUONTRG",idDDL));
  fFile[0] = fopen(name,"w");

  idDDL = 1;// MUTR
  strcpy(name,AliDAQ::DdlFileName("MUONTRG",idDDL));
  fFile[1] = fopen(name,"w");

   WriteTriggerDDL();
  
  // reset and close
  fclose(fFile[0]);
  fclose(fFile[1]);

  AliDebug(1,"Trigger written");

  fMUONData->ResetDigits();
  fMUONData->ResetTrigger();  
  fMUONData->GetLoader()->UnloadDigits();

  AliDebug(1,"muondata reset");
  
  return kTRUE;
}

//____________________________________________________________________
Int_t AliMUONRawWriter::WriteTrackerDDL(Int_t iSt)
{
  /// writing DDL for tracker
  /// used inverse mapping

  fTimers[kWriteTracker].Start(kFALSE);
  
  static const Int_t kMAXADC = (1<<12)-1; // We code the charge on a 12 bits ADC.

  // resets
  TClonesArray* muonDigits = 0;

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

  // buffer size (max'ed out)
  // (((43 manus max per bus patch *64 channels + 4 bus patch words) * 5 bus patch 
  //   + 10 dsp words)*5 dsps + 8 block words)*2 blocks 
  static const Int_t kBufferSize = (((43*64 + 4)*5 + 10)*5 + 8)*2;
  
  Int_t nDigits;

  AliMpExMap busPatchMap(kTRUE);
  
  fTimers[kDigitLoop].Start(kFALSE);
  
  for (Int_t iCh = iSt*2; iCh <= iSt*2 + 1; ++iCh) {

    muonDigits = fMUONData->Digits(iCh);
    
    nDigits = muonDigits->GetEntriesFast();
    
    // loop over digit
    for (Int_t idig = 0; idig < nDigits; ++idig) {
      
      AliMUONDigit* digit = static_cast<AliMUONDigit*>(muonDigits->UncheckedAt(idig));
      
      charge = digit->ADC();
      if ( charge > kMAXADC )
      {
        // This is most probably an error in the digitizer (which should insure
        // the adc is below kMAXADC), so make it a (non-fatal) error indeed.
        AliError(Form("adc value %d above %x for ch %d . Setting to %x. Digit is:",iCh,
                      charge,kMAXADC,kMAXADC));
        StdoutToAliError(digit->Print());
        charge = kMAXADC;
      }
      
      // inverse mapping
      fTimers[kGetBusPatch].Start(kFALSE);
      busPatchId = GetBusPatch(*digit);
      fTimers[kGetBusPatch].Stop();
      if (busPatchId<0) continue;
      
      if ( digit->ManuId() > 0x7FF || digit->ManuId() < 0 ||
           digit->ManuChannel() > 0x3F || digit->ManuChannel() < 0 )
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
        parity ^=  ((word >> i) & 0x1);
      AliBitPacking::PackWord((UInt_t)parity,word,31,31);
      
      AliMUONBusStruct* busStruct = 
        static_cast<AliMUONBusStruct*>(busPatchMap.GetValue(busPatchId));
      
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
      
    } // idig
  } // loop over chamber in station
    
  fTimers[kDigitLoop].Stop();
  
  // getting info for the number of buspatches
  Int_t iBusPatch;
  Int_t length;
  Int_t iBusPerDSP[5];//number of bus patches per DSP
  Int_t iDspMax; //number max of DSP per block
  Int_t iFile = 0;

  AliMUONBusStruct* busStructPtr(0x0);

  // open DDL files, 4 per station
  for (Int_t iDDL = iSt*4; iDDL < 4 + iSt*4; ++iDDL) {

    fBusPatchManager->ResetBusItr(iDDL);
    fBusPatchManager->GetDspInfo(iDDL, iDspMax, iBusPerDSP);

    Int_t buffer[kBufferSize];
    
    totalDDLLength = 0;

    indexBlk = 0;
    indexDsp = 0;
    index = 0;

    // two blocks A and B per DDL
    for (Int_t iBlock = 0; iBlock < 2; ++iBlock) {
      
      // block header
      length = fBlockHeader->GetHeaderLength();
      memcpy(&buffer[index],fBlockHeader->GetHeader(),length*4);
      indexBlk = index;
      index += length; 
      
      // 5 DSP's max per block
      for (Int_t iDsp = 0; iDsp < iDspMax; ++iDsp) {
        
        // DSP header
        length = fDspHeader->GetHeaderLength();
        memcpy(&buffer[index],fDspHeader->GetHeader(),length*4);
        indexDsp = index;
        index += length; 
        
        // 5 buspatches max per DSP
        for (Int_t i = 0; i < iBusPerDSP[iDsp]; i++) {
          
          iBusPatch = fBusPatchManager->NextBusInDDL(iDDL);
          
          // iteration over bus patch in DDL
          if (iBusPatch == -1) {
            AliWarning(Form("Error in bus itr in DDL %d\n", iDDL));
            continue;
          }
          
          // 4 DDL's per station, condition needed for station 3
          iFile = iDDL - iSt*4; // works only if DDL begins at zero (as it should be) !!!
          
          busStructPtr = static_cast<AliMUONBusStruct*>(busPatchMap.GetValue(iBusPatch));
          
          // check if buspatchid has digit
          if (busStructPtr) {
            // add bus patch structure header
            length = busStructPtr->GetHeaderLength();
            memcpy(&buffer[index],busStructPtr->GetHeader(),length*4);
            index += length;
            
            // add bus patch data
            length = busStructPtr->GetLength();
            memcpy(&buffer[index],busStructPtr->GetData(),length*4);
            index += length;
            
            if (AliLog::GetGlobalDebugLevel() == 3) {
              for (Int_t j = 0; j < busStructPtr->GetLength(); j++) {
                printf("busPatchId %d, manuId %d channelId %d\n", busStructPtr->GetBusPatchId(), 
                       busStructPtr->GetManuId(j), busStructPtr->GetChannelId(j));
              }
            }
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
    // write DDL 1 - 4
    // total length in bytes
    fHeader.fSize = (totalDDLLength + headerSize) * 4;
      
    fwrite((char*)(&fHeader),headerSize*4,1,fFile[iFile]);
    fwrite(buffer,sizeof(int),index,fFile[iFile]);
  }
  
  fTimers[kWriteTracker].Stop();
  return kTRUE;
}

//____________________________________________________________________
Int_t AliMUONRawWriter::GetBusPatch(Int_t detElemId, Int_t manuId) const
{
  /// Determine the BusPatch this digit belongs to.
  
  Int_t* ptr = 0;
    
  AliMpPlaneType plane = 
    (manuId & AliMpConstants::ManuMask(kNonBendingPlane)) ? 
    kNonBendingPlane : kBendingPlane; 
  
  AliMpStationType stationType = AliMpDEManager::GetStationType(detElemId);
  
  if ( stationType == kStation1)
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
  else if ( stationType == kStation2)
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
  
  Int_t m = ( manuId & 0x3FF ); // remove bit 10
                                        //FIXME : how can we remove that condition
                                        // on the 10-th bit ? All the rest need not any knowledge about it,
                                        // can't we find a way to get manu<->buspatch transparent to this too ?
  
  if ( stationType == kStation1 || stationType == kStation2 )
  {
    for (pos = 11; pos >=0 ; --pos)
    {
      if (m >= ptr[pos]) break;
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
                  pos, (Int_t)vec->GetSize(), manuId, detElemId));
    AliError(Form("Chamber %s Plane %s manuId %d m %d",
                  StationTypeName(stationType).Data(),
                  PlaneTypeName(plane).Data(),
                  manuId,
                  m));
    return -1;
  }
  
  Int_t busPatchId = vec->At(pos);
  
  if ( ( stationType == kStation1 || stationType == kStation2 ) &&
       ( plane == kNonBendingPlane ) )
  {
    busPatchId += 12;
  }
  
  return busPatchId;
}

//____________________________________________________________________
Int_t AliMUONRawWriter::GetBusPatch(const AliMUONDigit& digit) const
{
  /// Determine the BusPatch this digit belongs to.

  return GetBusPatch(digit.DetElemId(),digit.ManuId());
}

//____________________________________________________________________
Int_t AliMUONRawWriter::WriteTriggerDDL()
{
  /// Write trigger DDL

  fTimers[kWriteTrigger].Start(kFALSE);
  
 // DDL event one per half chamber

  // stored local id number 
  TArrayI isFired(256);
  isFired.Reset();


 // DDL header size
  Int_t headerSize = sizeof(AliRawDataHeader)/4;

  TClonesArray* localTrigger;
  TClonesArray* globalTrigger;
  TClonesArray* regionalTrigger;

  AliMUONGlobalTrigger* gloTrg;
  AliMUONLocalTrigger* locTrg = 0x0;
  AliMUONRegionalTrigger* regTrg = 0x0;

  // global trigger for trigger pattern
  globalTrigger = fMUONData->GlobalTrigger(); 
  gloTrg = (AliMUONGlobalTrigger*)globalTrigger->UncheckedAt(0);
  if (!gloTrg) 
  {
    fTimers[kWriteTrigger].Stop();
    return 0;
  }
  
  Int_t gloTrigResp = gloTrg->GetGlobalResponse();

  // local trigger 
  localTrigger = fMUONData->LocalTrigger();   


  // regional trigger
  regionalTrigger = fMUONData->RegionalTrigger();   


  UInt_t word;
  Int_t* buffer = 0;
  Int_t index;
  Int_t iEntries = 0;
  Int_t iLocCard, locCard;
  UChar_t locDec, trigY, posY, posX, regOut;
  UInt_t regInpLpt;
  UInt_t regInpHpt;

  UInt_t devX;
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
  if(fScalerEvent)
    eventPhys = 0; //set to generate scaler events

  Int_t nEntries = (Int_t) (localTrigger->GetEntries());// 234 local cards
  // stored the local card id that's fired
  for (Int_t i = 0; i <  nEntries; i++) {
    locTrg = (AliMUONLocalTrigger*)localTrigger->At(i);
    isFired[locTrg->LoCircuit()] = 1; // storing local boards with informations
  }

  if (!nEntries)
    AliDebug(1, "No Trigger information available");

  if(fScalerEvent)
    buffer = new Int_t [kScalerBufferSize];
  else
    buffer = new Int_t [kBufferSize];

  // reset crate

  // open DDL file, on per 1/2 chamber
  for (Int_t iDDL = 0; iDDL < 2; iDDL++) {

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
    if (iDDL == 0)
     fDarcHeader->SetGlobalOutput(gloTrigResp);
    else 
     fDarcHeader->SetGlobalOutput(0);

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
      memcpy(fDarcHeader->GetGlobalScalers(), &buffer[index], kGlobalScalerLength*4);
      index += kGlobalScalerLength;
    }

    // end of global word
    buffer[index++] = fDarcHeader->GetEndOfGlobal();

    // 8 regional cards per DDL
    for (Int_t iReg = 0; iReg < 8; iReg++) {

      // crate info
      AliMUONTriggerCrate* crate = fCrateManager->Crate(iDDL, iReg);

      if (!crate) 
	AliWarning(Form("Missing crate number %d in DDL %d\n", iReg, iDDL));

      // regional info tree, make sure that no reg card missing
      for (Int_t i = 0; i < 16; ++i) {
	regTrg  = (AliMUONRegionalTrigger*)regionalTrigger->At(i);
	if (regTrg)
	  if (regTrg->GetId() == (iReg + iDDL*8)) break;
      }

      // Regional card header
      word = 0;

      // set darc status word
      fRegHeader->SetDarcWord(word);

      regOut    = regTrg->GetOutput();
      regInpHpt = regTrg->GetLocalOutput(0);
      regInpLpt = regTrg->GetLocalOutput(1);

      // fill darc word, not darc status for the moment (empty)
      //see  AliMUONRegHeader.h for details
      AliBitPacking::PackWord((UInt_t)eventPhys,word,31,31); 
      AliBitPacking::PackWord((UInt_t)serialNb,word,19,24); 
      AliBitPacking::PackWord((UInt_t)version,word,16,23);
      AliBitPacking::PackWord((UInt_t)iReg,word,15,18);
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
      
      TObjArray *boards = crate->Boards();


      // 16 local card per regional board
      //      UShort_t localMask = 0x0;

      for (Int_t iLoc = 0; iLoc < 16; iLoc++) {

	// slot zero for Regional card
	AliMUONLocalTriggerBoard* localBoard = (AliMUONLocalTriggerBoard*)boards->At(iLoc+1);

	if (localBoard) { // if not empty slot

	  if ((iLocCard = localBoard->GetNumber()) != 0) {// if notified board

	    if (isFired[iLocCard]) { // if card has triggered
	      locTrg  = (AliMUONLocalTrigger*)localTrigger->At(iEntries++);
	      locCard = locTrg->LoCircuit();
	      locDec  = locTrg->GetLoDecision();
	      trigY = 0;
	      posY  = locTrg->LoStripY();
	      posX  = locTrg->LoStripX();
	      devX  = locTrg->LoDev();

	      AliDebug(4,Form("loctrg %d, posX %d, posY %d, devX %d\n", 
			      locTrg->LoCircuit(),locTrg->LoStripX(),locTrg->LoStripY(),locTrg->LoDev()));
	    } else { //no trigger (see PRR chpt 3.4)
	      locDec = 0;
	      trigY = 1;
	      posY = 15;
	      posX = 0;
	      devX = 0x8;
	      // set local card id to -1
	      locCard = -1; 
	    }
	   
	    //packing word
	    word = 0;
	    AliBitPacking::PackWord((UInt_t)iLoc,word,19,22); //card id number in crate
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

	    } else {
	      buffer[index++] = 0; // 4 words for x1, x2, y1, y2
	      buffer[index++] = 0; 
	      buffer[index++] = 0; 
	      buffer[index++] = 0; 
	      buffer[index++] = (Int_t)word; // data word

	    }
	  } else {// number!=0
	  // fill with 10CDEAD word for 'non-notified' slots
	  for (Int_t i = 0; i < fLocalStruct->GetLength(); i++)
	    buffer[index++] = fLocalStruct->GetDisableWord(); 
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
      fRegHeader->SetInput(regInpHpt, 0);
      fRegHeader->SetInput(regInpHpt, 1);
      memcpy(&buffer[indexReg],fRegHeader->GetHeader(),kRegHeaderLength*4);

    } // Regional card
    

    // writting onto disk
    // write DDL's
    fHeader.fSize = (index + headerSize) * 4;// total length in bytes
    fwrite((char*)(&fHeader),headerSize*4,1,fFile[iDDL]);
    fwrite(buffer,sizeof(int),index,fFile[iDDL]);
  
  }
  delete[] buffer;

  fTimers[kWriteTrigger].Stop();
  
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

