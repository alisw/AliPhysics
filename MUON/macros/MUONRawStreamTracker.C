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

// $Id$

/// \ingroup macros
/// \file MUONRawStreamTracker.C
/// \brief Macro for reading tracker raw data
///
/// \author Ch. Finck, Subatech Febuary
///
/// Added example routines to show how to use the interface of the high
/// performance decoder AliMUONRawStreamTrackerHP.
/// by  Artur Szostak <artursz@iafrica.com>
///
/// This macro is interfaced with AliRawReader for RAW
///
/// There are 2 ways of reading the data: 
/// - one where each intermediate structure (block, dsp, buspatch) is looped over
/// - and one, using an iterator, where we're directly accessing the pad informations 
/// (charge).
///


#if !defined(__CINT__) || defined(__MAKECINT__)

// RAW includes
#include "AliRawReader.h"

// MUON includes
#include "AliMUONRawStreamTrackerHP.h"
#include "TStopwatch.h"
#include "AliRawDataErrorLog.h"
#include "Riostream.h"
#include "AliMUONLogger.h"

#endif

void MUONRawStreamTrackerExpert(TString fileName = "./", Int_t maxEvent = 1000,  
                                Int_t minDDL = 0, Int_t maxDDL = 19)
{
  /// This routine shows how to use the decoder's expert interface.
  
  TStopwatch timer;
  timer.Start(kTRUE);
  
  AliRawReader* rawReader = AliRawReader::Create(fileName.Data());
  
  // raw stream
  AliMUONRawStreamTrackerHP rawStream(rawReader);
  
  // light weight interfaces to headers
  const AliMUONRawStreamTrackerHP::AliBlockHeader*      blkHeader  = 0x0;
  const AliMUONRawStreamTrackerHP::AliDspHeader*        dspHeader  = 0x0;
  const AliMUONRawStreamTrackerHP::AliBusPatch*         busStruct  = 0x0;
  
  //   Loop over events  
  Int_t iEvent = 0;
  Int_t dataSize;
  
  while (rawReader->NextEvent()) {
    
    if (iEvent == maxEvent)
      break;
    
    printf("Event %d\n",iEvent++);
    
    // read DDL while < 20 DDL
    while(rawStream.NextDDL()) {
      
      if (rawStream.GetDDL() < minDDL || rawStream.GetDDL() > maxDDL)
        continue;
      
      printf("\niDDL %d\n", rawStream.GetDDL());
      
      // loop over block structure
      Int_t nBlock = rawStream.GetBlockCount();
      for(Int_t iBlock = 0; iBlock < nBlock ;iBlock++){
        
        blkHeader = rawStream.GetBlockHeader(iBlock);
        printf("Block %d Total length %d\n",iBlock,blkHeader->GetTotalLength());
       
        // loop over DSP structure
        Int_t nDsp = rawStream.GetDspCount(iBlock);
        for(Int_t iDsp = 0; iDsp < nDsp ;iDsp++){   //DSP loop
          
          dspHeader =  blkHeader->GetDspHeader(iDsp);
          printf("Dsp %d length %d error word %d\n",iDsp,dspHeader->GetTotalLength(), dspHeader->GetErrorWord());
          
          // loop over BusPatch structure
          Int_t nBusPatch = rawStream.GetBusPatchCount(iBlock, iDsp);
          for(Int_t iBusPatch = 0; iBusPatch < nBusPatch; iBusPatch++) {  
            
            busStruct = dspHeader->GetBusPatch(iBusPatch);
            
            // loop over data
            dataSize = busStruct->GetLength();
            for (Int_t iData = 0; iData < dataSize; iData++) {
              
              Int_t  manuId    = busStruct->GetManuId(iData);
              Int_t  channelId = busStruct->GetChannelId(iData);
              Int_t  charge    = busStruct->GetCharge(iData);
              printf("buspatch %5d manuI %4d channel %3d charge %4d\n", 
                     busStruct->GetBusPatchId(),
                     manuId, 
                     channelId, charge);
            } // iData
          } // iBusPatch
        } // iDsp
      } // iBlock
    } // NextDDL
  }// NextEvent
  
  delete rawReader;
  timer.Print();
}


void MUONRawStreamTrackerExpert2(TString fileName = "./", Int_t maxEvent = 1000,  
                                Int_t minDDL = 0, Int_t maxDDL = 19)
{
  /// This routine shows an alternate way to iterate over the DDL structures
  /// compared to MUONRawStreamTrackerExpert().
  
  TStopwatch timer;
  timer.Start(kTRUE);
  
  AliRawReader* rawReader = AliRawReader::Create(fileName.Data());
  
  // raw stream
  AliMUONRawStreamTrackerHP rawStream(rawReader);
  
  // light weight interfaces to headers
  const AliMUONRawStreamTrackerHP::AliBlockHeader*      blkHeader  = 0x0;
  const AliMUONRawStreamTrackerHP::AliDspHeader*        dspHeader  = 0x0;
  const AliMUONRawStreamTrackerHP::AliBusPatch*         busStruct  = 0x0;
  
  //   Loop over events  
  Int_t iEvent = 0;
  Int_t dataSize;
  
  while (rawReader->NextEvent()) {
    
    if (iEvent == maxEvent)
      break;
    
    printf("Event %d\n",iEvent++);
    
    // read DDL while < 20 DDL
    while(rawStream.NextDDL()) {
      
      if (rawStream.GetDDL() < minDDL || rawStream.GetDDL() > maxDDL)
        continue;
      
      printf("\niDDL %d\n", rawStream.GetDDL());
      
      // loop over block structure
      Int_t nBlock = rawStream.GetBlockCount();
      for(Int_t iBlock = 0; iBlock < nBlock ;iBlock++){
        
        blkHeader = rawStream.GetBlockHeader(iBlock);
        printf("Block %d Total length %d\n",iBlock,blkHeader->GetTotalLength());
       
        // loop over DSP structure
        Int_t nDsp = rawStream.GetDspCount(iBlock);
        for(Int_t iDsp = 0; iDsp < nDsp ;iDsp++){   //DSP loop
          
          dspHeader =  rawStream.GetDspHeader(iBlock, iDsp);
          printf("Dsp %d length %d error word %d\n",iDsp,dspHeader->GetTotalLength(), dspHeader->GetErrorWord());
          
          // loop over BusPatch structure
          Int_t nBusPatch = rawStream.GetBusPatchCount(iBlock, iDsp);
          for(Int_t iBusPatch = 0; iBusPatch < nBusPatch; iBusPatch++) {  
            
            busStruct = rawStream.GetBusPatch(iBlock, iDsp, iBusPatch);
            
            // loop over data
            dataSize = busStruct->GetLength();
            for (Int_t iData = 0; iData < dataSize; iData++) {
              
              Int_t  manuId    = busStruct->GetManuId(iData);
              Int_t  channelId = busStruct->GetChannelId(iData);
              Int_t  charge    = busStruct->GetCharge(iData);
              printf("buspatch %5d manuI %4d channel %3d charge %4d\n", 
                     busStruct->GetBusPatchId(),
                     manuId, 
                     channelId, charge);
            } // iData
          } // iBusPatch
        } // iDsp
      } // iBlock
    } // NextDDL
  }// NextEvent
  
  delete rawReader;
  timer.Print();
}


void MUONRawStreamTrackerExpert3(TString fileName = "./", Int_t maxEvent = 1000,  
                                Int_t minDDL = 0, Int_t maxDDL = 19)
{
  /// This routine shows yet another alternate way to iterate over the DDL
  /// structures compared to MUONRawStreamTrackerExpert().
  
  TStopwatch timer;
  timer.Start(kTRUE);
  
  AliRawReader* rawReader = AliRawReader::Create(fileName.Data());
  
  // raw stream
  AliMUONRawStreamTrackerHP rawStream(rawReader);
  
  // light weight interfaces to headers
  const AliMUONRawStreamTrackerHP::AliBlockHeader*      blkHeader  = 0x0;
  const AliMUONRawStreamTrackerHP::AliDspHeader*        dspHeader  = 0x0;
  const AliMUONRawStreamTrackerHP::AliBusPatch*         busStruct  = 0x0;
  
  //   Loop over events  
  Int_t iEvent = 0;
  Int_t dataSize;
  
  while (rawReader->NextEvent()) {
    
    if (iEvent == maxEvent)
      break;
    
    printf("Event %d\n",iEvent++);
    
    // read DDL while < 20 DDL
    while(rawStream.NextDDL()) {
      
      if (rawStream.GetDDL() < minDDL || rawStream.GetDDL() > maxDDL)
        continue;
      
      printf("\niDDL %d\n", rawStream.GetDDL());
      
      // loop over block structure
      Int_t iBlock = 0;
      blkHeader = rawStream.GetFirstBlockHeader();
      while (blkHeader != NULL)
      {
        printf("Block %d Total length %d\n",iBlock,blkHeader->GetTotalLength());
       
        // loop over DSP structure
        Int_t iDsp = 0;
        dspHeader = blkHeader->GetFirstDspHeader();
        while (dspHeader != NULL)
        {
          printf("Dsp %d length %d error word %d\n",iDsp,dspHeader->GetTotalLength(), dspHeader->GetErrorWord());
          
          // loop over BusPatch structure
          Int_t iBusPatch = 0;
          busStruct = dspHeader->GetFirstBusPatch();
          while (busStruct != NULL)
          {
            // loop over data
            dataSize = busStruct->GetLength();
            for (Int_t iData = 0; iData < dataSize; iData++) {
              
              Int_t  manuId    = busStruct->GetManuId(iData);
              Int_t  channelId = busStruct->GetChannelId(iData);
              Int_t  charge    = busStruct->GetCharge(iData);
              printf("buspatch %5d manuI %4d channel %3d charge %4d\n", 
                     busStruct->GetBusPatchId(),
                     manuId, 
                     channelId, charge);
            } // iData
            busStruct = busStruct->Next();
            iBusPatch++;
          } // iBusPatch
          dspHeader = dspHeader->Next();
          iDsp++;
        } // iDsp
        blkHeader = blkHeader->Next();
        iBlock++;
      } // iBlock
    } // NextDDL
  }// NextEvent
  
  delete rawReader;
  timer.Print();
}


void MUONRawStreamTrackerSimple(TString fileName = "./", Int_t maxEvent = 1000)
{
  /// This routine shows how to use the high performance decoder's simple interface.

  TStopwatch timer;
  timer.Start(kTRUE);
  
  AliRawReader* rawReader = AliRawReader::Create(fileName.Data());
  
  // raw stream
  AliMUONRawStreamTrackerHP rawStream(rawReader);    
  
  //   Loop over events  
  Int_t iEvent = 0;
  
  while (rawReader->NextEvent()) {
    
    if (iEvent == maxEvent)
      break;
    
    printf("Event %d\n",iEvent++);
    
    Int_t busPatch;
    UShort_t manuId, adc;
    UChar_t manuChannel;
    
    rawStream.First();
    
    while ( rawStream.Next(busPatch,manuId,manuChannel,adc) )
    {      
      printf("buspatch %5d manuI %4d channel %3d charge %4d\n", 
             busPatch,manuId,manuChannel, adc);
    }
  }
  
  delete rawReader;
  timer.Print();
}


void ShowErrors(const AliRawReader& rawReader)
{
  for ( Int_t i = 0; i < rawReader.GetNumberOfErrorLogs(); ++i )
  {
    AliRawDataErrorLog* error = rawReader.GetErrorLog(i);
    error->Print();
  }

  cout << Form("Number of error logs : %d (%d events)",
               rawReader.GetNumberOfErrorLogs(),
               rawReader.GetNumberOfEvents()) << endl;
}


void MUONRawStreamTrackerSimple2(TString fileName = "./", Int_t maxEvent = 1000)
{
  /// This routine is an alternative to MUONRawStreamTrackerSimple() which is even faster.

  TStopwatch timer;
  timer.Start(kTRUE);
  
  AliRawReader* rawReader = AliRawReader::Create(fileName.Data());
  
  // raw stream
  AliMUONRawStreamTrackerHP rawStream(rawReader);
  rawStream.EnableRawReaderErrorLogger();
  
  //   Loop over events  
  Int_t iEvent = 0;
  
  while (rawReader->NextEvent()) {
    
    if (iEvent == maxEvent)
      break;
    
    printf("Event %d\n",iEvent++);
    
    UShort_t manuId, adc;
    UChar_t manuChannel;
    
    rawStream.First();
    const AliMUONRawStreamTrackerHP::AliBusPatch* buspatch = NULL;
    while ((buspatch = rawStream.Next()) != NULL)
    {
      for (UInt_t i = 0; i < buspatch->GetDataCount(); i++)
      {
        buspatch->GetData(i, manuId, manuChannel, adc);
        printf("buspatch %5d manuI %4d channel %3d charge %4d\n", 
               buspatch->GetBusPatchId(), manuId, manuChannel, adc);
      }
    }
  }

  ShowErrors(*rawReader);
  
  delete rawReader;
  timer.Print();
}

void MUONRawStreamTrackerErrorCount(TString fileName = "collection://filelist", Int_t maxEvent = -1)
{
  /// This routine is just a loop to get the error log at the end
  
  TStopwatch timer;
  timer.Start(kTRUE);
  
  AliRawReader* rawReader = AliRawReader::Create(fileName.Data());
  
  // raw stream
  AliMUONRawStreamTrackerHP rawStream(rawReader);
  rawStream.DisableWarnings();    

  AliMUONLogger logger;
  
  rawStream.EnableMUONErrorLogger();  
  rawStream.SetMUONErrorLogger(&logger);    
  rawStream.SetLoggingDetailLevel(AliMUONRawStreamTrackerHP::kMediumErrorDetail);
  
  //   Loop over events  
  Int_t iEvent = 0;
  
  while (rawReader->NextEvent()) 
  {    
    if (iEvent == maxEvent) break;
    
    rawStream.First();

    const AliMUONRawStreamTrackerHP::AliBusPatch* buspatch = NULL;
    
    while ((buspatch = rawStream.Next()) != NULL)
    {
    }
    
    ++iEvent;
  }
  
  logger.Print();
  
//  ShowErrors(*rawReader);

  delete rawReader;
  timer.Print();
}

