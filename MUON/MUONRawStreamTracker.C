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

//
// Macro for reading tracker raw data
// Ch. Finck, Subatech Febuary
//

#if !defined(__CINT__) || defined(__MAKECINT__)

// RAW includes
#include "AliRawReaderDate.h"
#include "AliRawReaderFile.h"
#include "AliRawReaderRoot.h"

// MUON includes
#include "AliMUONRawStreamTracker.h"
#include "AliMUONDspHeader.h"
#include "AliMUONBlockHeader.h"
#include "AliMUONBusStruct.h"
#include "AliMUONDDLTracker.h"

#endif


// Macro to read rawdata for tracker
// Ch. Finck, Subatech, April. 2006
//
// This macro is interface with AliRawReader for RAW
// The different stucture of the patload are readout and stored in TClonesArray
// with AliMUONRawStreamTracker classe.
// The macro just simpy read again the TClonesArray contents.
// The parameter of each structure could be seen in the container classes
// AliMUONBlockHeader, AliMUONBlockHeader, AliMUONBusStruct.
// The class AliMUONDDLTracker manages the structure containers.
// The number of structures in the rawdata file could be set.

void MUONRawStreamTracker(Int_t maxEvent = 1000, Int_t minDDL = 0, Int_t maxDDL = 19, TString fileName = "./")
{
   AliRawReader* rawReader = 0x0;

   // check extention to choose the rawdata file format
   if (fileName.EndsWith("/")) {
     rawReader = new AliRawReaderFile(fileName); // DDL files
   } else if (fileName.EndsWith(".root")) {
     rawReader = new AliRawReaderRoot(fileName);
   } else if (!fileName.IsNull()) {
     rawReader = new AliRawReaderDate(fileName); // DATE file
   }

   // raw stream
   AliMUONRawStreamTracker* rawStream  = new AliMUONRawStreamTracker(rawReader);    

   // set the number of DDL block Dsp & buspatch structures that are PRESENT in the rawdata file
   // it's NOT the number to be read.
   // default wise set to 20, 2, 5 ans 5 respectively.
   //    rawStream->SetMaxDDL(xx);
   //    rawStream->SetMaxBlock(xx);
   //    rawStream->SetMaxDsp(xx);
   //    rawStream->SetMaxBus(xx);

   // containers
   AliMUONDDLTracker*       ddlTracker = 0x0;
   AliMUONBlockHeader*      blkHeader  = 0x0;
   AliMUONDspHeader*        dspHeader  = 0x0;
   AliMUONBusStruct*        busStruct  = 0x0;

   //   Loop over events  
   Int_t iEvent = 0;
   Int_t dataSize;

   while (rawReader->NextEvent()) {

     if (iEvent == maxEvent)
       break;

     printf("Event %d\n",iEvent++);

     // read DDL while < 20 DDL
     while(rawStream->NextDDL()) {

       if (rawStream->GetDDL() < minDDL || rawStream->GetDDL() > maxDDL)
	 continue;

       printf("\niDDL %d\n", rawStream->GetDDL());

       ddlTracker =  rawStream->GetDDLTracker();

       // loop over block structure
       Int_t nBlock = ddlTracker->GetBlkHeaderEntries();
       for(Int_t iBlock = 0; iBlock < nBlock ;iBlock++){

	 blkHeader = ddlTracker->GetBlkHeaderEntry(iBlock);
	 printf("Block Total length %d\n",blkHeader->GetTotalLength());

	 // loop over DSP structure
	 Int_t nDsp = blkHeader->GetDspHeaderEntries();
	 for(Int_t iDsp = 0; iDsp < nDsp ;iDsp++){   //DSP loop

	   dspHeader =  blkHeader->GetDspHeaderEntry(iDsp);
	   //   printf("Dsp length %d even word %d\n",dspHeader->GetTotalLength(), dspHeader->GetEventWord());

	   // loop over BusPatch structure
	   Int_t nBusPatch = dspHeader->GetBusPatchEntries();
	   for(Int_t iBusPatch = 0; iBusPatch < nBusPatch; iBusPatch++) {  

	     busStruct = dspHeader->GetBusPatchEntry(iBusPatch);

	     printf("busPatchId %d", busStruct->GetBusPatchId());
	     printf(" BlockId %d", busStruct->GetBlockId());
	     printf(" DspId %d\n", busStruct->GetDspId());

	   // loop over data
	     dataSize = busStruct->GetLength();
	     for (Int_t iData = 0; iData < dataSize; iData++) {

	       Int_t  manuId    = busStruct->GetManuId(iData);
	       Int_t  channelId = busStruct->GetChannelId(iData);
	       Int_t  charge    = busStruct->GetCharge(iData);
	       printf("manuId: %d, channelId: %d charge: %d\n", manuId, 
		      channelId, charge);

	     } // iData
	   } // iBusPatch
	 } // iDsp
       } // iBlock
     } // NextDDL
   }// NextEvent

   delete rawReader;
   delete rawStream;
}
