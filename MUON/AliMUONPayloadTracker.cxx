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


///////////////////////////////////////////////////////////////////////////////
///
/// Class Payload
///
/// Decodes rawdata from buffer and stores in TClonesArray.
/// 
/// First version implement for Tracker
///
///////////////////////////////////////////////////////////////////////////////

#include "AliMUONPayloadTracker.h"

#include "AliRawReader.h"
#include "AliRawDataHeader.h"
#include "AliLog.h"

#include "AliMUONDspHeader.h"
#include "AliMUONBlockHeader.h"
#include "AliMUONBusStruct.h"
#include "AliMUONDDLTracker.h"

ClassImp(AliMUONPayloadTracker)

AliMUONPayloadTracker::AliMUONPayloadTracker()
  : TObject(),
    fMaxBlock(2),
    fMaxDsp(5),
    fMaxBus(5)
{
  //
  // create an object to read MUON raw digits
  // Default ctor for monitoring purposes
  //
  fBusPatchManager = new AliMpBusPatch();
  fBusPatchManager->ReadBusPatchFile();

  fDDLTracker      = new AliMUONDDLTracker();
  fBusStruct       = new AliMUONBusStruct();
  fBlockHeader     = new AliMUONBlockHeader();
  fDspHeader       = new AliMUONDspHeader();
}

//_________________________________________________________________
AliMUONPayloadTracker::AliMUONPayloadTracker(const AliMUONPayloadTracker& stream) :
  TObject(stream)
{ 
  //
  // copy ctor
  //
  AliFatal("copy constructor not implemented");
}

//______________________________________________________________________
AliMUONPayloadTracker& AliMUONPayloadTracker::operator = (const AliMUONPayloadTracker& 
					      /* stream */)
{
  // 
  // assignment operator
  //
  AliFatal("assignment operator not implemented");
  return *this;
}


//___________________________________
AliMUONPayloadTracker::~AliMUONPayloadTracker()
{
  //
  // clean up
  //
  delete fBusPatchManager;
  delete fDDLTracker;
  delete fBusStruct;
  delete fBlockHeader;
  delete fDspHeader;
}

//______________________________________________________
Bool_t AliMUONPayloadTracker::Decode(UInt_t* buffer, Int_t ddl)
{
  // reading tracker DDL
  // store buspatch info into Array
  // store only non-empty structures (buspatch info with datalength !=0)


  //Read Header Size of DDL,Block,DSP and BusPatch
  static Int_t kBlockHeaderSize    = fBlockHeader->GetHeaderLength();
  static Int_t kDspHeaderSize      = fDspHeader->GetHeaderLength();
  static Int_t kBusPatchHeaderSize = fBusStruct->GetHeaderLength();

  //  Int_t totalDDLSize;
  Int_t totalBlockSize;
  Int_t totalDspSize;
  Int_t totalBusPatchSize;
  Int_t dataSize; 


  Int_t iBusPerDSP[5]; // number of bus patches per DSP
  Int_t iDspMax;       // number max of DSP per block

  // minimum data size (only header's)
  //  Int_t blankDDLSize;
  Int_t blankBlockSize;
  Int_t blankDspSize; 

  AliDebug(3, Form("DDL Number %d\n", ddl ));

  // getting DSP info
  fBusPatchManager->GetDspInfo(ddl/2, iDspMax, iBusPerDSP);

  // Each DDL is made with 2 Blocks each of which consists of 5 DSP's at most 
  // and each of DSP has at most 5 buspatches.
  // This information is used to calculate the size of headers (Block and DSP) 
  // which has empty data.

  //  blankDDLSize   = 2*kBlockHeaderSize + 2*iDspMax*kDspHeaderSize;
  blankBlockSize = kBlockHeaderSize + iDspMax*kDspHeaderSize;
  // totalDDLSize   = sizeof(buffer)/4;

  for (Int_t i = 0; i < iDspMax; i++) {
    //  blankDDLSize   += 2*iBusPerDSP[i]*kBusPatchHeaderSize;
    blankBlockSize +=   iBusPerDSP[i]*kBusPatchHeaderSize;
  }

  // Compare the DDL header with an empty DDL header size to read the file
  //  if(totalDDLSize > blankDDLSize) {  //should not happen in real life    
 
    // indexes
    Int_t indexDsp;
    Int_t indexBusPatch;
    Int_t index = 0;

    for(Int_t iBlock = 0; iBlock < 2 ;iBlock++){  // loop over 2 blocks

      // copy within padding words
      memcpy(fBlockHeader->GetHeader(),&buffer[index], (kBlockHeaderSize+1)*4);

      totalBlockSize = fBlockHeader->GetTotalLength();

      fDDLTracker->AddBlkHeader(*fBlockHeader);

      if (fBlockHeader->GetPadding() == 0xDEAD) // skipping padding word
	index++;

      if(totalBlockSize > blankBlockSize) {        // compare block header
	index += kBlockHeaderSize;

	for(Int_t iDsp = 0; iDsp < iDspMax ;iDsp++){   //DSP loop

	  if (iDsp > fMaxDsp) break;
	  
	  memcpy(fDspHeader->GetHeader(),&buffer[index], kDspHeaderSize*4);

	  totalDspSize = fDspHeader->GetTotalLength();
	  indexDsp = index;

	  blankDspSize =  kDspHeaderSize + iBusPerDSP[iDsp]*kBusPatchHeaderSize; // no data just header

	  fDDLTracker->AddDspHeader(*fDspHeader, iBlock);

	  if(totalDspSize > blankDspSize) {       // Compare DSP Header
	    index += kDspHeaderSize;
		
	    for(Int_t iBusPatch = 0; iBusPatch < iBusPerDSP[iDsp]; iBusPatch++) {  

	      if (iBusPatch > fMaxBus) break; 

	      //copy buffer into header structure
	      memcpy(fBusStruct->GetBusPatchHeader(), &buffer[index], kBusPatchHeaderSize*4);

	      totalBusPatchSize = fBusStruct->GetTotalLength();
	      indexBusPatch     = index;
		
	      //Check Buspatch header, not empty events
	      if(totalBusPatchSize > kBusPatchHeaderSize) {    

		index   += kBusPatchHeaderSize;
		dataSize = fBusStruct->GetLength();
		Int_t bufSize = fBusStruct->GetBufSize();

		if(dataSize > 0) { // check data present
		  if (dataSize > bufSize) fBusStruct->SetAlloc(dataSize);

		  //copy buffer into data structure
		  memcpy(fBusStruct->GetData(), &buffer[index], dataSize*4);
		  fBusStruct->SetBlockId(iBlock); // could be usefull in future applications ?
		  fBusStruct->SetDspId(iDsp);
		  fDDLTracker->AddBusPatch(*fBusStruct, iBlock, iDsp);

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

    //  } //loop checking the header size of DDL

  return kTRUE;
}

//______________________________________________________
void AliMUONPayloadTracker::ResetDDL()
{
  // reseting TClonesArray
  // after each DDL
  //
  fDDLTracker->GetBlkHeaderArray()->Delete();
}

//______________________________________________________
void AliMUONPayloadTracker::SetMaxBlock(Int_t blk) 
{
  // set regional card number
  if (blk > 2) blk = 2;
  fMaxBlock = blk;
}
