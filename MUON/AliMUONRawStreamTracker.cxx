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
/// This class provides access to MUON digits in raw data.
///
/// It loops over all MUON digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE (under develpment)
/// It can loop also over DDL and store the decoded rawdata in TClonesArray.
/// 
/// First version implement for Tracker
///
///////////////////////////////////////////////////////////////////////////////

#include "AliMUONRawStreamTracker.h"

#include "AliRawReader.h"
#include "AliRawDataHeader.h"
#include "AliLog.h"

#include "AliMUONDspHeader.h"
#include "AliMUONBlockHeader.h"
#include "AliMUONBusStruct.h"
#include "AliMUONDDLTracker.h"

ClassImp(AliMUONRawStreamTracker)

AliMUONRawStreamTracker::AliMUONRawStreamTracker()
  : TObject(),
    fRawReader(0x0),
    fDDL(0),
    fBusPatchId(0),
    fDspId(0),
    fBlkId(0),
    fNextDDL(kTRUE),
    fMaxDDL(20),
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
AliMUONRawStreamTracker::AliMUONRawStreamTracker(AliRawReader* rawReader)
  : TObject(),
    fDDL(0),
    fBusPatchId(0),
    fDspId(0),
    fBlkId(0),
    fNextDDL(kTRUE),
    fMaxDDL(20),
    fMaxBlock(2),
    fMaxDsp(5),
    fMaxBus(5)
{
  //
  // ctor with AliRawReader as argument
  // for reconstruction purpose
  //

  fRawReader       = rawReader;

  fBusPatchManager = new AliMpBusPatch();
  fBusPatchManager->ReadBusPatchFile();

  fDDLTracker      = new AliMUONDDLTracker();
  fBusStruct       = new AliMUONBusStruct();
  fBlockHeader     = new AliMUONBlockHeader();
  fDspHeader       = new AliMUONDspHeader();
}

//_________________________________________________________________
AliMUONRawStreamTracker::AliMUONRawStreamTracker(const AliMUONRawStreamTracker& stream) :
  TObject(stream)
{ 
  //
  // copy ctor
  //
  AliFatal("copy constructor not implemented");
}

//______________________________________________________________________
AliMUONRawStreamTracker& AliMUONRawStreamTracker::operator = (const AliMUONRawStreamTracker& 
					      /* stream */)
{
  // 
  // assignment operator
  //
  AliFatal("assignment operator not implemented");
  return *this;
}


//___________________________________
AliMUONRawStreamTracker::~AliMUONRawStreamTracker()
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

//_____________________________________________________________
Bool_t AliMUONRawStreamTracker::Next()
{
  //
  // read the next raw digit (buspatch structure)
  // returns kFALSE if there is no digit left
  // (under development)

     AliMUONDDLTracker*       ddlTracker = 0x0;
     AliMUONBlockHeader*      blkHeader  = 0x0;
     AliMUONDspHeader*        dspHeader  = 0x0;
     Int_t nBusPatch;
     Int_t nDsp;
     Int_t nBlock;

 next:  
     if (fNextDDL){
       printf("iDDL %d\n", fDDL+1);
       fBlkId = 0;
       fDspId = 0;
       fBusPatchId = 0;
       if(!NextDDL()) 
	 return kFALSE;
     }
     fNextDDL = kFALSE;

     ddlTracker = GetDDLTracker();

     nBlock = ddlTracker->GetBlkHeaderEntries();
     if (fBlkId <  nBlock) {

       blkHeader = ddlTracker->GetBlkHeaderEntry(fBlkId);
       nDsp      = blkHeader->GetDspHeaderEntries();

       if( fDspId < nDsp) {
	 dspHeader = blkHeader->GetDspHeaderEntry(fDspId);
	 nBusPatch = dspHeader->GetBusPatchEntries();

	 if (fBusPatchId < nBusPatch) {
	   fBusStructPtr = dspHeader->GetBusPatchEntry(fBusPatchId++);
	   return kTRUE;

	 } else {// iBusPatch
	   fDspId++;
	   fBusPatchId = 0;
	   goto next;
	   //	Next();
	 }

       } else {// iDsp
	 fBlkId++;
	 fDspId = 0;
	 fBusPatchId = 0;
	 goto next;
	 //      Next();
       }

     } else {// iBlock
       fBlkId = 0;
       fDspId = 0;
       fBusPatchId = 0;
       fNextDDL = kTRUE;
       //return kTRUE;
       goto next; 
     }

     return kFALSE;
}

//______________________________________________________
Bool_t AliMUONRawStreamTracker::NextDDL()
{
  // reading tracker DDL
  // store buspatch info into Array
  // store only non-empty structures (buspatch info with datalength !=0)

  fDDLTracker->GetBlkHeaderArray()->Delete();
  // fDDLTracker->GetBlkHeaderArray()->Clear("C");

  //Read Header Size of DDL,Block,DSP and BusPatch

  Int_t kDDLHeaderSize      = sizeof(AliRawDataHeader)/4;;
  Int_t kBlockHeaderSize    = fBlockHeader->GetHeaderLength();
  Int_t kDspHeaderSize      = fDspHeader->GetHeaderLength();
  Int_t kBusPatchHeaderSize = fBusStruct->GetHeaderLength();

  Int_t totalDDLSize, totalBlockSize, totalDspSize , totalBusPatchSize, dataSize; 


  Int_t iBusPerDSP[5]; // number of bus patches per DSP
  Int_t iDspMax;       // number max of DSP per block

  // minimum data size (only header's)
  Int_t blankDDLSize;
  Int_t blankBlockSize;
  Int_t blankDspSize; 

  if (fDDL >= 20) {
    fDDL = 0;
    return kFALSE;
  }
  AliDebug(3, Form("DDL Number %d\n", fDDL ));

  // getting DSP info
  fBusPatchManager->GetDspInfo(fDDL/2, iDspMax, iBusPerDSP);

  // Each DDL is made with 2 Blocks each of which consists of 5 DSP's at most 
  // and each of DSP has at most 5 buspatches.
  // This information is used to calculate the size of headers (DDL,Block and DSP) 
  // which has no interesting data.
  blankDDLSize   = kDDLHeaderSize + 2*kBlockHeaderSize + 2*iDspMax*kDspHeaderSize;
  blankBlockSize = kBlockHeaderSize + iDspMax*kDspHeaderSize;

  for (Int_t i = 0; i < iDspMax; i++) {
    blankDDLSize   += 2*iBusPerDSP[i]*kBusPatchHeaderSize;
    blankBlockSize +=   iBusPerDSP[i]*kBusPatchHeaderSize;
  }
  fRawReader->Reset();
  fRawReader->Select(0X9, fDDL, fDDL);  //Select the DDL file to be read  

  fRawReader->ReadHeader();

  // 4 is multiplied to convert 2 bytes words
  totalDDLSize = (fRawReader->GetDataSize() + sizeof(AliRawDataHeader))/4; 

  // Compare the DDL header with an empty DDL header size to read the file
  if(totalDDLSize > blankDDLSize) {      

    Int_t totalDataWord = fRawReader->GetDataSize(); // in bytes
    UInt_t *buffer = new UInt_t[totalDataWord/4];
  
    fRawReader->ReadNext((UChar_t*)buffer, totalDataWord); 
 
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

    delete[] buffer;
  } //loop checking the header size of DDL


  fDDL++;

  return kTRUE;
}
//______________________________________________________
void AliMUONRawStreamTracker::ResetDDL()
{
  // reseting TClonesArray
  // after each DDL
  //
  fDDLTracker->GetBlkHeaderArray()->Clear("C");
}

//______________________________________________________
void AliMUONRawStreamTracker::SetMaxDDL(Int_t ddl) 
{
  // set DDL number
  if (ddl > 20) ddl = 20;
  fMaxDDL = ddl;
}

//______________________________________________________
void AliMUONRawStreamTracker::SetMaxBlock(Int_t blk) 
{
  // set regional card number
  if (blk > 2) blk = 2;
  fMaxBlock = blk;
}
