// $Id: AliHLTTPCHWClusterFinderEmulatorComponent.cxx 48710 2011-03-24 12:14:53Z richterm $

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Timm Steinbeck, Matthias Richter                      *
//* Developers:      Kenneth Aamodt <kenneth.aamodt@student.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTPCHWClusterFinderSupport.cxx
    @author Sergey Gorbunov <gorbunov@fias.uni-frankfurt.de>
    @date   
    @brief  Input interface for AliHLTTPCHWClusterFinderEmulator
*/


#include "AliHLTTPCHWClusterFinderSupport.h"
#include "AliHLTDataTypes.h"
#include "AliHLTTPCMapping.h"

#include "AliHLTTPCDigitReaderUnpacked.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCDefinitions.h"
#include "AliRawDataHeader.h"

#if __GNUC__>= 3
using namespace std;
#endif

#include <cstdlib>
#include <cerrno>
#include <sys/time.h>


AliHLTTPCHWClusterFinderSupport::AliHLTTPCHWClusterFinderSupport()
  : 
  AliHLTLogging(),
  fEventMemory(0),
  fEventMCMemory(0)
{
  // see header file for class documentation
}

AliHLTTPCHWClusterFinderSupport::~AliHLTTPCHWClusterFinderSupport()
{
  // see header file for class documentation
  delete[] fEventMemory;
  delete[] fEventMCMemory;
}

AliHLTTPCHWClusterFinderSupport::AliHLTTPCHWClusterFinderSupport(const AliHLTTPCHWClusterFinderSupport&)
  : 
  AliHLTLogging(),
  fEventMemory(0),
  fEventMCMemory(0)
{
  // dummy
}

AliHLTTPCHWClusterFinderSupport& AliHLTTPCHWClusterFinderSupport::operator=(const AliHLTTPCHWClusterFinderSupport&){
  // dummy
  return *this;
}



int AliHLTTPCHWClusterFinderSupport::CreateRawEvent
( const AliHLTComponentBlockData* block, 
  const AliHLTUInt32_t *&rawEvent, AliHLTUInt64_t &rawEventSize32, 
  const AliHLTInt32_t *&mcLabels,  AliHLTUInt64_t &mcLabelsSize32 
  )
{
  // the method creates TPC raw data out of the input block
  // MC labels are provided if possible  
  //

  delete[] fEventMemory;
  delete[] fEventMCMemory;
  
  rawEvent = 0;
  rawEventSize32 = 0;
  mcLabels = 0;
  mcLabelsSize32 = 0;
      
  if( block->fPtr==NULL ){
    HLTWarning("NULL pointer to the data block");
    return 0;
  }

  Int_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *block );
  Int_t patch = AliHLTTPCDefinitions::GetMinPatchNr( *block );
  AliHLTTPCMapping mapping(patch);

  const char *str=Form("slice %d patch %d:", slice, patch);

  if ( block->fDataType == (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC) )
    {    
      // already raw format -> only set the pointers and estimate the size

      // read CDH header, estimate size of the data 
      
      AliHLTUInt64_t headerSize = sizeof(AliRawDataHeader);                  
 
      AliRawDataHeader *cdhHeader = reinterpret_cast<AliRawDataHeader*>(block->fPtr);
      
      AliHLTUInt64_t blockSize = block->fSize; // size of the raw data in bytes      

      if( cdhHeader->fSize!=0xFFFFFFFF ){ // use size information from the header
	blockSize = cdhHeader->fSize;
	if( blockSize > block->fSize ){
	  HLTWarning("%s Could not find a valid DDL header!",str);
	  return 0;
	}
      }
      
      if( blockSize < headerSize ){
	HLTWarning("%s Buffer size is smaller than CDH header size", str);
	return 0;
      }
      
      rawEvent = reinterpret_cast<AliHLTUInt32_t*> (reinterpret_cast<UChar_t*>(block->fPtr)+headerSize);
      rawEventSize32 = ( blockSize - headerSize )/sizeof(AliHLTUInt32_t);

    }
  else if ( block->fDataType == AliHLTTPCDefinitions::fgkUnpackedRawDataType )
    {

      AliHLTTPCDigitReaderUnpacked digitReader; 	 
      digitReader.SetUnsorted(kTRUE);      
      
      if( digitReader.InitBlock(block->fPtr,block->fSize,patch,slice)<0 ) {
	HLTWarning("failed setting up digit reader (InitBlock)");
	return 0;
      }

      int nDigitsTotal = 0;
      int nBunchesTotal = 0;
     
      while( digitReader.NextChannel() ){
	while(digitReader.NextBunch()){
	  nBunchesTotal++;
	  nDigitsTotal+=digitReader.GetBunchSize();
	}
      }
      
      digitReader.Reset();

      if( digitReader.InitBlock(block->fPtr,block->fSize,patch,slice)<0) {
	HLTWarning("failed setting up digit reader (InitBlock)");
	return 0;
      }
      
      Int_t nPadsTotal = 0;
      Int_t firstRow = AliHLTTPCTransform::GetFirstRow(patch);
      Int_t nRows = AliHLTTPCTransform::GetNRows(patch);

      for( int i=0; i<nRows; i++ ){
	nPadsTotal += AliHLTTPCTransform::GetNPads(firstRow+i);	
      }

      int totalSize32 = (nDigitsTotal + nBunchesTotal*2)/3+2*nPadsTotal + 10;
      int totalSizeMC32 = 3*nDigitsTotal + 10;
      
      fEventMemory = new AliHLTUInt32_t[totalSize32];
      if( !fEventMemory ){
	HLTWarning("Not enougth memory: can not allocate %d bytes",totalSize32*8);
	return 0;
      }

      fEventMCMemory = new AliHLTInt32_t[totalSizeMC32];
      if( !fEventMCMemory ){
	HLTWarning("Not enougth memory: can not allocate %d bytes",totalSizeMC32*8);
	delete[] fEventMemory;
	return 0;
      }

      AliHLTUInt32_t nWords32 = 0;
      int mcIndex = 0;
      int err=0;

      AliHLTTPCDigitData tmpDigit;
      tmpDigit.fTrackID[0] = -1;
      tmpDigit.fTrackID[1] = -1;
      tmpDigit.fTrackID[2] = -1;

      while( !err && digitReader.NextChannel() ){
	
	Int_t row=digitReader.GetRow();
	Int_t pad=digitReader.GetPad();

	// create header

	if( nWords32 >= totalSize32){ err = 1; break; }	
	
	AliHLTUInt32_t *header = fEventMemory + nWords32;
	nWords32++;

	int seek10 = 2;
	int prevTime = 10000000;
	int nWords10 = 0;
	while(digitReader.NextBunch()){
	  
	  Int_t nSignals = digitReader.GetBunchSize();
	  if( nSignals <=0 ){
	    HLTWarning("Empty bunch received");
	    continue;
	  }

	  Int_t time = digitReader.GetTime() + nSignals-1;
	  
	  if( time-nSignals+1<0 || time>=AliHLTTPCTransform::GetNTimeBins() ){
	    HLTWarning("Wrong time bins received: %d-%d for row %d pad %d", time-nSignals+1, time, row, pad);
	    break;
	  }

	  if( time >= prevTime ){
	    HLTWarning("Unexpected order of TPC bunches in row %d, pad %d", row, pad);	    
	    break;
	  }

	  prevTime = time-nSignals+1;
	  
	  if( nWords32+( 2+nSignals)/3+1 >= totalSize32 ){ err = 1; break; }
	  if( mcIndex + nSignals*3   >= totalSizeMC32 ){ err = 1; break; }

	  if( nWords10 + 2 + nSignals > 0x2FF ){
	    HLTWarning("Too much data in row %d, pad %d", row, pad);	    
	    break;
	  }

	  nWords10 += 2 + nSignals;

	  Add10Word( nWords32, seek10, nSignals + 2 );
	  Add10Word( nWords32, seek10, time );
	  
	  const UInt_t *bunchData = digitReader.GetSignals();
	  const AliHLTTPCDigitData *mcDigits = digitReader.GetBunchDigits();
	  if( !mcDigits ){
	    HLTWarning("No MC labels found for a bunch of digits");
	  }

	  for(Int_t i=nSignals-1; i>=0; i--){
	    Add10Word( nWords32, seek10, bunchData[i] );	    
	    const AliHLTTPCDigitData &digit = mcDigits ?mcDigits[i] :tmpDigit;
	    fEventMCMemory[mcIndex++] = digit.fTrackID[0];
	    fEventMCMemory[mcIndex++] = digit.fTrackID[1];
	    fEventMCMemory[mcIndex++] = digit.fTrackID[2];	    	    
	  }	  

	} // bunches
	
      *header = (1<<30) | ((nWords10&0x2FF)<<16) | (mapping.GetHwAddress(row, pad) & 0xFFF);

      }// channels (pads)

      if( err ){
	HLTError("Internal error: too less memory allocated");	
      } else {
	rawEvent = fEventMemory;
	rawEventSize32 = nWords32;
	mcLabels = fEventMCMemory;
	mcLabelsSize32 = mcIndex;
      }

    } // unpacked data type

  return 0;
}


void AliHLTTPCHWClusterFinderSupport::Add10Word( AliHLTUInt32_t &nWords32, int &seek10, UInt_t data )
{
  // add 10-bit data to the 32-bit word
  // fEventMemory [nWords32] --- current 32-bit word
  // *seek10 --- 10-bit position withing the word
  // pointers are increased, a new word is first initialised to 0

  data = data & 0x3FF; // truncate to 10 bits
  AliHLTUInt32_t *seek = fEventMemory + nWords32;

  if( seek10 == 2 ){
    nWords32++;
    fEventMemory[nWords32-1] = data<<20;
    seek10 = 1;
  } else if( seek10 == 1 ){
    fEventMemory[nWords32-1] &= 0xFFF003FF;
    fEventMemory[nWords32-1] |= (data<<10);
    seek10 = 0;
  } else if( seek10 == 0 ){
    fEventMemory[nWords32-1] &= 0xFFFFFC00;
    fEventMemory[nWords32-1] |= data;
    seek10 = 2;
  } 
}



int AliHLTTPCHWClusterFinderSupport::CheckRawData( const AliHLTUInt32_t *buffer,
						   unsigned long bufferSize32, int patch, int slice )
{
  //
  // The procedure checks consistency of the data
  //

  const unsigned int headerSize32 = 8;

  if (!buffer) return 0;

  const char *str=Form("slice %d patch %d:", slice, patch);
  
  if( bufferSize32 < headerSize32 ){
    HLTWarning("%s Buffer size is smaller than CDH header size", str);
    return kFALSE;
  }    
  
  // read data header 
 
  AliHLTUInt32_t blockSize32 = bufferSize32; // size of the raw data in words

  if( buffer[0]!=0xFFFFFFFF ) blockSize32 = buffer[0]/4; // use size information from the header  
  if( blockSize32 > bufferSize32 ){  
    HLTWarning(Form("%s Could not find a valid DDL header!",str));
    return 0;
  }
  
  UChar_t rcuVer = (UChar_t)( (buffer[1] >> 24) & 0xFF ); 

  if (rcuVer < 2) {
    HLTWarning("%s Old data format, RCU version %d", str,rcuVer);
    return 0;
  }

  // is the block valid
  //AliHLTUInt32_t blockAttributes = buffer[3]; // block attributes (bits 24-31) and participating sub detectors 
  //cout<<blockAttributes<<" "<<(blockAttributes >> 24)<<endl;
  //if ( !( (blockAttributes >> 24) & 1) ) return 0; 
     

  const AliHLTUInt32_t* fData = buffer + headerSize32;
  unsigned long  dataSize32 = blockSize32 - headerSize32;       
  
  // Read the RCU trailer according to the RCU formware version specified in CDH
  // Cross-check with version found in the trailer
  // The two major bit should be 11 (identifies the end of the trailer)    
  
  AliHLTUInt32_t word = fData[dataSize32 - 1];
  
  if ((word >> 30) != 3) {
    HLTWarning("%s Last RCU trailer word not found!", str);
    return 0;
  }
  
  UChar_t ver = (word >> 16) & 0xFF;
  Int_t rcuId = (Int_t)((word >> 7) & 0x1FF);
  Int_t rcuTrailerSize32 = (word & 0x7F); // size of RCU trailer data in words
  
  if (ver != rcuVer) {
    HLTWarning("%s Wrong RCU firmware version detected: %d != %d",
	       str,ver,rcuVer);
    return 0;
  }  

  if (rcuTrailerSize32 < 2) {
    HLTWarning(Form("Invalid trailer size found (%d bytes) !",
		    rcuTrailerSize32*4));
    return 0;
  }
  
  if( rcuTrailerSize32 > dataSize32 ){
    HLTWarning(Form("%s Invalid trailer size found (%d bytes) ! The size is bigger than the raw data size (%d bytes)!",
		    str, rcuTrailerSize32*4,dataSize32*4));
    return 0;    
  }

  // check the trailer

  Int_t trailerIndex = dataSize32 - rcuTrailerSize32;

  for( int i=trailerIndex; i<dataSize32-1; i++){
    if ((fData[i] >> 30) != 2) {
      HLTWarning("%s Missing RCU trailer identifier pattern!",str);
      continue;
    }
  }

  // Read the payload size  
 
  Int_t  rcuPayloadSize32 = fData[trailerIndex] & 0x3FFFFFF;

  if ( rcuPayloadSize32 + rcuTrailerSize32  != dataSize32) {
    HLTWarning(Form("%s Inconsistent raw data size ! Raw data size - %d bytes (from CDH), RCU trailer - %d bytes, raw data size (from RCU trailer) - %d bytes !",
		    str, dataSize32*4,
		    (rcuTrailerSize32)*4,
		    rcuPayloadSize32*4));
    return 0;
  }
    
  
  AliHLTTPCMapping *fMapping = new AliHLTTPCMapping(patch);
  const int kMaxNTimeBins = 2000;

  UShort_t  *channelData10 = new UShort_t[kMaxNTimeBins];    // cache for the decoded altro payload

  Int_t position = 0; // current position (32-bit words) in fData
 
  while(1){

    // Search for the next Altro channel

    word = 0;
    while( position < rcuPayloadSize32 ){
      word = fData[position++];
      if( (word >> 30) == 1) break;
    }
    if (position >= rcuPayloadSize32 ) break; // no next channel found
    
    // extract channel payload and hw address

    Int_t channelPayloadSize10 = (word >> 16) & 0x3FF; // payload size in 10-bit words 
    Int_t channelPayloadSize32 = (channelPayloadSize10+2)/3;
    Bool_t channelIOErrors = (word >> 29) & 0x1; // check for readout errors    
    Short_t  channelHWAddress = word & 0xFFF;

    if( position + channelPayloadSize32-1> rcuPayloadSize32 ){
      HLTWarning(Form("%s Inconsistent channel payload data size: expected <= %d bytes from RCU trailer, found %d bytes in the channel header!",
		      str,(rcuPayloadSize32 - position)*4, channelPayloadSize32*4 ));
      continue;
    }

    bool channelBranchAB = ((channelHWAddress >> 11) & 0x1);
    int channelFEC       = ((channelHWAddress >> 7) & 0xF); // front-end card index
    int channelAltro = ((channelHWAddress >> 4) & 0x7); // altro chip index
    int channelIndex = (channelHWAddress & 0xF); // channel index
    int channelRow = fMapping->GetRow(channelHWAddress);
    int channelPad = fMapping->GetPad(channelHWAddress);

    // Now unpack the Altro data: 10-bit words to 16 bit-words    
    
    Int_t channelData10Index = 0;// current position in the payload

    for (Int_t iword = 0; iword < channelPayloadSize32; iword++) {
      word = fData[position++];
      if ((word >> 30) != 0) {
	HLTWarning(Form("%s Unexpected end of payload in altro channel payload! Address=0x%x, word=0x%x",
			str, channelHWAddress,word));
	channelIOErrors = 1;
	position--;
	break;
      }
      channelData10[channelData10Index++] = (word >> 20) & 0x3FF;
      channelData10[channelData10Index++] = (word >> 10) & 0x3FF;
      channelData10[channelData10Index++] = word & 0x3FF;
    }
  
    if ( channelIOErrors ) continue;    

    // read bunches

    Int_t prevTimeBin =  1024;
    channelData10Index = 0;

    while(1){
      
      // Read next Altro bunch 
  
      if ((channelData10Index+1 >= channelPayloadSize10) ) break;
    
      Int_t bunchLength = channelData10[channelData10Index++];
      Int_t bunchStartTimeBin = channelData10[channelData10Index++];

      if (bunchLength <= 2) {
	// Invalid bunch size
	HLTWarning(Form("%s Too short bunch length (%d) in Address=0x%x!",
			str, bunchLength,channelHWAddress));	
	break;
      }
      if( channelData10Index + bunchLength - 2 > channelPayloadSize10 ){
	// Too long bunch detected
	HLTWarning(Form("%s Too long bunch detected in Address=0x%x! Expected <= %d 10-bit words, found %d !",
			str,channelHWAddress,channelPayloadSize10-channelData10Index,bunchLength));
	break;
      }
            
      if( bunchStartTimeBin-bunchLength+1 < 0) {
	HLTWarning(Form("%s Invalid start time-bin in Address=0x%x ! (%d-%d+1) < 0",
			str,channelHWAddress,bunchStartTimeBin,bunchLength));
	break;
      }
      if (bunchStartTimeBin >= prevTimeBin) {
	HLTWarning(Form("%s Invalid start time-bin in Address=0x%x! (%d>=%d)",
			str,channelHWAddress,bunchStartTimeBin,prevTimeBin));
	break;
      }

      prevTimeBin = bunchStartTimeBin-bunchLength+1;
  
      bunchLength -= 2;
	
      UShort_t* bunchData = &channelData10[channelData10Index];   // pointer to the current bunch samples
      channelData10Index += bunchLength;            
    }
  }
  
  delete[] channelData10;
  delete fMapping;

  return 1;
}

