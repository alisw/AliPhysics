// $Id$
//****************************************************************************
//* This file is property of and copyright by the ALICE HLT Project          * 
//* ALICE Experiment at CERN, All rights reserved.                           *
//*                                                                          *
//* Primary Authors: Sergey Gorbunov, Torsten Alt                            *
//* Developers:      Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de> *
//*                  Torsten Alt <talt@cern.ch>                              *
//*                  for The ALICE HLT Project.                              *
//*                                                                          *
//* Permission to use, copy, modify and distribute this software and its     *
//* documentation strictly for non-commercial purposes is hereby granted     *
//* without fee, provided that the above copyright notice appears in all     *
//* copies and that both the copyright notice and this permission notice     *
//* appear in the supporting documentation. The authors make no claims       *
//* about the suitability of this software for any purpose. It is            *
//* provided "as is" without express or implied warranty.                    *
//****************************************************************************

//  @file   AliHLTTPCHWCFSupport.cxx
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @brief  Input interfaces for FPGA ClusterFinder Emulator for TPC
//  @brief  ( see AliHLTTPCHWCFEmulator class )
//  @note


#include "AliHLTTPCHWCFSupport.h"
#include "AliHLTDataTypes.h"
#include "AliHLTTPCMapping.h"
#include "AliHLTTPCDigitReaderUnpacked.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCDefinitions.h"
#include "AliRawDataHeader.h"
#include "AliHLTTPCHWCFEmulator.h"
#include "AliTPCcalibDB.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include "TMath.h"

#if __GNUC__>= 3
using namespace std;
#endif

#include <cstdlib>
#include <algorithm>
#include <cerrno>
#include <sys/time.h>


AliHLTTPCHWCFSupport::AliHLTTPCHWCFSupport()
  : 
  AliHLTLogging(),
  fEventMemory(0),
  fEventMCMemory(0)
{
  // see header file for class documentation
  for( int i=0; i<fgkNSlices; i++ )
    for( int j=0; j<fgkNPatches; j++ ) fMapping[i][j] = 0;
}


AliHLTTPCHWCFSupport::~AliHLTTPCHWCFSupport()
{
  // see header file for class documentation
  for( int i=0; i<fgkNSlices; i++ )
    for( int j=0; j<fgkNPatches; j++ ) delete[] fMapping[i][j];
  ReleaseEventMemory(); 
}

AliHLTTPCHWCFSupport::AliHLTTPCHWCFSupport(const AliHLTTPCHWCFSupport&)
  : 
  AliHLTLogging(),
  fEventMemory(0),
  fEventMCMemory(0)
{
  // dummy
}

AliHLTTPCHWCFSupport& AliHLTTPCHWCFSupport::operator=(const AliHLTTPCHWCFSupport&){
  // dummy
  return *this;
}


void AliHLTTPCHWCFSupport::ReleaseEventMemory()
{
  // clean up 
  if( fEventMemory ) delete[] fEventMemory;
  if( fEventMCMemory )delete[] fEventMCMemory;
  fEventMemory = 0;
  fEventMCMemory = 0;
}


const AliHLTUInt32_t *AliHLTTPCHWCFSupport::GetMapping( int slice, int patch )
{ 
  // see header file for class documentation
  if( slice<0 || slice>=fgkNSlices ){
    HLTFatal("Wrong slice number %d, no mapping is provided.", slice);
    return 0;
  }
  if( patch<0 || patch>= fgkNPatches ){
    HLTFatal("Wrong patch number %d, no mapping is provided.", patch);
    return 0;
  }
  if( !fMapping[slice][patch] ) fMapping[slice][patch] = ReadMapping(slice,patch);
return fMapping[slice][patch];
}


AliHLTUInt32_t *AliHLTTPCHWCFSupport::ReadMapping( int slice, int patch, const char *mappingFileName ) const
{
  // Create mapping array for one patch 
  // If no mapping file provided, reads from default file
  // Output: mapping [] array of type AliHLTUInt32_t, where :
  //
  // mapping[0] == N hardware adresses in the array (mapping size is maping[0] + 1 )
  // mapping[hwAddress] == configWord
  //
  // configWord consist of:
  //
  // bits 0-7: pad number
  // bits 8-13: row number
  // bit  14 : flag for border pad
  // bit  15 : is the pad active
  // bits 16->28 : gain calibration as 13 bit fixed point,
  //               with 1 bit position before decimal point

  const AliHLTUInt32_t  kBorderFlag = (1 << 14); 
  const AliHLTUInt32_t  kActiveFlag = (1 << 15); 
  
  if( slice<0 || slice>=fgkNSlices ){
     HLTFatal("Wrong slice number %d, no mapping is provided.", slice);
     return 0;
  }

  if( patch<0 || patch>5 ){
     HLTFatal("Wrong patch number %d, no mapping is provided.", patch);
     return 0;
  }

  // AliHLTTPCTransform::GetFirstRow returns first row in scheme A.
  // We have to transform to scheme B by AliHLTTPCTransform::Slice2Sector.

  UInt_t offsetSchemeB=0;
  Int_t sector = 0;
  {
    Int_t tmp=0;
    AliHLTTPCTransform::Slice2Sector(slice, AliHLTTPCTransform::GetFirstRow(patch),
				     sector, tmp);
    offsetSchemeB = (UInt_t) tmp;
  }
  

  AliTPCcalibDB *calib = AliTPCcalibDB::Instance();  
  AliTPCCalPad * gainTPC = 0;
  AliTPCCalROC * gainROC = 0;
  if( calib ) gainTPC = calib->GetPadGainFactor();
  if( gainTPC ) gainROC = gainTPC->GetCalROC(sector);  // pad gains per given sector
  else{      
    HLTWarning("No TPC gain calibration found");  
  }

  TString filename;
  
  if( mappingFileName ){
    filename = mappingFileName;
  } else {
    const char* basePath=getenv("ALICE_ROOT");
    if (basePath) filename.Form("%s/TPC/mapping/Patch%d.data", basePath,patch);    
  } 
  
  ifstream inFile;
  inFile.open(filename.Data());
  if (!inFile) {
    HLTFatal("Unable to open mapping file: %s   This means no mapping is provided.", filename.Data());
    return 0;
  }


  AliHLTUInt32_t *mapping = 0; 
  AliHLTUInt32_t *rowBranchPadHw = 0;
  bool err = 1;
  do{

    UInt_t nHWAdd=0;
    UInt_t maxHWAdd=0;

    if( !(inFile >> nHWAdd ) || !(inFile >> maxHWAdd)  ){
      HLTError("Mapping file for patch %d corrupted &s", patch,filename.Data());
      break;
    }

    if( maxHWAdd > 0xFFF ){
      HLTError("Max hardware address exceeded for patch %d, max number is %d, number from mapping file is %d.",patch, 0xFFF, maxHWAdd+1);	
      break;
    }

    if(nHWAdd > maxHWAdd ){
      HLTError("Too large number of hardware addresses for patch %d: max number is %d, number from mapping file is %d.",patch, maxHWAdd, nHWAdd );
      break;
    }
      
    mapping = new AliHLTUInt32_t[maxHWAdd+2];
    rowBranchPadHw = new AliHLTUInt32_t[nHWAdd];
    if( !mapping || !rowBranchPadHw ){
      HLTError("Can not allocate &d bytes of memory", (maxHWAdd+1+nHWAdd)*sizeof(AliHLTUInt32_t));
      break;
    }

    for( unsigned int i=0; i<maxHWAdd+2; i++ ) mapping[i] = 0;
    for( unsigned int i=0; i<nHWAdd; i++ ) rowBranchPadHw[i] = 0;    
    mapping[0] = maxHWAdd+1;
    UInt_t nRead = 0;
    err = 0;
    while(!err ){
      UInt_t hwAdd=0;
      UInt_t row=0;
      UInt_t pad=0;
      if( !(inFile>>hwAdd) || !(inFile>>row) || !(inFile>>pad) ) break;      

      err = 1;

      if ( nRead >= nHWAdd ){
	HLTError("Too many hardware addresses: %d, expected %d, mapping file %s corrupted?", nRead+1,  nHWAdd, filename.Data());
	break;
      }
      if (hwAdd>maxHWAdd) {
	HLTError("hardware address exceeds max hwAddress %d, mapping file %s corrupted?", maxHWAdd, filename.Data());	
	break;
      }

      if( row < offsetSchemeB ){
	HLTError("row number %d below minimum %d for patch %d, mapping file %s corrupted?", row, offsetSchemeB, patch, filename.Data());	
	break;	  
      }	

      row -= offsetSchemeB;
	
      if( row > 0x3F ){
	HLTError("row number %d withing patch exceed the maximum %d for patch %d, mapping file %s corrupted?", row, 0x3F, patch, filename.Data());	
	break;	  
      }

      if( pad > 0xFF ){
	HLTError("pad number %d exceed the maximum %d for patch %d, mapping file %s corrupted?", pad, 0xFF, patch, filename.Data());	
	break;	  
      }

      bool active = true; // Currently all channels are always active	
      //

      AliHLTFloat64_t gain = 1.;
      if( gainROC ){
	gain = gainROC->GetValue(offsetSchemeB+row,pad);
	if( gain>1.e-4 ) gain = 1./gain;
	else gain = 0;
      }
      gain*= (1<<12);
      AliHLTUInt32_t  gainCalib = TMath::Nint(gain); 
      if( gainCalib > 0x1FFF ) gainCalib = 0x1FFF;

      AliHLTUInt32_t configWord = ( (row & 0x3F) << 8 ) | (pad & 0xFF);
      if ( active ) configWord |= kActiveFlag;
      configWord |= (gainCalib & 0x1FFF) << 16;	

      mapping[1+hwAdd] = configWord;
	
      AliHLTUInt32_t branch = (hwAdd >> 11) & 0x1;	
      rowBranchPadHw[nRead] = (row<<25) | (branch<<24) | (pad<<16) | hwAdd;

      nRead++;
      err = 0;
    }
    
    if( err ) break;
    
    if ( nRead!= nHWAdd ){
      HLTError("Too less hardware addresses: %d, expected %d, mapping file %s corrupted?", nRead,  nHWAdd, filename.Data());
      err = 1;
      break;
    }
    
    // mark pads at borders of A/B branches 
      
    std::sort(rowBranchPadHw, rowBranchPadHw + nHWAdd);
    int rowBranchLast = -1;
    for( unsigned int i=0; i<nHWAdd; i++ ){
      int rowBranch = rowBranchPadHw[i]>>24;
      if( rowBranch != rowBranchLast ){
	mapping[1+(rowBranchPadHw[i] & 0xFFF)] |= kBorderFlag;
	rowBranchLast = rowBranch;
	if( i>0 ) mapping[1+(rowBranchPadHw[i-1] & 0xFFF)] |= kBorderFlag;	  
      }
    }
    mapping[1+(rowBranchPadHw[nRead-1] & 0xFFF)] |= kBorderFlag;
    
  } while(0);
  
  inFile.close();

  delete[] rowBranchPadHw;

  if( err ){
    delete[] mapping; 
    return 0;
  }
  return mapping;
}


int AliHLTTPCHWCFSupport::CreateRawEvent
( const AliHLTComponentBlockData* block, 
  const AliHLTUInt32_t *&rawEvent, AliHLTUInt32_t &rawEventSize32, 
  const AliHLTTPCClusterMCLabel *&mcLabels,  AliHLTUInt32_t &nMCLabels 
)
{
  // the method creates TPC raw data out of the input block
  // MC labels are provided if possible  
  //

  ReleaseEventMemory();
  
  rawEvent = 0;
  rawEventSize32 = 0;
  mcLabels = 0;
  nMCLabels = 0;
      
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

      AliHLTUInt32_t totalSize32 = (nDigitsTotal + nBunchesTotal*2)/3+2*nPadsTotal + 10;
      AliHLTUInt32_t totalNMC = nDigitsTotal + 10;
      
      fEventMemory = new AliHLTUInt32_t[totalSize32];
      if( !fEventMemory ){
	HLTWarning("Not enougth memory: can not allocate %d bytes",totalSize32*8);
	return 0;
      }

      fEventMCMemory = new AliHLTTPCClusterMCLabel[totalNMC];
      if( !fEventMCMemory ){
	HLTWarning("Not enougth memory: can not allocate %d bytes",totalNMC*sizeof(AliHLTTPCClusterMCLabel));
	delete[] fEventMemory;
	fEventMemory = 0;
	return 0;
      }

      AliHLTUInt32_t nWords32 = 0;
      AliHLTUInt32_t mcIndex = 0;
      int err=0;

      AliHLTTPCDigitData tmpDigit;
      tmpDigit.fTrackID[0] = -1;
      tmpDigit.fTrackID[1] = -1;
      tmpDigit.fTrackID[2] = -1;

      while( !err && digitReader.NextChannel() ){
	
	Int_t row=digitReader.GetRow();
	Int_t pad=digitReader.GetPad();

	AliHLTUInt32_t hwAddr = mapping.GetHwAddress(row, pad);

	// create header

	if( nWords32 >= totalSize32){ err = 1; break; }	
	
	AliHLTUInt32_t *header = fEventMemory + nWords32;
	nWords32++;

	int seek10 = 2;
	int prevTime = 10000000;
	int nWords10 = 0;
	while(digitReader.NextBunch()){

	  if( hwAddr > 0xFFF ) continue;

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
	  if( mcIndex + nSignals   >= totalNMC ){ err = 1; break; }

	  if( nWords10 + 2 + nSignals > 0x3FF ){
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

	  for(Int_t is=nSignals-1; is>=0; is--){
	    Add10Word( nWords32, seek10, bunchData[is] );	    
	    const AliHLTTPCDigitData &digit = mcDigits ?mcDigits[is] :tmpDigit;
	    int nmc = 0;
	    for( int i=0; i<3; i++ ) if( digit.fTrackID[i] >=0 ) nmc++;      
	    for( int i=0; i<3; i++ ){
	      fEventMCMemory[mcIndex].fClusterID[i].fMCID = digit.fTrackID[i];
	      fEventMCMemory[mcIndex].fClusterID[i].fWeight = 0;
	      if( digit.fTrackID[i] >=0 ){      		
		fEventMCMemory[mcIndex].fClusterID[i].fWeight  = ((float)bunchData[is])/nmc;
	      }
	    }	  
	    mcIndex++;
	  }	  

	} // bunches
	
      *header = (1<<30) | ((nWords10&0x3FF)<<16) | (hwAddr & 0xFFF);

      }// channels (pads)

      if( err ){
	HLTError("Internal error: too less memory allocated");	
      } else {
	for( AliHLTUInt32_t i=0; i<nWords32; i++ ) fEventMemory[i] = AliHLTTPCHWCFEmulator::WriteBigEndian(fEventMemory[i]);
	rawEvent = fEventMemory;
	rawEventSize32 = nWords32;
	mcLabels = fEventMCMemory;
	nMCLabels = mcIndex;
      }

    } // unpacked data type

  return 0;
}


void AliHLTTPCHWCFSupport::Add10Word( AliHLTUInt32_t &nWords32, int &seek10, UInt_t data )
{
  // add 10-bit data to the 32-bit word
  // fEventMemory [nWords32] --- current 32-bit word
  // *seek10 --- 10-bit position withing the word
  // pointers are increased, a new word is first initialised to 0

  data = data & 0x3FF; // truncate to 10 bits

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



int AliHLTTPCHWCFSupport::CheckRawData( const AliHLTUInt32_t *buffer,
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
  //AliHLTUInt32_t rcuId = (Int_t)((word >> 7) & 0x1FF);
  AliHLTUInt32_t rcuTrailerSize32 = (word & 0x7F); // size of RCU trailer data in words
  
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
    HLTWarning(Form("%s Invalid trailer size found (%d bytes) ! The size is bigger than the raw data size (%ld bytes)!",
		    str, rcuTrailerSize32*4,dataSize32*4));
    return 0;    
  }

  // check the trailer

  Int_t trailerIndex = dataSize32 - rcuTrailerSize32;

  for( unsigned int i=trailerIndex; i<dataSize32-1; i++){
    if ((fData[i] >> 30) != 2) {
      HLTWarning("%s Missing RCU trailer identifier pattern!",str);
      continue;
    }
  }

  // Read the payload size  
 
  Int_t  rcuPayloadSize32 = fData[trailerIndex] & 0x3FFFFFF;

  if ( rcuPayloadSize32 + rcuTrailerSize32  != dataSize32) {
    HLTWarning(Form("%s Inconsistent raw data size ! Raw data size - %ld bytes (from CDH), RCU trailer - %d bytes, raw data size (from RCU trailer) - %d bytes !",
		    str, dataSize32*4,
		    (rcuTrailerSize32)*4,
		    rcuPayloadSize32*4));
    return 0;
  }
    
  
  //AliHLTTPCMapping *mapping = new AliHLTTPCMapping(patch);
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

    //bool channelBranchAB = ((channelHWAddress >> 11) & 0x1);
    // int channelFEC       = ((channelHWAddress >> 7) & 0xF); // front-end card index
    //int channelAltro = ((channelHWAddress >> 4) & 0x7); // altro chip index
    //int channelIndex = (channelHWAddress & 0xF); // channel index
    //int channelRow = mapping->GetRow(channelHWAddress);
    //int channelPad = mapping->GetPad(channelHWAddress);

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
	
      //UShort_t* bunchData = &channelData10[channelData10Index];   // pointer to the current bunch samples
      channelData10Index += bunchLength;            
    }
  }
  
  delete[] channelData10;
  //delete mapping;

  return 1;
}

