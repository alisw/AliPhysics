// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
 *          Jochen Thaeder <thaeder@kip.uni-heidelberg.de>                *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCDigitReaderRaw.cxx
    @author Timm Steinbeck
    @date   
    @brief  A digit reader implementation for the RAW data coming from the RCU.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCDigitReaderRaw.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCRootTypes.h"
#include "AliHLTStdIncludes.h"
#include "AliHLTTPCLogging.h"

ClassImp(AliHLTTPCDigitReaderRaw)

AliHLTTPCDigitReaderRaw::AliHLTTPCDigitReaderRaw( unsigned formatVersion )
  :
  fBuffer(NULL),
  fBufferSize(0),
  fPatch(-1),
  fSlice(-1),
  fDataFormatVersion(formatVersion),
  fCurrentRow(0),
  fCurrentPad(0),
  fCurrentBin(-1),
  fVerify(false),
  
  // For sorting
  fNRows(0),
  fRowOffset(0),
  fNMaxRows(0),
  fNMaxPads(0),
  fNTimeBins(0),
  fData(NULL),
  fMapErrThrown(0)
{
#ifndef HAVE_TPC_MAPPING
  memset(fMapping0, 0, fMapping0Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping1, 0, fMapping1Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping2, 0, fMapping2Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping3, 0, fMapping3Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping4, 0, fMapping4Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping5, 0, fMapping5Size*fMappingDimension*sizeof(Int_t));
#endif //#ifndef HAVE_TPC_MAPPING

    if ( fDataFormatVersion==0 || fDataFormatVersion==2 || fDataFormatVersion==4 )
      {
	
	// get max number of rows
	for (Int_t ii=0; ii < 6; ii++)
	  if (AliHLTTPCTransform::GetNRows(ii) > fNMaxRows) 
	    fNMaxRows = AliHLTTPCTransform::GetNRows(ii);
	
	// get max number of pads
	for (Int_t ii=0; ii < AliHLTTPCTransform::GetNRows();ii++ )
	  if (AliHLTTPCTransform::GetNPads(ii) > fNMaxPads) 
	    fNMaxPads = AliHLTTPCTransform::GetNPads(ii);
	
	// get max number of bins
	fNTimeBins = AliHLTTPCTransform::GetNTimeBins();
	
	//	HLTDebug("Array Borders ||| MAXPAD=%d ||| MAXROW=%d ||| MAXBIN=%d ||| MAXMUL=%d", 
	//	 fNMaxPads, fNMaxRows, fNTimeBins, fNTimeBins*fNMaxRows*fNMaxPads);
	
	// init Data array
	fData = new Int_t[ fNMaxRows*fNMaxPads*fNTimeBins ];
      }
}

AliHLTTPCDigitReaderRaw::AliHLTTPCDigitReaderRaw(const AliHLTTPCDigitReaderRaw& src)
  :
  fBuffer(NULL),
  fBufferSize(0),
  fPatch(-1),
  fSlice(-1),
  fDataFormatVersion(src.fDataFormatVersion),
  fCurrentRow(0),
  fCurrentPad(0),
  fCurrentBin(-1),
  fVerify(false),
  
  // For sorting
  fNRows(0),
  fRowOffset(0),
  fNMaxRows(0),
  fNMaxPads(0),
  fNTimeBins(0),
  fData(NULL),
  fMapErrThrown(0)
{
#ifndef HAVE_TPC_MAPPING
  memset(fMapping0, 0, fMapping0Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping1, 0, fMapping1Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping2, 0, fMapping2Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping3, 0, fMapping3Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping4, 0, fMapping4Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping5, 0, fMapping5Size*fMappingDimension*sizeof(Int_t));
#endif //#ifndef HAVE_TPC_MAPPING
  HLTFatal("copy constructor not for use");
}

AliHLTTPCDigitReaderRaw& AliHLTTPCDigitReaderRaw::operator=(const AliHLTTPCDigitReaderRaw& src)
{
#ifndef HAVE_TPC_MAPPING
  memset(fMapping0, 0, fMapping0Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping1, 0, fMapping1Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping2, 0, fMapping2Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping3, 0, fMapping3Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping4, 0, fMapping4Size*fMappingDimension*sizeof(Int_t));
  memset(fMapping5, 0, fMapping5Size*fMappingDimension*sizeof(Int_t));
#endif //#ifndef HAVE_TPC_MAPPING
  fBuffer=NULL;
  fBufferSize=0;
  fPatch=-1;
  fSlice=-1;
  fDataFormatVersion=src.fDataFormatVersion;
  fCurrentRow=0;
  fCurrentPad=0;
  fCurrentBin=-1;
  fVerify=false;
  
  // For sorting
  fNRows=0;
  fRowOffset=0;
  fNMaxRows=0;
  fNMaxPads=0;
  fNTimeBins=0;
  fData=NULL;
  fMapErrThrown=0;
  HLTFatal("assignment operator not for use");
  return (*this);
}

AliHLTTPCDigitReaderRaw::~AliHLTTPCDigitReaderRaw(){
  if ( fDataFormatVersion==0 || fDataFormatVersion==2 || fDataFormatVersion==4 )
    {
      if ( fData )
	delete [] fData;
      fData = NULL;
    }
}

int AliHLTTPCDigitReaderRaw::InitBlock(void* ptr,unsigned long size,Int_t firstrow,Int_t lastrow, Int_t patch, Int_t slice) {
  return AliHLTTPCDigitReader::InitBlock(ptr, size, firstrow, lastrow, patch, slice);
}

int AliHLTTPCDigitReaderRaw::InitBlock(void* ptr,unsigned long size, Int_t patch, Int_t slice){

    fBuffer = (AliHLTUInt8_t*) ptr;
    if (fBuffer==NULL) {
      HLTError("invalid data buffer");
      return -EINVAL;
    }
    fBufferSize = size;
    if (fBufferSize<=0) HLTWarning("no data available: zero length buffer");
    fPatch = patch;
    fSlice = slice;
    fPad = -1;
    fRow = -1;
    
    fAltroBlockPositionBytes = 0;
    fAltroBlockLengthBytes = 0;
    fAltroBlock10BitWordCnt = 0xFFFFU;
    fAltroBlockHWAddress = 0xFFFFU;
    fBunchPosition = 0xFFFFU;
    fBunchTimebinStart = ~0U;
    fBunchLength = 0;
    fWordInBunch = (unsigned)-1;

    Int_t firstrow=AliHLTTPCTransform::GetFirstRow(patch);
    Int_t lastrow=AliHLTTPCTransform::GetLastRow(patch);

    if ( fDataFormatVersion==0 || fDataFormatVersion==2 || fDataFormatVersion==4 )
      {
	fCurrentRow = 0;
	fCurrentPad = 0;
	fCurrentBin = -1;
	
	fNRows = lastrow - firstrow + 1;
	
	Int_t offset=0;
	if (patch > 1) offset =  AliHLTTPCTransform::GetFirstRow( 2 );
	
	fRowOffset = firstrow - offset;
	firstrow -= offset;
	lastrow  -= offset;
	
	// Init array with -1
	memset( fData, 0xFF, sizeof(Int_t)*(fNMaxRows*fNMaxPads*fNTimeBins) );

	const Int_t maxErrorPrintout=20;
	Int_t errorCount=0;
	Int_t entryCount=0;
	// read data and fill in array

	while( RealNext()){

	  entryCount++;
	  Int_t row = GetRealRow();
	  Int_t pad = GetRealPad();
	  Int_t bin = GetRealTime();
	  
//	  HLTFatal("Index out of array range: PAD=%d ||| ROW=%d ||| BIN=%d ||| OFFSET=%d ||| ROWOFFSET=%d", pad, row, bin, offset, fRowOffset);

	  if ( row < firstrow || row > lastrow || pad > AliHLTTPCTransform::GetNPads(row + offset) || bin > fNTimeBins || pad<0 || bin<0){
//	  if ( row < firstrow || row > lastrow || pad > AliHLTTPCTransform::GetNPads(row + offset) || bin > fNTimeBins){
	    if (errorCount++<maxErrorPrintout) {
	      HLTFatal("Index out of range. Probably wrong patch! slice %d - patch %d", slice, patch);
	      HLTFatal("PAD=%d out of %d ||| ROW=%d (%d to %d)  ||| BIN=%d out of %d  ||| OFFSET=%d ||| ROWOFFSET=%d",
		       pad, AliHLTTPCTransform::GetNPads(row + offset), row, firstrow, lastrow, bin, fNTimeBins,
		       offset, fRowOffset);

	      if ( row < firstrow || row > lastrow ) 
		HLTFatal("Row out of range: %d  ( %d to %d)", row, firstrow, lastrow);
	      if ( pad > AliHLTTPCTransform::GetNPads(row + offset) ) 
		HLTFatal("Pad out of range: %d  (pad count %d)", pad, AliHLTTPCTransform::GetNPads(row + offset));
	      if ( bin > fNTimeBins )
		HLTFatal("Time bin out of range: %d (bin count %d)", bin, fNTimeBins);
	    }
	    // stop at the fist error message in order to avoid endless messages and
	    // to handle corrupted events
	    //continue;
	    break;
	  } else if ((row-fRowOffset)*fNMaxPads*fNTimeBins+ pad*fNTimeBins + bin >=  fNMaxRows*fNMaxPads*fNTimeBins ) {
	    if (errorCount++<maxErrorPrintout) {
	      HLTFatal("index out of range: PAD=%d ||| ROW=%d ||| BIN=%d ||| OFFSET=%d ||| ROWOFFSET=%d", pad, row, bin, offset, fRowOffset);
	    }
	    // stop at the fist error message in order to avoid endless messages and
	    // to handle corrupted events
	    //continue;
	    break;
	  } else {
	    fData[ (row-fRowOffset)*fNMaxPads*fNTimeBins+ pad*fNTimeBins + bin ] = GetRealSignal() ;
	  }
	}
	if (errorCount>0) {
	  HLTFatal("%d of %d entries out of range", errorCount, entryCount);
	}
      }

    return 0;
}

bool AliHLTTPCDigitReaderRaw::Next(){
  if ( fDataFormatVersion==0 || fDataFormatVersion==2 || fDataFormatVersion==4 )
    {
      Bool_t readvalue = kTRUE;
      while (1) {
	fCurrentBin++;
	if (fCurrentBin >= fNTimeBins){
	  fCurrentBin = 0;
	  fCurrentPad++;
	  
	  if (fCurrentPad >=fNMaxPads){
	    fCurrentPad = 0;
	    fCurrentRow++;
	    
	    if (fCurrentRow >= fNMaxRows){
	      readvalue = kFALSE;
	      break;
	    }
	  }
	}
	
	if (fCurrentRow*fNMaxPads*fNTimeBins+ fCurrentPad*fNTimeBins + fCurrentBin >=  fNMaxRows*fNMaxPads*fNTimeBins ) {
	  HLTFatal("Overflow: fCurrentRow=%d fCurrentPad=%d fCurrentBin=%d", fCurrentRow, fCurrentPad, fCurrentBin);
	  readvalue = kFALSE;
	  break;
	}
	
	if (fData[ fCurrentRow*fNMaxPads*fNTimeBins + fCurrentPad*fNTimeBins + fCurrentBin  ] != -1) break;
      }
      return readvalue;
    }
  else
    return RealNext();
}

int AliHLTTPCDigitReaderRaw::GetRow(){
  if ( fDataFormatVersion==0 || fDataFormatVersion==2 || fDataFormatVersion==4 )
    {
      return (fCurrentRow + fRowOffset);
    }
  else
    return GetRealRow();
}
int AliHLTTPCDigitReaderRaw::GetPad(){
  if ( fDataFormatVersion==0 || fDataFormatVersion==2 || fDataFormatVersion==4 )
    {
      return fCurrentPad;
    }
  else
    return GetRealPad();
}
int AliHLTTPCDigitReaderRaw::GetSignal(){
  if ( fDataFormatVersion==0 || fDataFormatVersion==2 || fDataFormatVersion==4 )
    {
      return fData[ fCurrentRow*fNMaxPads*fNTimeBins+ fCurrentPad*fNTimeBins + fCurrentBin ];
    }
  else
    return GetRealSignal();
}
int AliHLTTPCDigitReaderRaw::GetTime(){
  if ( fDataFormatVersion==0 || fDataFormatVersion==2 || fDataFormatVersion==4 )
    {
      return fCurrentBin;
    }
  else
    return GetRealTime();
}


bool AliHLTTPCDigitReaderRaw::RealNext(){
//    printf( "%u %u %u %u %u\n", fBunchPosition, fBunchLength, fBunchTimebinStart, fWordInBunch, (unsigned)fAltroBlock10BitWordCnt );
    fWordInBunch++; // use next word in bunch
    if ( fWordInBunch==fBunchLength ) { // we have a bunch at all but have reached its end (or do not have an altro block yet)
	if ( fBunchPosition+fBunchLength==fAltroBlock10BitWordCnt ) { // We were at the last bunch of this altro block (or do not have an altro block yet)
	    if ( !NextAltroBlock() )
		return false;
	    fBunchPosition = 0;
	}
	else {
	    fBunchPosition += fBunchLength;
	}
	fBunchLength = GetAltroBlock10BitWord( fBunchPosition );
	fBunchTimebinStart = GetAltroBlock10BitWord( fBunchPosition+1 );
	fWordInBunch = 2;
    }
    //HLTDebug( "%u %u %u %u %u\n", fBunchPosition, fBunchLength, fBunchTimebinStart, fWordInBunch, (unsigned)fAltroBlock10BitWordCnt );
    return true;
}
int AliHLTTPCDigitReaderRaw::GetRealRow(){
    return fRow;
}
int AliHLTTPCDigitReaderRaw::GetRealPad(){
    return fPad;
}
int AliHLTTPCDigitReaderRaw::GetRealSignal(){
    return GetAltroBlock10BitWord( fBunchPosition+fWordInBunch );
}
int AliHLTTPCDigitReaderRaw::GetRealTime(){
  //HLTDebug( "GetRealTime: %u - %u\n", fBunchTimebinStart, fWordInBunch );
    return fBunchTimebinStart-(fWordInBunch-2);
}

AliHLTUInt32_t AliHLTTPCDigitReaderRaw::GetRCUTrailer( unsigned offset ){
  if (fBufferSize<=0) return 0;
  unsigned rcuDataBlockLen = GetRCUDataBlockLength(); 
  if ( offset >= rcuDataBlockLen ) return 0;
  return ((AliHLTUInt32_t*)(fBuffer+fBufferSize-rcuDataBlockLen))[offset];
}

bool AliHLTTPCDigitReaderRaw::NextAltroBlock()
    {
    if (fBufferSize<=0) return 0;
    bool first = false;
    if ( !fAltroBlockLengthBytes )
	{
	// First block in back linked list (last block in memory)
	fAltroBlockPositionBytes = fBufferSize-GetRCUDataBlockLength();
	first = true;
	}
    else
	{
	if ( fAltroBlockPositionBytes<fAltroBlockLengthBytes+GetCommonDataHeaderSize() )
	  {
	    HLTFatal("Inconsistent Data: fAltroBlockPositionBytes=%d fAltroBlockLengthBytes=%d", fAltroBlockPositionBytes, fAltroBlockLengthBytes);
	  }
	if ( fAltroBlockPositionBytes<=fAltroBlockLengthBytes+GetCommonDataHeaderSize() )
	    return false; // We have reached the end of the back linked list
	fAltroBlockPositionBytes -= fAltroBlockLengthBytes;
	}

      AliHLTUInt64_t altroTrailerWord = GetAltroBlock40BitWord( 0 );
      // Undefined hack from experience to match fill words appearing in simulated data
      // Seem to be between 0 and 3 fill words, most likely to bring the number of 40bit words
      // to a multiple of four / to bring the total number of bytes to a common multiple of 4 and 5.
      // (RCU sends 40 bit (5 byte) words, DDL uses 32 bit (4 bytes) words.
      unsigned short tmpCnt=0;
      //HLTDebug( "Altro trailer word 0: 0x%016LX\n", altroTrailerWord );
      while ( first && altroTrailerWord==0x000000AAAAAAAAAAULL && tmpCnt++<4 ) // Allow up to 4 fill values
	{
	  altroTrailerWord = GetAltroBlock40BitWord( tmpCnt );
	  //HLTDebug( "Altro trailer word %hu: 0x%016LX\n", tmpCnt, altroTrailerWord );
	}

      fAltroBlockPositionBytes -= 5*tmpCnt;
      if ( fVerify && ((altroTrailerWord & 0xFFFC000000ULL)!=0xAAA8000000ULL) )
	{
	  HLTFatal("Data inconsistency in Altro Block at byte position %#x (%d): Expected 0x2AAA in high 14 bits of altro trailer word; Found %#llx (%#llx)",
		   fAltroBlockPositionBytes, fAltroBlockPositionBytes, 
		   ((altroTrailerWord & 0xFFFC000000ULL) >> 26), altroTrailerWord);


	  return false;
	}

      if ( fVerify && ((altroTrailerWord & 0x000000F000ULL)!=0x000000A000ULL) )
	{
	  HLTFatal("Data inconsistency in Altro Block at byte position %#x (%d): Expected 0xA in bits 12-15 of altro trailer word; Found %#llx .",
		   fAltroBlockPositionBytes, fAltroBlockPositionBytes,  ((altroTrailerWord & 0x000000F000ULL) >> 12)); 

	  return false;
	}

      fAltroBlock10BitWordCnt = (altroTrailerWord >> 16) & 0x3FF;
      fAltroBlockHWAddress = altroTrailerWord & 0xFFF;

      // ApplyMapping
      if (!ApplyMapping())
	{
	  HLTFatal("Mapping failed Patch %d HWA %#x (%d) - maxHWA %#x (%d)",
		   fPatch, fAltroBlockHWAddress, fAltroBlockHWAddress, fMaxHWA[fPatch], fMaxHWA[fPatch]);

	}

      unsigned words40Bit = fAltroBlock10BitWordCnt/4;
      if ( fAltroBlock10BitWordCnt % 4 )
	  words40Bit++;
      words40Bit++;
      fAltroBlockLengthBytes = words40Bit*5;
    if ( fAltroBlock10BitWordCnt % 4 )
	fAltroBlock10BitFillWordCnt = 4-(fAltroBlock10BitWordCnt % 4);
    else
	fAltroBlock10BitFillWordCnt=0;
    if ( fVerify )
      {
	for ( unsigned b = 0; b < fAltroBlock10BitFillWordCnt; b++ )
	  {
	    if ( GetAltroBlockReal10BitWord(b)!=0x2AA )
	      {
		HLTFatal("Data inconsistency in trailing 10 bit fill word of Altro Block at byte position %#x (%d): Expected 0x2AA; Found %#x",
			 fAltroBlockPositionBytes, fAltroBlockPositionBytes, GetAltroBlockReal10BitWord(b));
		
		return false;
	      }
	  }
      }
    return true;
    }

AliHLTUInt32_t AliHLTTPCDigitReaderRaw::GetAltroBlockHWaddr(){
return fAltroBlockHWAddress;
}
unsigned AliHLTTPCDigitReaderRaw::GetAltroBlock10BitWordCnt(){
return fAltroBlock10BitWordCnt;
}
AliHLTUInt64_t AliHLTTPCDigitReaderRaw::GetAltroBlock40BitWord( unsigned long ndx ){
AliHLTUInt64_t val=0;
unsigned wordOffset32Bit = (ndx / 4)*5;
switch ( ndx % 4 ) // 40 bit word index in a 4*40 bit=5*32 bit group
    {
    case 0:
	val = (*(AliHLTUInt32_t*)(fBuffer+fAltroBlockPositionBytes-(wordOffset32Bit+1)*sizeof(AliHLTUInt32_t)));
	val <<= 8;
	val |= (*(AliHLTUInt32_t*)(fBuffer+fAltroBlockPositionBytes-(wordOffset32Bit+2)*sizeof(AliHLTUInt32_t))) >> 24;
	break;
    case 1:
	val = ((*(AliHLTUInt32_t*)(fBuffer+fAltroBlockPositionBytes-(wordOffset32Bit+2)*sizeof(AliHLTUInt32_t))) & 0x00FFFFFF);
	val <<= 16;
	val |= ((*(AliHLTUInt32_t*)(fBuffer+fAltroBlockPositionBytes-(wordOffset32Bit+3)*sizeof(AliHLTUInt32_t))) >> 16) & 0xFFFF;
	break;
    case 2:
	val = ((*(AliHLTUInt32_t*)(fBuffer+fAltroBlockPositionBytes-(wordOffset32Bit+3)*sizeof(AliHLTUInt32_t))) & 0xFFFF);
	val <<= 24;
	val |= ((*(AliHLTUInt32_t*)(fBuffer+fAltroBlockPositionBytes-(wordOffset32Bit+4)*sizeof(AliHLTUInt32_t))) >> 8);
	break;
    case 3:
	val = ((*(AliHLTUInt32_t*)(fBuffer+fAltroBlockPositionBytes-(wordOffset32Bit+4)*sizeof(AliHLTUInt32_t))) & 0xFF);
	val <<= 32;
	val |= *(AliHLTUInt32_t*)(fBuffer+fAltroBlockPositionBytes-(wordOffset32Bit+5)*sizeof(AliHLTUInt32_t));
	break;
    }
return val;
}
AliHLTUInt16_t AliHLTTPCDigitReaderRaw::GetAltroBlock10BitWord( unsigned long ndx ){
unsigned long realNdx = ndx+fAltroBlock10BitFillWordCnt;
unsigned long word40BitNdx = (realNdx / 4)+1;
AliHLTUInt64_t word40Bit = GetAltroBlock40BitWord( word40BitNdx );
switch ( realNdx % 4 )
    {
    case 3:
	return word40Bit & 0x3FF;
    case 2:
	return (word40Bit>>10) & 0x3FF;
    case 1:
	return (word40Bit>>20) & 0x3FF;
    case 0:
	return (word40Bit>>30) & 0x3FF;
    }

 return 0xFFFF; 
}

AliHLTUInt16_t AliHLTTPCDigitReaderRaw::GetAltroBlockReal10BitWord( unsigned long ndx ){
unsigned long word40BitNdx = (ndx / 4)+1;
AliHLTUInt64_t word40Bit = GetAltroBlock40BitWord( word40BitNdx );
switch ( ndx % 4 )
    {
    case 3:
	return word40Bit & 0x3FF;
    case 2:
	return (word40Bit>>10) & 0x3FF;
    case 1:
	return (word40Bit>>20) & 0x3FF;
    case 0:
	return (word40Bit>>30) & 0x3FF;
    }

 return 0xFFFF; 
}

// Return length of trailing RCU data block in bytes
unsigned AliHLTTPCDigitReaderRaw::GetRCUDataBlockLength() const
    {
    switch ( fDataFormatVersion )
	{
	case 0:
	case 1:
	    return 4;
	    break;
	case 2:
	case 3:
	    return 12;
	    break;
	case 4:
	case 5:
	    return 8;
	    break;
	default:
	    return fBufferSize;
	}
    }

unsigned AliHLTTPCDigitReaderRaw::GetCommonDataHeaderSize() const
    {
    return 32;
    }


Bool_t AliHLTTPCDigitReaderRaw::ApplyMapping(){

#ifndef HAVE_TPC_MAPPING
  if (fMapErrThrown++==0) {
    HLTFatal("mapping not available, you must compile with HAVE_TPC_MAPPING");
  }
  return -1;
#endif //#ifndef HAVE_TPC_MAPPING
    if ( (unsigned)fAltroBlockHWAddress > fMaxHWA[fPatch]){
	fPad = -1;
	fRow = -1;
	return kFALSE;
    }

    switch(fPatch){
	case 0:
	    fRow = fMapping0[(unsigned)fAltroBlockHWAddress][0];
	    fPad = fMapping0[(unsigned)fAltroBlockHWAddress][1];
	    break;
        case 1:
	    fRow = AliHLTTPCDigitReaderRaw::fMapping1[(unsigned)fAltroBlockHWAddress][0];
	    fPad = AliHLTTPCDigitReaderRaw::fMapping1[(unsigned)fAltroBlockHWAddress][1];
#if 0
	    printf ("pad %d # row %d (hwa: %u / 0x%08X\n", fMapping1[(unsigned)fAltroBlockHWAddress][0],fMapping1[(unsigned)fAltroBlockHWAddress][1], (unsigned)fAltroBlockHWAddress, (unsigned)fAltroBlockHWAddress);
	    printf ("pad %d # row %d (hwa: %u / 0x%08X\n", fMapping1[(unsigned)fAltroBlockHWAddress-1][0],fMapping1[(unsigned)fAltroBlockHWAddress-1][1], (unsigned)fAltroBlockHWAddress-1, (unsigned)fAltroBlockHWAddress-1);
	    printf ("pad %d # row %d (hwa: %u / 0x%08X\n", fMapping1[(unsigned)fAltroBlockHWAddress+1][0],fMapping1[(unsigned)fAltroBlockHWAddress+1][1], (unsigned)fAltroBlockHWAddress+1, (unsigned)fAltroBlockHWAddress+1);
#endif
	    break;
	case 2:
	    fRow = fMapping2[(unsigned)fAltroBlockHWAddress][0];
	    fPad = fMapping2[(unsigned)fAltroBlockHWAddress][1];
	    break;
        case 3:
	    fRow = fMapping3[(unsigned)fAltroBlockHWAddress][0];
	    fPad = fMapping3[(unsigned)fAltroBlockHWAddress][1];
	    break;
	case 4:
	    fRow = fMapping4[(unsigned)fAltroBlockHWAddress][0];
	    fPad = fMapping4[(unsigned)fAltroBlockHWAddress][1];
	    break;
        case 5:
	    fRow = fMapping5[(unsigned)fAltroBlockHWAddress][0];
	    fPad = fMapping5[(unsigned)fAltroBlockHWAddress][1];
	    break;
	default:
	    fRow = -1;
	    fPad = -1;
	    return kFALSE;
    }
    return kTRUE;
}


Int_t AliHLTTPCDigitReaderRaw::GetRow( unsigned patch, unsigned hwAddr )
{
#ifndef HAVE_TPC_MAPPING
  if (fMapErrThrown++==0) {
    HLTFatal("mapping not available, you must compile with HAVE_TPC_MAPPING");
  }
  return -1;
#endif //#ifndef HAVE_TPC_MAPPING
    if ( (unsigned)hwAddr > fMaxHWA[fPatch]){
	return -1;
    }

    switch(fPatch){
	case 0:
	    return fMapping0[hwAddr][0];
        case 1:
	    return fMapping1[hwAddr][0];
	case 2:
	    return fMapping2[hwAddr][0];
        case 3:
	    return fMapping3[hwAddr][0];
	case 4:
	    return fMapping4[hwAddr][0];
        case 5:
	    return fMapping5[hwAddr][0];
	default:
	  return -1;
    }
}
Int_t AliHLTTPCDigitReaderRaw::GetPad( unsigned patch, unsigned hwAddr )
{
#ifndef HAVE_TPC_MAPPING
  if (fMapErrThrown++==0) {
    HLTFatal("mapping not available, you must compile with HAVE_TPC_MAPPING");
  }
  return -1;
#endif //#ifndef HAVE_TPC_MAPPING
    if ( (unsigned)hwAddr > fMaxHWA[fPatch]){
	return -1;
    }

    switch(fPatch){
	case 0:
	    return fMapping0[hwAddr][1];
        case 1:
	    return fMapping1[hwAddr][1];
	case 2:
	    return fMapping2[hwAddr][1];
        case 3:
	    return fMapping3[hwAddr][1];
	case 4:
	    return fMapping4[hwAddr][1];
        case 5:
	    return fMapping5[hwAddr][1];
	default:
	  return -1;
    }
}

unsigned AliHLTTPCDigitReaderRaw::GetMaxHWA( unsigned patch )
{
  if ( patch>=fNofPatches )
    return 0;
  return fMaxHWA[patch];
}

Int_t AliHLTTPCDigitReaderRaw::DecodeMode(Int_t mode) {
  Int_t decodedMode;

  if ( mode >= kNofRawReaderModes ) 
    decodedMode = -1;
  else
    decodedMode = mode;

  return decodedMode;
}

Int_t AliHLTTPCDigitReaderRaw::DecodeMode(const Char_t *mode) {
  Int_t decodedMode;
  Char_t *cpErr;

  // Check if String is convertible to Int_t
  // if not decode the string, otherwise, check if Int_t is valid
  Int_t intMode = strtoul( mode, &cpErr ,0);

  if ( *cpErr ) {
    if ( !strcmp( mode, "sorted_3_trailerword" ) ) 
      decodedMode = kSorted3Trailerword;
    
    else if ( !strcmp( mode, "sorted_2_trailerword" ) ) 
      decodedMode = kSorted2Trailerword;
    
    else if ( !strcmp( mode, "sorted_1_trailerword" ) ) 
      decodedMode = kSorted1Trailerword;
    
    else if ( !strcmp( mode, "unsorted_3_trailerword" ) ) 
      decodedMode = kUnsorted3Trailerword;
    
    else if ( !strcmp( mode, "unsorted_2_trailerword" ) ) 
      decodedMode = kUnsorted2Trailerword;
    
    else if ( !strcmp( mode, "unsorted_1_trailerword" ) ) 
      decodedMode = kUnsorted1Trailerword;
    
    else if ( ! strcmp( mode, "offline" ) )
      decodedMode = -2;
    
    else 
      decodedMode = -1;
  }  // END if ( *cpErr ) {
  else {
    if ( intMode >= kNofRawReaderModes ) 
      decodedMode = -1;
    else
      decodedMode = intMode;
  }

  return decodedMode;
}


// ----- MAPPING ARRAYS
#if defined(HAVE_TPC_MAPPING)
#include "mapping_array_out.inc"
#else
// dummy definitions in case of missing mapping
Int_t AliHLTTPCDigitReaderRaw::fMapping0[fMapping0Size][fMappingDimension];
Int_t AliHLTTPCDigitReaderRaw::fMapping1[fMapping1Size][fMappingDimension];
Int_t AliHLTTPCDigitReaderRaw::fMapping2[fMapping2Size][fMappingDimension];
Int_t AliHLTTPCDigitReaderRaw::fMapping3[fMapping3Size][fMappingDimension];
Int_t AliHLTTPCDigitReaderRaw::fMapping4[fMapping4Size][fMappingDimension];
Int_t AliHLTTPCDigitReaderRaw::fMapping5[fMapping5Size][fMappingDimension];
#endif //#if defined(HAVE_TPC_MAPPING)
