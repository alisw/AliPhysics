// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
 *                  Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
 *                  for The ALICE HLT Project.                            *
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
  fRow(-1),
  fPad(-1),
  fAltroBlockPositionBytes(0),
  fAltroBlockLengthBytes(0),
  fAltroBlockHWAddress(0),
  fAltroBlock10BitWordCnt(0),
  fAltroBlock10BitFillWordCnt(0),
  fDataFormatVersion(formatVersion),
  fBunchPosition(0xFFFFU),
  fBunchTimebinStart(~0U),
  fBunchLength(0),
  fWordInBunch((unsigned)-1),
  fVerify(false),
  
  fCurrentRow(0),
  fCurrentPad(0),
  fCurrentBin(-1),
  fRowOffset(0),
  fNRows(0),
  fNMaxRows(0),
  fNMaxPads(0),
  fNTimeBins(0),
  fData(NULL),
  fMapErrThrown(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
#ifndef HAVE_TPC_MAPPING
  memset(fgMapping0, 0, fgkMapping0Size*fgkMappingDimension*sizeof(Int_t));
  memset(fgMapping1, 0, fgkMapping1Size*fgkMappingDimension*sizeof(Int_t));
  memset(fgMapping2, 0, fgkMapping2Size*fgkMappingDimension*sizeof(Int_t));
  memset(fgMapping3, 0, fgkMapping3Size*fgkMappingDimension*sizeof(Int_t));
  memset(fgMapping4, 0, fgkMapping4Size*fgkMappingDimension*sizeof(Int_t));
  memset(fgMapping5, 0, fgkMapping5Size*fgkMappingDimension*sizeof(Int_t));
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

AliHLTTPCDigitReaderRaw::~AliHLTTPCDigitReaderRaw()
{
  // see header file for class documentation
  if ( fDataFormatVersion==0 || fDataFormatVersion==2 || fDataFormatVersion==4 )
    {
      if ( fData )
	delete [] fData;
      fData = NULL;
    }
}

int AliHLTTPCDigitReaderRaw::InitBlock(void* ptr,unsigned long size,Int_t firstrow,Int_t lastrow, Int_t patch, Int_t slice) 
{
  // see header file for class documentation
  return AliHLTTPCDigitReader::InitBlock(ptr, size, firstrow, lastrow, patch, slice);
}

int AliHLTTPCDigitReaderRaw::InitBlock(void* ptr,unsigned long size, Int_t patch, Int_t slice)
{
  // see header file for class documentation

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

	const Int_t kMaxErrorPrintout=20;
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
	    if (errorCount++<kMaxErrorPrintout) {
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
	    if (errorCount++<kMaxErrorPrintout) {
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

bool AliHLTTPCDigitReaderRaw::Next()
{
  // see header file for class documentation

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

int AliHLTTPCDigitReaderRaw::GetRow()
{
  // see header file for class documentation

  if ( fDataFormatVersion==0 || fDataFormatVersion==2 || fDataFormatVersion==4 )
    {
      return (fCurrentRow + fRowOffset);
    }
  else
    return GetRealRow();
}

int AliHLTTPCDigitReaderRaw::GetPad()
{
  // see header file for class documentation

  if ( fDataFormatVersion==0 || fDataFormatVersion==2 || fDataFormatVersion==4 )
    {
      return fCurrentPad;
    }
  else
    return GetRealPad();
}

int AliHLTTPCDigitReaderRaw::GetSignal()
{
  // see header file for class documentation

  if ( fDataFormatVersion==0 || fDataFormatVersion==2 || fDataFormatVersion==4 )
    {
      return fData[ fCurrentRow*fNMaxPads*fNTimeBins+ fCurrentPad*fNTimeBins + fCurrentBin ];
    }
  else
    return GetRealSignal();
}

int AliHLTTPCDigitReaderRaw::GetTime()
{
  // see header file for class documentation

  if ( fDataFormatVersion==0 || fDataFormatVersion==2 || fDataFormatVersion==4 )
    {
      return fCurrentBin;
    }
  else
    return GetRealTime();
}

bool AliHLTTPCDigitReaderRaw::RealNext()
{
  // see header file for class documentation

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

int AliHLTTPCDigitReaderRaw::GetRealRow() const
{
  // see header file for class documentation
    return fRow;
}

int AliHLTTPCDigitReaderRaw::GetRealPad() const
{
  // see header file for class documentation
    return fPad;
}

int AliHLTTPCDigitReaderRaw::GetRealSignal()
{
  // see header file for class documentation
    return GetAltroBlock10BitWord( fBunchPosition+fWordInBunch );
}

int AliHLTTPCDigitReaderRaw::GetRealTime() const
{
  // see header file for class documentation
  //HLTDebug( "GetRealTime: %u - %u\n", fBunchTimebinStart, fWordInBunch );
    return fBunchTimebinStart-(fWordInBunch-2);
}

AliHLTUInt32_t AliHLTTPCDigitReaderRaw::GetRCUTrailer( unsigned offset ) const
{
  // see header file for class documentation
  if (fBufferSize<=0) return 0;
  unsigned rcuDataBlockLen = GetRCUDataBlockLength(); 
  if ( offset >= rcuDataBlockLen ) return 0;
  return ((AliHLTUInt32_t*)(fBuffer+fBufferSize-rcuDataBlockLen))[offset];
}

bool AliHLTTPCDigitReaderRaw::NextAltroBlock()
{
  // see header file for class documentation
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
		   fPatch, fAltroBlockHWAddress, fAltroBlockHWAddress, fgMaxHWA[fPatch], fgMaxHWA[fPatch]);

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

AliHLTUInt32_t AliHLTTPCDigitReaderRaw::GetAltroBlockHWaddr() const
{
  // see header file for class documentation
  return fAltroBlockHWAddress;
}

unsigned AliHLTTPCDigitReaderRaw::GetAltroBlock10BitWordCnt() const
{
  // see header file for class documentation
  return fAltroBlock10BitWordCnt;
}

AliHLTUInt64_t AliHLTTPCDigitReaderRaw::GetAltroBlock40BitWord( unsigned long ndx ) const
{
  // see header file for class documentation
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

AliHLTUInt16_t AliHLTTPCDigitReaderRaw::GetAltroBlock10BitWord( unsigned long ndx )
{
  // see header file for class documentation
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

AliHLTUInt16_t AliHLTTPCDigitReaderRaw::GetAltroBlockReal10BitWord( unsigned long ndx )
{
  // see header file for class documentation
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

unsigned AliHLTTPCDigitReaderRaw::GetRCUDataBlockLength() const
{
  // see header file for class documentation
  // Return length of trailing RCU data block in bytes
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
  // see header file for class documentation
  return 32;
}


Bool_t AliHLTTPCDigitReaderRaw::ApplyMapping()
{
  // see header file for class documentation

#ifndef HAVE_TPC_MAPPING
  if (fMapErrThrown++==0) {
    HLTFatal("mapping not available, you must compile with HAVE_TPC_MAPPING");
  }
  return -1;
#endif //#ifndef HAVE_TPC_MAPPING
    if ( (unsigned)fAltroBlockHWAddress > fgMaxHWA[fPatch]){
	fPad = -1;
	fRow = -1;
	return kFALSE;
    }

    switch(fPatch){
	case 0:
	    fRow = fgMapping0[(unsigned)fAltroBlockHWAddress][0];
	    fPad = fgMapping0[(unsigned)fAltroBlockHWAddress][1];
	    break;
        case 1:
	    fRow = AliHLTTPCDigitReaderRaw::fgMapping1[(unsigned)fAltroBlockHWAddress][0];
	    fPad = AliHLTTPCDigitReaderRaw::fgMapping1[(unsigned)fAltroBlockHWAddress][1];
#if 0
	    printf ("pad %d # row %d (hwa: %u / 0x%08X\n", fgMapping1[(unsigned)fAltroBlockHWAddress][0],fgMapping1[(unsigned)fAltroBlockHWAddress][1], (unsigned)fAltroBlockHWAddress, (unsigned)fAltroBlockHWAddress);
	    printf ("pad %d # row %d (hwa: %u / 0x%08X\n", fgMapping1[(unsigned)fAltroBlockHWAddress-1][0],fgMapping1[(unsigned)fAltroBlockHWAddress-1][1], (unsigned)fAltroBlockHWAddress-1, (unsigned)fAltroBlockHWAddress-1);
	    printf ("pad %d # row %d (hwa: %u / 0x%08X\n", fgMapping1[(unsigned)fAltroBlockHWAddress+1][0],fgMapping1[(unsigned)fAltroBlockHWAddress+1][1], (unsigned)fAltroBlockHWAddress+1, (unsigned)fAltroBlockHWAddress+1);
#endif
	    break;
	case 2:
	    fRow = fgMapping2[(unsigned)fAltroBlockHWAddress][0];
	    fPad = fgMapping2[(unsigned)fAltroBlockHWAddress][1];
	    break;
        case 3:
	    fRow = fgMapping3[(unsigned)fAltroBlockHWAddress][0];
	    fPad = fgMapping3[(unsigned)fAltroBlockHWAddress][1];
	    break;
	case 4:
	    fRow = fgMapping4[(unsigned)fAltroBlockHWAddress][0];
	    fPad = fgMapping4[(unsigned)fAltroBlockHWAddress][1];
	    break;
        case 5:
	    fRow = fgMapping5[(unsigned)fAltroBlockHWAddress][0];
	    fPad = fgMapping5[(unsigned)fAltroBlockHWAddress][1];
	    break;
	default:
	    fRow = -1;
	    fPad = -1;
	    return kFALSE;
    }
    return kTRUE;
}


Int_t AliHLTTPCDigitReaderRaw::GetRow( unsigned /*patch*/, unsigned hwAddr )
{
  // see header file for class documentation

#ifndef HAVE_TPC_MAPPING
  if (fMapErrThrown++==0) {
    HLTFatal("mapping not available, you must compile with HAVE_TPC_MAPPING");
  }
  return -1;
#endif //#ifndef HAVE_TPC_MAPPING
    if ( (unsigned)hwAddr > fgMaxHWA[fPatch]){
	return -1;
    }

    switch(fPatch){
	case 0:
	    return fgMapping0[hwAddr][0];
        case 1:
	    return fgMapping1[hwAddr][0];
	case 2:
	    return fgMapping2[hwAddr][0];
        case 3:
	    return fgMapping3[hwAddr][0];
	case 4:
	    return fgMapping4[hwAddr][0];
        case 5:
	    return fgMapping5[hwAddr][0];
	default:
	  return -1;
    }
}

Int_t AliHLTTPCDigitReaderRaw::GetPad( unsigned /*patch*/, unsigned hwAddr )
{
  // see header file for class documentation

#ifndef HAVE_TPC_MAPPING
  if (fMapErrThrown++==0) {
    HLTFatal("mapping not available, you must compile with HAVE_TPC_MAPPING");
  }
  return -1;
#endif //#ifndef HAVE_TPC_MAPPING
    if ( (unsigned)hwAddr > fgMaxHWA[fPatch]){
	return -1;
    }

    switch(fPatch){
	case 0:
	    return fgMapping0[hwAddr][1];
        case 1:
	    return fgMapping1[hwAddr][1];
	case 2:
	    return fgMapping2[hwAddr][1];
        case 3:
	    return fgMapping3[hwAddr][1];
	case 4:
	    return fgMapping4[hwAddr][1];
        case 5:
	    return fgMapping5[hwAddr][1];
	default:
	  return -1;
    }
}

unsigned AliHLTTPCDigitReaderRaw::GetMaxHWA( unsigned patch ) const
{
  // see header file for class documentation

  if ( (int)patch>=fgkNofPatches )
    return 0;
  return fgMaxHWA[patch];
}

Int_t AliHLTTPCDigitReaderRaw::DecodeMode(Int_t mode) 
{
  // see header file for class documentation

  Int_t decodedMode;

  if ( mode >= kNofRawReaderModes ) 
    decodedMode = -1;
  else
    decodedMode = mode;

  return decodedMode;
}

Int_t AliHLTTPCDigitReaderRaw::DecodeMode(const Char_t *mode) 
{
  // see header file for class documentation

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
unsigned AliHLTTPCDigitReaderRaw::fgMaxHWA[fgkNofPatches];
Int_t AliHLTTPCDigitReaderRaw::fgMapping0[fgkMapping0Size][fgkMappingDimension];
Int_t AliHLTTPCDigitReaderRaw::fgMapping1[fgkMapping1Size][fgkMappingDimension];
Int_t AliHLTTPCDigitReaderRaw::fgMapping2[fgkMapping2Size][fgkMappingDimension];
Int_t AliHLTTPCDigitReaderRaw::fgMapping3[fgkMapping3Size][fgkMappingDimension];
Int_t AliHLTTPCDigitReaderRaw::fgMapping4[fgkMapping4Size][fgkMappingDimension];
Int_t AliHLTTPCDigitReaderRaw::fgMapping5[fgkMapping5Size][fgkMappingDimension];
#endif //#if defined(HAVE_TPC_MAPPING)
