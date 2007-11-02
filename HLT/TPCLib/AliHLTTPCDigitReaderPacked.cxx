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

/** @file   AliHLTTPCDigitReaderPacked.cxx
    @author Timm Steinbeck, Jochen Thaeder, Matthias Richter, Kenneth Aamodt
    @date   
    @brief  A digit reader implementation for simulated, packed TPC 'raw' data.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCDigitReaderPacked.h"

#include "AliTPCRawStream.h"
#include "AliRawReaderMemory.h"
#include "AliRawDataHeader.h"

#if ENABLE_PAD_SORTING
#include "AliHLTTPCTransform.h"
#endif // ENABLE_PAD_SORTING
#include "AliHLTStdIncludes.h"

ClassImp(AliHLTTPCDigitReaderPacked)

#if defined(HAVE_ALIRAWDATA) && defined(HAVE_ALITPCRAWSTREAM_H)
AliHLTTPCDigitReaderPacked::AliHLTTPCDigitReaderPacked()
  :
#if ENABLE_PAD_SORTING
  fCurrentRow(0),
  fCurrentPad(0),
  fCurrentBin(-1),
  fNRows(0),
  fRowOffset(0),
  fNMaxRows(0),
  fNMaxPads(0),
  fNTimeBins(0),
  fData(NULL),
#endif // ENABLE_PAD_SORTING  
  fRawMemoryReader(NULL),
  fTPCRawStream(NULL),
  fOldRCUFormat(kFALSE)
{
  fRawMemoryReader = new AliRawReaderMemory;
  
  fTPCRawStream = new AliTPCRawStream( fRawMemoryReader );

#if ENABLE_PAD_SORTING
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

  HLTDebug("Array Borders  ||| MAXPAD=%d ||| MAXROW=%d ||| MAXBIN=%d ||| MAXMUL=%d", 
	   fNMaxPads, fNMaxRows, fNTimeBins, fNTimeBins*fNMaxRows*fNMaxPads);

  // init Data array
  fData = new Int_t[ fNMaxRows*fNMaxPads*fNTimeBins ];
#endif // ENABLE_PAD_SORTING
}

AliHLTTPCDigitReaderPacked::~AliHLTTPCDigitReaderPacked(){
  if ( fRawMemoryReader )
    delete fRawMemoryReader;
  fRawMemoryReader = NULL;
  if ( fTPCRawStream )
      delete fTPCRawStream;
  fTPCRawStream = NULL;
#if ENABLE_PAD_SORTING 
  if ( fData )
      delete [] fData;
  fData = NULL;
#endif // ENABLE_PAD_SORTING
}

Int_t AliHLTTPCDigitReaderPacked::InitBlock(void* ptr,unsigned long size, Int_t patch, Int_t slice){

  fRawMemoryReader->SetMemory( reinterpret_cast<UChar_t*>( ptr ), size );

  //get DDL ID in order to tell the memory reader which slice/patch to use
  Int_t DDLid= 0;
  if (patch < 2)
    DDLid = 768 + 2*slice + patch;
  else 
    DDLid = 840 + 4*slice + patch-2;

  fRawMemoryReader->SetEquipmentID(DDLid);
  //fRawMemoryReader->SetEquipmentID(1);
  if(fOldRCUFormat)
    fTPCRawStream->SetOldRCUFormat(kTRUE);

#if ENABLE_PAD_SORTING

  fCurrentRow = 0;
  fCurrentPad = 0;
  fCurrentBin = -1;

  Int_t firstrow=AliHLTTPCTransform::GetFirstRow(patch);
  Int_t lastrow=AliHLTTPCTransform::GetLastRow(patch);
  fNRows = lastrow - firstrow + 1;

  Int_t offset=0;
  if (patch > 1) offset =  AliHLTTPCTransform::GetFirstRow( 2 );

  fRowOffset = firstrow - offset;
  firstrow -= offset;
  lastrow  -= offset;

  // Init array with -1
  memset( fData, 0xFF, sizeof(Int_t)*(fNMaxRows*fNMaxPads*fNTimeBins) );

  // read data and fill in array
  while( fTPCRawStream->Next()){

      Int_t row = fTPCRawStream->GetRow();
      Int_t pad = fTPCRawStream->GetPad();
      Int_t bin = fTPCRawStream->GetTime();

      if ( row < firstrow || row > lastrow || pad > AliHLTTPCTransform::GetNPads(row + offset) || bin > fNTimeBins){
	HLTFatal("Index out of Range Probably wrong patch! %d - %d", slice, patch);
	if ( row < firstrow || row > lastrow ) 
	  HLTFatal("Row out of Range %d < %d < %d",firstrow, row, lastrow);
	if ( pad > AliHLTTPCTransform::GetNPads(row + offset) ) 
	  HLTFatal("Pad out of Range %d < %d < %d",pad, AliHLTTPCTransform::GetNPads(row + offset));
	if ( bin > fNTimeBins )
	  HLTFatal("Bin out of Range %d < %d < %d",bin, fNTimeBins);
      }
      else {  
	  if ((row-fRowOffset)*fNMaxPads*fNTimeBins+ pad*fNTimeBins + bin >=  fNMaxRows*fNMaxPads*fNTimeBins ) {
	      HLTFatal("Index out of array range PAD=%d ||| ROW=%d ||| BIN=%d ||| OFFSET=%d ||| ROWOFFSET=%d", pad, row, bin, offset, fRowOffset);
	      continue;
	  }
	  else {
	      fData[ (row-fRowOffset)*fNMaxPads*fNTimeBins+ pad*fNTimeBins + bin ] = fTPCRawStream->GetSignal() ;
	  }
      }
  }
#endif // ENABLE_PAD_SORTING

  return 0;
}

Bool_t AliHLTTPCDigitReaderPacked::Next(){
  Bool_t readvalue = kTRUE;

#if ENABLE_PAD_SORTING
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
	  HLTFatal("Overflow: row=%d pad=%d bin=%d", fCurrentRow, fCurrentPad, fCurrentBin);
	  readvalue = kFALSE;
	  break;
      }

      if (fData[ fCurrentRow*fNMaxPads*fNTimeBins + fCurrentPad*fNTimeBins + fCurrentBin  ] != -1) break;
  }
#else // !ENABLE_PAD_SORTING
  readvalue = fTPCRawStream->Next();
#endif // ENABLE_PAD_SORTING

  return readvalue;
}

Int_t AliHLTTPCDigitReaderPacked::GetRow(){
#if ENABLE_PAD_SORTING
  return (fCurrentRow + fRowOffset);
#else // !ENABLE_PAD_SORTING
  return (Int_t) fTPCRawStream->GetRow();
#endif // ENABLE_PAD_SORTING
}

int AliHLTTPCDigitReaderPacked::GetPad(){
#if ENABLE_PAD_SORTING
  return fCurrentPad;
#else // !ENABLE_PAD_SORTING
  return fTPCRawStream->GetPad();
#endif // ENABLE_PAD_SORTING
}

Int_t AliHLTTPCDigitReaderPacked::GetSignal(){ 
#if ENABLE_PAD_SORTING
  return fData[ fCurrentRow*fNMaxPads*fNTimeBins+ fCurrentPad*fNTimeBins + fCurrentBin ];
#else // !ENABLE_PAD_SORTING
  return fTPCRawStream->GetSignal();
#endif // ENABLE_PAD_SORTING
}

Int_t AliHLTTPCDigitReaderPacked::GetTime(){
#if ENABLE_PAD_SORTING
  return fCurrentBin;
#else // !ENABLE_PAD_SORTING
  return fTPCRawStream->GetTime();
#endif // ENABLE_PAD_SORTING
}
#endif //defined(HAVE_ALIRAWDATA) && defined(HAVE_ALITPCRAWSTREAM_H)
