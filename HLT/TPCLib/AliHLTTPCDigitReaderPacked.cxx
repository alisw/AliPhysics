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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for reading packed data for the HLT                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if __GNUC__== 3
using namespace std;
#endif

#if defined(HAVE_ALIRAWDATA) && defined(HAVE_ALITPCRAWSTREAM_H)
#include "AliHLTTPCDigitReaderPacked.h"
#include "AliTPCRawStream.h"
#include "AliRawReaderMemory.h"
#include "AliRawDataHeader.h"

#include "AliHLTStdIncludes.h"

ClassImp(AliHLTTPCDigitReaderPacked)

AliHLTTPCDigitReaderPacked::AliHLTTPCDigitReaderPacked(){
  fRawMemoryReader = new AliRawReaderMemory;
  fTPCRawStream = new AliTPCRawStream( fRawMemoryReader );
}

AliHLTTPCDigitReaderPacked::~AliHLTTPCDigitReaderPacked(){
  if ( fRawMemoryReader )
    delete fRawMemoryReader;
  fRawMemoryReader = NULL;
  if ( fTPCRawStream )
    delete fTPCRawStream;
  fTPCRawStream = NULL;
}

int AliHLTTPCDigitReaderPacked::InitBlock(void* ptr,unsigned long size,Int_t firstrow, Int_t lastrow){
  fRawMemoryReader->SetMemory( reinterpret_cast<UChar_t*>( ptr ), size );
  return 0;
}

bool AliHLTTPCDigitReaderPacked::Next(){
  bool rreadvalue;
  rreadvalue = fTPCRawStream->Next();
  return rreadvalue;
}

int AliHLTTPCDigitReaderPacked::GetRow(){
  int rrow;
  rrow = (int)fTPCRawStream->GetRow();
  return rrow;
}

int AliHLTTPCDigitReaderPacked::GetPad(){
  int rpad;
  rpad = fTPCRawStream->GetPad();
  return rpad   ;
}

int AliHLTTPCDigitReaderPacked::GetSignal(){ 
  int rsignal;
  rsignal = fTPCRawStream->GetSignal();
  return rsignal;
}

int AliHLTTPCDigitReaderPacked::GetTime(){
  int rtime;
  rtime = fTPCRawStream->GetTime();
  return rtime;
}
#endif //defined(HAVE_ALIRAWDATA) && defined(HAVE_ALITPCRAWSTREAM_H)
