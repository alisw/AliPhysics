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

/* $Id$ */
 
#include "AliMUONBusStruct.h"
#include "AliLog.h"


/// 
/// Bus patch structure for tracker raw data
/// each Dsp contains at most 5 bus patch structure
/// Beside the total length and length of the below data
/// the header of the block contains the bus patch id, trigger words 
/// and data structure itself (11bits for manu id, 6 bits for channel id and
/// 12 bits for charge)
///

/// \cond CLASSIMP
ClassImp(AliMUONBusStruct)
/// \endcond

const Int_t  AliMUONBusStruct::fgkHeaderLength = 4;
const UInt_t AliMUONBusStruct::fgkDefaultDataKey = 0xB000000B;

//___________________________________________
AliMUONBusStruct::AliMUONBusStruct()
  :  TObject(),
     fDataKey(0),
     fTotalLength(0),
     fLength(0),
     fBusPatchId(0),
     fBufSize(43*64), 
  /* assuming 43 manus max per bustpatch.
  Anyway fData is resized where needed (though it makes it slower) */
     fData(new UInt_t[fBufSize]),
     fDspId(0),
     fBlkId(0)
{
  //
  // ctor
  // 

}
//___________________________________________
AliMUONBusStruct::~AliMUONBusStruct()
{
  //
  // dtor
  //
  delete[] fData;
}

//___________________________________________
void AliMUONBusStruct::SetAlloc(Int_t size)
{
  //
  // Allocate size 
  // return if size < fBufSize 
  //
  if (size < fBufSize) 
    return;
  else 
    ResizeData(size);
}
//___________________________________________
void AliMUONBusStruct::AddData(UInt_t data)
{
  // could have used class from ROOT
  // but the structure must be as simple as possible
  // to be written on disc blockwise, not so sure ?
  if (fLength == fBufSize) 
    ResizeData();
  fData[fLength++] = data;
  fTotalLength = fLength + fgkHeaderLength;
}

//___________________________________________
void AliMUONBusStruct::ResizeData(Int_t size)
{
  // In case of resizing the vector
  // the most simplest way to do it
  //
  AliInfo("reallocating");
  if (size == 0)
    fBufSize *= 2;
  else
    fBufSize = size;
  UInt_t* newData = new UInt_t[fBufSize];
  for (Int_t i = 0; i < fLength; i++)
    newData[i] = fData[i];
  delete[] fData;
  fData = newData;
}
//___________________________________________
AliMUONBusStruct::
AliMUONBusStruct(const AliMUONBusStruct& event)
  : TObject(event),
    fDataKey(event.fDataKey),
    fTotalLength(event.fTotalLength),
    fLength(event.fLength),
    fBusPatchId(event.fBusPatchId),
    fBufSize(event.fBufSize),
    fData(new UInt_t[event.fBufSize]),
    fDspId(event.fDspId),
    fBlkId(event.fBlkId)
{
  //
  // copy ctor
  //

  for (int i = 0; i < event.fBufSize; i++)
    fData[i] = event.fData[i];
}
//___________________________________________
AliMUONBusStruct&
AliMUONBusStruct::operator=(const AliMUONBusStruct& event)
{
  //
  // assignment operator
  //
  if (this == &event) return *this;
  fDataKey     = event.fDataKey;
  fTotalLength = event.fTotalLength;
  fLength      = event.fLength;
  fBusPatchId  = event.fBusPatchId;
  fBufSize     = event.fBufSize;

  fBlkId = event.fBlkId;
  fDspId = event.fDspId;

  delete [] fData;  
  fData =  new UInt_t[event.fBufSize];
  for (int i = 0; i < event.fLength; i++)
    fData[i] = event.fData[i];

  return *this;
}
//___________________________________________
Int_t AliMUONBusStruct::Compare(const TObject *obj) const
{
  // 
  // sort bus patch by bus patch number
  // important for AliMUONRawWriter
  //
  AliMUONBusStruct* event = (AliMUONBusStruct*) obj;
  return (fBusPatchId > event->GetBusPatchId()) ? 1 : -1;
}

//___________________________________________
void AliMUONBusStruct::Clear(Option_t *)
{
  // clear
  // delete the allocated memory 
  //
  AliInfo("here");
  delete[] fData;
}
//___________________________________________
UInt_t AliMUONBusStruct::GetData(Int_t n) const 
{
  //
  // get data
  //
  if ( n>=0 && n<fLength ) return fData[n];

  AliError("Index outside limits."); 
  return 0; 
}

//___________________________________________
Char_t AliMUONBusStruct::GetParity(Int_t n) const   
{
  //
  // get parity
  //
  if ( n>=0 && n<fLength ) return (Char_t)(fData[n] >> 31) &  0x1;

  AliError("Index outside limits."); 
  return 0; 
}

//___________________________________________
UShort_t AliMUONBusStruct::GetManuId(Int_t n) const     
{
  //
  // get manu Id
  //
  if ( n>=0 && n<fLength ) return (UShort_t)(fData[n] >> 18) &  0x7FF;

  AliError("Index outside limits."); 
  return 0; 
}

//___________________________________________
Char_t AliMUONBusStruct::GetChannelId(Int_t n) const  
{
  // 
  // get channel Id
  //
  if ( n>=0 && n<fLength ) return (Char_t)(fData[n] >> 12) & 0x3F;

  AliError("Index outside limits."); 
  return 0; 
}

//___________________________________________
UShort_t AliMUONBusStruct::GetCharge(Int_t n) const     
{
  //
  // get charge (in ADC)
  //
  if ( n>=0 && n<fLength ) return (UShort_t)(fData[n] & 0xFFF);

  AliError("Index outside limits."); 
  return 0; 
}
