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

// $Id$
 
#include "AliMUONSubEventTracker.h"

#include "AliLog.h"

/// 
/// Bus patch structure for tracker raw data
/// each Dsp contains at most 5 bus patch structure
/// Beside the total length and length of the below data
/// the header of the block contains the bus patch id, trigger words 
/// and data structure itself (11bits for manu id, 6 bits for channel id and
/// 12 bits for charge)
///


const Int_t AliMUONSubEventTracker::fgkHeaderLength = 4;

ClassImp(AliMUONSubEventTracker)

//___________________________________________
AliMUONSubEventTracker::AliMUONSubEventTracker()
  :  TObject(),
     fTotalLength(0),
     fLength(0),
     fBusPatchId(0),
     fTriggerWord(0),
     fBufSize(1024)
{
  //
  //ctor
  //
  fData = new UInt_t[fBufSize];
}
//___________________________________________
AliMUONSubEventTracker::~AliMUONSubEventTracker()
{
  //
  // dtor
  //
  delete[] fData;
}

//___________________________________________
void AliMUONSubEventTracker::SetAlloc(Int_t size)
{
  //
  // Allocate size per default 1024;
  // return if size < 1024
  //
  if (size < fBufSize) 
    return;
  else 
    ResizeData(size);
}
//___________________________________________
void AliMUONSubEventTracker::AddData(UInt_t data)
{
  // could have used class from ROOT
  // but the structure must be as simple as possible
  // to be written on disc blockwise, not so sure ?
  if (fLength == fBufSize) 
    ResizeData();
  fData[fLength++] = data;
  fTotalLength = fLength + 4;
}

//___________________________________________
void AliMUONSubEventTracker::ResizeData(Int_t size)
{
  // In case of resizing the vector
  // the most simplest way to do it
  //
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
AliMUONSubEventTracker::
AliMUONSubEventTracker(const AliMUONSubEventTracker& event): TObject(event)
{
  //
  // copy ctor
  //
  fTotalLength = event.fTotalLength;
  fLength      = event.fLength;
  fBusPatchId  = event.fBusPatchId;
  fTriggerWord = event.fTriggerWord;
  fBufSize     = event.fBufSize;

  fData =  new UInt_t[event.fBufSize];
  for (int i = 0; i < event.fBufSize; i++)
    fData[i] = event.fData[i];
}
//___________________________________________
AliMUONSubEventTracker&
AliMUONSubEventTracker::operator=(const AliMUONSubEventTracker& event)
{
  //
  // assignment operator
  //
  if (this == &event) return *this;
  fTotalLength = event.fTotalLength;
  fLength      = event.fLength;
  fBusPatchId  = event.fBusPatchId;
  fTriggerWord = event.fTriggerWord;
  fBufSize     = event.fBufSize;

  delete [] fData;  
  fData =  new UInt_t[event.fBufSize];
  for (int i = 0; i < event.fBufSize; i++)
    fData[i] = event.fData[i];

  return *this;
}
//___________________________________________
Int_t AliMUONSubEventTracker::Compare(const TObject *obj) const
{
  // 
  // sort bus patch by bus patch number
  // important for AliMUONRawWriter
  //
  AliMUONSubEventTracker* event = (AliMUONSubEventTracker*) obj;
  return (fBusPatchId > event->GetBusPatchId()) ? 1 : -1;
}

//___________________________________________
UInt_t  AliMUONSubEventTracker::GetData(Int_t n) const 
{
  //
  // get data
  //
  if ( n>=0 && n<fLength ) return fData[n];

  AliError("Index outside limits."); 
  return 0; 
}

//___________________________________________
Char_t   AliMUONSubEventTracker::GetParity(Int_t n) const   
{
  //
  // get pariry
  //
  if ( n>=0 && n<fLength ) return (Char_t)(fData[n] >> 29) &  0x7;

  AliError("Index outside limits."); 
  return 0; 
}

//___________________________________________
UShort_t AliMUONSubEventTracker::GetManuId(Int_t n) const     
{
  //
  // get manu Id
  //
  if ( n>=0 && n<fLength ) return (UShort_t)(fData[n] >> 18) &  0x7FF;

  AliError("Index outside limits."); 
  return 0; 
}

//___________________________________________
Char_t   AliMUONSubEventTracker::GetChannelId(Int_t n) const  
{
  // 
  // get channel Id
  //
  if ( n>=0 && n<fLength ) return (Char_t)(fData[n] >> 12) & 0x3F;

  AliError("Index outside limits."); 
  return 0; 
}

//___________________________________________
UShort_t AliMUONSubEventTracker::GetCharge(Int_t n) const     
{
  //
  // get charge (in ADC)
  //
  if ( n>=0 && n<fLength ) return (UShort_t)(fData[n] & 0xFFF);

  AliError("Index outside limits."); 
  return 0; 
}
