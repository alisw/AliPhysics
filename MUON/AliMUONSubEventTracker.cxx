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

#include "AliMUONSubEventTracker.h"

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
  //ctor
  fData = new UInt_t[fBufSize];
}
//___________________________________________
AliMUONSubEventTracker::~AliMUONSubEventTracker()
{
  if(fData)
    delete[] fData;
}

//___________________________________________
void AliMUONSubEventTracker::SetAlloc(Int_t size)
{
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
  // to be written on disc blockwise
  if (fLength == fBufSize) 
    ResizeData();
  fData[fLength++] = data;
  fTotalLength = fLength + 4;
}

//___________________________________________
void AliMUONSubEventTracker::ResizeData(Int_t size)
{
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
  fTotalLength = event.fTotalLength;
  fLength = event.fLength;
  fBusPatchId = event.fBusPatchId;
  fTriggerWord = event.fTriggerWord;

  fData =  new UInt_t[event.fBufSize];
  for (int i = 0; i < event.fBufSize; i++)
    fData[i] = event.fData[i];
}
//___________________________________________
AliMUONSubEventTracker&
AliMUONSubEventTracker::operator=(const AliMUONSubEventTracker& event)
{
  if (this == &event) return *this;
  fTotalLength = event.fTotalLength;
  fLength = event.fLength;
  fBusPatchId = event.fBusPatchId;
  fTriggerWord = event.fTriggerWord;
  
  for (int i = 0; i < event.fLength; i++)
    fData[i] = event.fData[i];

  return *this;
}
//___________________________________________
Int_t AliMUONSubEventTracker::Compare(const TObject *obj) const
{
  AliMUONSubEventTracker* event = (AliMUONSubEventTracker*) obj;
  return (fBusPatchId > event->GetBusPatchId()) ? 1 : -1;
}
