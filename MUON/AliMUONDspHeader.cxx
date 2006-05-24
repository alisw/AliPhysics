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
 
#include "AliMUONDspHeader.h"
#include "AliMUONBusStruct.h"

/// 
/// DSP structure for tracker raw data.
/// Each block contains at most 5 Dsp structures.
/// Beside the total length and length of the below data
/// the header of the Dsp contains the front end DSP id, trigger words
/// and event word (1 for nb of word is odd and 0 if not 
///

/// \cond CLASSIMP
ClassImp(AliMUONDspHeader)
/// \endcond
  
  const Int_t AliMUONDspHeader::fgkHeaderLength = 8;

//___________________________________________
AliMUONDspHeader::AliMUONDspHeader()
  :  TObject(),
     fTotalLength(0),
     fLength(0),
     fDspId(0),
     fEventWord(0)
{
  //
  //ctor
  //
  for (Int_t i = 0; i < 4; i++)
    fTriggerWord[i] = 0;

  fBusPatchArray  = new TClonesArray("AliMUONBusStruct",5);

}

//___________________________________________
AliMUONDspHeader::~AliMUONDspHeader()
{
  //
  // dtr
  //
  fBusPatchArray->Delete();
  delete fBusPatchArray;
}

//___________________________________________
AliMUONDspHeader::AliMUONDspHeader(const AliMUONDspHeader& event)
  :  TObject(event)
{
  // 
  // copy constructor
  //
  fTotalLength = event.fTotalLength;
  fLength      = event.fLength;
  fDspId       = event.fDspId;
  fEventWord   = event.fEventWord;

  //ctor
  for (Int_t i = 0; i < 4; i++)
    fTriggerWord[i] = event.fTriggerWord[i];

  fBusPatchArray = new TClonesArray("AliMUONBusStruct", 5);
  for (Int_t index = 0; index < (event.fBusPatchArray)->GetEntriesFast(); index++) {
    {new ((*fBusPatchArray)[fBusPatchArray->GetEntriesFast()]) 
        AliMUONBusStruct(*(AliMUONBusStruct*)(event.fBusPatchArray)->At(index));}
  }
  // fBusPatchArray->SetOwner();
 
}

//___________________________________________
AliMUONDspHeader& AliMUONDspHeader::operator=(const AliMUONDspHeader& event)
{
  //
  // assignemnt constructor
  //
  if (this == &event) return *this;

  fTotalLength = event.fTotalLength;
  fLength      = event.fLength;
  fDspId       = event.fDspId;
  fEventWord   = event.fEventWord;

  //ctor
  for (Int_t i = 0; i < 4; i++)
    fTriggerWord[i] = event.fTriggerWord[i];

  fBusPatchArray = new TClonesArray("AliMUONBusStruct", 5);
  for (Int_t index = 0; index < (event.fBusPatchArray)->GetEntriesFast(); index++) {
    {new ((*fBusPatchArray)[fBusPatchArray->GetEntriesFast()]) 
        AliMUONBusStruct(*(AliMUONBusStruct*)(event.fBusPatchArray)->At(index));}
  }
  return *this;
}
//___________________________________________
void AliMUONDspHeader::AddBusPatch(const AliMUONBusStruct& busPatch)
{
  // 
  // adding buspatch info
  // into TClonesArray
  //
  TClonesArray &eventArray = *fBusPatchArray;
  new(eventArray[eventArray.GetEntriesFast()]) AliMUONBusStruct(busPatch);
}
//___________________________________________
void AliMUONDspHeader::Clear(Option_t* )
{
  // Clear TClones arrays
  // instead of deleting
  //
  fBusPatchArray->Clear("C");
 
}
