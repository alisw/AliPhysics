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
 
#include "AliMUONBlockHeader.h"
#include "AliMUONDspHeader.h"

/// \class AliMUONBlockHeader
/// Block structure for tracker raw data
/// each DDL contains two blocks,
/// each block contains at most 5 dsp structure.
/// Beside the total length and length of the below data
/// the header of the block contains the front end DSP id, trigger words and paddind word
///
/// \author Christian Finck

/// \cond CLASSIMP
ClassImp(AliMUONBlockHeader)
/// \endcond

const Int_t  AliMUONBlockHeader::fgkHeaderLength = 8;
const UInt_t AliMUONBlockHeader::fgkDefaultDataKey = 0xFC0000FC;
//___________________________________________
AliMUONBlockHeader::AliMUONBlockHeader()
  :  TObject(),
     fDataKey(0),
     fTotalLength(0),
     fLength(0),
     fDspId(0),
     fL0Trigger(0),
     fMiniEventId(0),
     fEventId1(0),
     fEventId2(0),
     fDspHeaderArray(new TClonesArray("AliMUONDspHeader", 5))

{
  ///
  /// ctor
  ///

}

//___________________________________________
AliMUONBlockHeader::~AliMUONBlockHeader()
{
  /// 
  /// dtor
  ///
  fDspHeaderArray->Delete();
  delete fDspHeaderArray;
}

//___________________________________________
AliMUONBlockHeader::AliMUONBlockHeader(const AliMUONBlockHeader& event)
  :  TObject(event),
     fDataKey(event.fDataKey),
     fTotalLength(event.fTotalLength),
     fLength(event.fLength),
     fDspId(event.fDspId),
     fL0Trigger(event.fL0Trigger),
     fMiniEventId(event.fMiniEventId),
     fEventId1(event.fEventId1),
     fEventId2(event.fEventId2),
     fDspHeaderArray(new TClonesArray("AliMUONDspHeader", 5))
{
  ///
  /// copy ctor
  ///

  for (Int_t index = 0; index < (event.fDspHeaderArray)->GetEntriesFast(); index++) {
    {new ((*fDspHeaderArray)[fDspHeaderArray->GetEntriesFast()]) 
        AliMUONDspHeader(*(AliMUONDspHeader*)(event.fDspHeaderArray)->At(index));}
  }
  //  fDspHeaderArray->SetOwner();
}

//___________________________________________
AliMUONBlockHeader&
AliMUONBlockHeader::operator=(const AliMUONBlockHeader &event)
{
  /// 
  /// assignment operator
  ///
  if (this == &event) return *this;

  fTotalLength = event.fTotalLength;
  fLength      = event.fLength;
  fDspId       = event.fDspId;
 
  fL0Trigger   = event.fL0Trigger;
  fMiniEventId = event.fMiniEventId;
  fEventId1    = event.fEventId1;
  fEventId2    = event.fEventId2;

  fDspHeaderArray = new TClonesArray("AliMUONDspHeader", 5);
  for (Int_t index = 0; index < (event.fDspHeaderArray)->GetEntriesFast(); index++) {
    new ((*fDspHeaderArray)[fDspHeaderArray->GetEntriesFast()]) 
        AliMUONDspHeader(*(AliMUONDspHeader*)(event.fDspHeaderArray)->At(index));
  }

  return *this;

}
//___________________________________________
void AliMUONBlockHeader::AddDspHeader(const AliMUONDspHeader& dspHeader)
{ 
  /// 
  /// adding the dsp structure
  /// into the TClonesArray
  ///
  TClonesArray &dspArray = *fDspHeaderArray;
  new(dspArray[dspArray.GetEntriesFast()]) AliMUONDspHeader(dspHeader);

}
//___________________________________________
void AliMUONBlockHeader::Clear(Option_t* )
{
  /// Clear TClones arrays
  /// instead of deleting
  ///
  fDspHeaderArray->Clear("C");
 
}
