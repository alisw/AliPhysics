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

 
#include "AliPMDDspHeader.h"



ClassImp(AliPMDDspHeader)


const Int_t   AliPMDDspHeader::fgkHeaderLength = 10;
const UInt_t  AliPMDDspHeader::fgkDefaultPaddingWord = 0xFFFFFFFF;

//------------------------------------------------------------
AliPMDDspHeader::AliPMDDspHeader() :
  TObject(),
  fDataKey(0),
  fTotalLength(0),
  fRawDataLength(0),
  fDspId(0),
  fBlkL1ATrigger(0),
  fMiniEventId(0),
  fL1ATrigger(0),
  fL1RTrigger(0),
  fPaddingWord(0),
  fErrorWord(0)
{
  //
  // ctor
  //

}

//___________________________________________
AliPMDDspHeader::~AliPMDDspHeader()
{
  // 
  // dtor
  //
}

//___________________________________________
AliPMDDspHeader::AliPMDDspHeader(const AliPMDDspHeader & dsph):
  TObject(),
  fDataKey(dsph.fDataKey),
  fTotalLength(dsph.fTotalLength),
  fRawDataLength(dsph.fRawDataLength),
  fDspId(dsph.fDspId),
  fBlkL1ATrigger(dsph.fBlkL1ATrigger),
  fMiniEventId(dsph.fMiniEventId),
  fL1ATrigger(dsph.fL1ATrigger),
  fL1RTrigger(dsph.fL1RTrigger),
  fPaddingWord(dsph.fPaddingWord),
  fErrorWord(dsph.fErrorWord)
{
  //
  // copy ctor
  //
}

//___________________________________________
AliPMDDspHeader& AliPMDDspHeader::operator=(const AliPMDDspHeader &dsph)
{
  // 
  // assignment operator
  //
  if (this != &dsph)
    {
      fDataKey       = dsph.fDataKey;
      fTotalLength   = dsph.fTotalLength;
      fRawDataLength = dsph.fRawDataLength;
      fDspId         = dsph.fDspId;
      fBlkL1ATrigger = dsph.fBlkL1ATrigger;
      fMiniEventId   = dsph.fMiniEventId;
      fL1ATrigger    = dsph.fL1ATrigger;
      fL1RTrigger    = dsph.fL1RTrigger;
      fPaddingWord   = dsph.fPaddingWord;
      fErrorWord     = dsph.fErrorWord;
    }
  return *this;
}
void AliPMDDspHeader::SetHeader(Int_t *header)
{
  fDataKey        = header[0];
  fTotalLength    = header[1];
  fRawDataLength  = header[2];
  fDspId          = header[3];
  fBlkL1ATrigger  = header[4];
  fMiniEventId    = header[5];
  fL1ATrigger     = header[6];
  fL1RTrigger     = header[7];
  fPaddingWord    = header[8];
  fErrorWord      = header[9];
}
      
