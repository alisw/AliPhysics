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

 
#include "AliPMDBlockHeader.h"



ClassImp(AliPMDBlockHeader)


const Int_t  AliPMDBlockHeader::fgkHeaderLength = 8;

//------------------------------------------------------------
AliPMDBlockHeader::AliPMDBlockHeader()
  :  TObject(),
     fTotalLength(0),
     fRawDataLength(0),
     fDspId(0),
     fTrWord1(0),
     fTrWord2(0),
     fTrWord3(0),
     fTrWord4(0),
     fPadWord(0)
{
  //
  // ctor
  //

}

//___________________________________________
AliPMDBlockHeader::~AliPMDBlockHeader()
{
  // 
  // dtor
  //
}

//___________________________________________
AliPMDBlockHeader::AliPMDBlockHeader(const AliPMDBlockHeader & blockh):
  TObject(),
  fTotalLength(blockh.fTotalLength),
  fRawDataLength(blockh.fRawDataLength),
  fDspId(blockh.fDspId),
  fTrWord1(blockh.fTrWord1),
  fTrWord2(blockh.fTrWord2),
  fTrWord3(blockh.fTrWord3),
  fTrWord4(blockh.fTrWord4),
  fPadWord(blockh.fPadWord)
{
  //
  // copy ctor
  //

}

//___________________________________________
AliPMDBlockHeader&
AliPMDBlockHeader::operator=(const AliPMDBlockHeader &blockh)
{
  // 
  // assignment operator
  //
  if (this != &blockh)
    {
      fTotalLength   = blockh.fTotalLength;
      fRawDataLength = blockh.fRawDataLength;
      fDspId         = blockh.fDspId;
      fTrWord1       = blockh.fTrWord1;
      fTrWord2       = blockh.fTrWord2;
      fTrWord3       = blockh.fTrWord3;
      fTrWord4       = blockh.fTrWord4;
      fPadWord       = blockh.fPadWord;
    }
  return *this;
}
void AliPMDBlockHeader::SetHeader(Int_t *header)
{
  fTotalLength   = header[0];
  fRawDataLength = header[1];
  fDspId         = header[2];
  fTrWord1       = header[3];
  fTrWord2       = header[4];
  fTrWord3       = header[5];
  fTrWord4       = header[6];
  fPadWord       = header[7];

}
      
