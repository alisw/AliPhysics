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
     fDataKey(0),
     fTotalLength(0),
     fRawDataLength(0),
     fDspId(0),
     fL0Trigger(0),
     fMiniEventId(0),
     fEventId1(0),
     fEventId2(0)
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
  fDataKey(blockh.fDataKey),
  fTotalLength(blockh.fTotalLength),
  fRawDataLength(blockh.fRawDataLength),
  fDspId(blockh.fDspId),
  fL0Trigger(blockh.fL0Trigger),
  fMiniEventId(blockh.fMiniEventId),
  fEventId1(blockh.fEventId1),
  fEventId2(blockh.fEventId2)
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
      fDataKey       = blockh.fDataKey;
      fTotalLength   = blockh.fTotalLength;
      fRawDataLength = blockh.fRawDataLength;
      fDspId         = blockh.fDspId;
      fL0Trigger     = blockh.fL0Trigger;
      fMiniEventId   = blockh.fMiniEventId;
      fEventId1      = blockh.fEventId1;
      fEventId2      = blockh.fEventId2;
    }
  return *this;
}
void AliPMDBlockHeader::SetHeader(Int_t *header)
{
  fDataKey       = header[0];
  fTotalLength   = header[1];
  fRawDataLength = header[2];
  fDspId         = header[3];
  fL0Trigger     = header[4];
  fMiniEventId   = header[5];
  fEventId1      = header[6];
  fEventId2      = header[7];
}
      
