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

 
#include "AliPMDPatchBusHeader.h"



ClassImp(AliPMDPatchBusHeader)


const Int_t  AliPMDPatchBusHeader::fgkHeaderLength = 4;

//------------------------------------------------------------
AliPMDPatchBusHeader::AliPMDPatchBusHeader()
  :  TObject(),
     fTotalLength(0),
     fRawDataLength(0),
     fPatchBusId(0),
     fTrWord(0)
{
  //
  // ctor
  //

}

//___________________________________________
AliPMDPatchBusHeader::~AliPMDPatchBusHeader()
{
  // 
  // dtor
  //
}

//___________________________________________
AliPMDPatchBusHeader::AliPMDPatchBusHeader(const AliPMDPatchBusHeader & pbush):
  TObject(),
  fTotalLength(pbush.fTotalLength),
  fRawDataLength(pbush.fRawDataLength),
  fPatchBusId(pbush.fPatchBusId),
  fTrWord(pbush.fTrWord)
{
  //
  // copy ctor
  //
}

//___________________________________________
AliPMDPatchBusHeader&
AliPMDPatchBusHeader::operator=(const AliPMDPatchBusHeader &pbush)
{
  // 
  // assignment operator
  //
  if (this != &pbush)
    {
      fTotalLength   = pbush.fTotalLength;
      fRawDataLength = pbush.fRawDataLength;
      fPatchBusId    = pbush.fPatchBusId;
      fTrWord        = pbush.fTrWord;
    }
  return *this;
}
void AliPMDPatchBusHeader::SetHeader(Int_t *header)
{
  fTotalLength   = header[0];
  fRawDataLength = header[1];
  fPatchBusId    = header[2];
  fTrWord        = header[3];
}

