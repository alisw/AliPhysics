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


const Int_t  AliPMDDspHeader::fgkHeaderLength = 8;

//------------------------------------------------------------
AliPMDDspHeader::AliPMDDspHeader()
  :  TObject(),
     fTotalLength(0),
     fRawDataLength(0),
     fTrWord1(0),
     fTrWord2(0),
     fTrWord3(0),
     fTrWord4(0),
     fDspId(0),
     fEvtWord(0)
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
AliPMDDspHeader::AliPMDDspHeader(const AliPMDDspHeader & dsph): TObject()
{
  //
  // copy ctor
  //

  fTotalLength   = dsph.fTotalLength;
  fRawDataLength = dsph.fRawDataLength;
  fTrWord1       = dsph.fTrWord1;
  fTrWord2       = dsph.fTrWord2;
  fTrWord3       = dsph.fTrWord3;
  fTrWord4       = dsph.fTrWord4;
  fDspId         = dsph.fDspId;
  fEvtWord       = dsph.fEvtWord;
}

//___________________________________________
AliPMDDspHeader& AliPMDDspHeader::operator=(const AliPMDDspHeader &dsph)
{
  // 
  // assignment operator
  //
  if (this == &dsph) return *this;

  fTotalLength   = dsph.fTotalLength;
  fRawDataLength = dsph.fRawDataLength;
  fTrWord1       = dsph.fTrWord1;
  fTrWord2       = dsph.fTrWord2;
  fTrWord3       = dsph.fTrWord3;
  fTrWord4       = dsph.fTrWord4;
  fDspId         = dsph.fDspId;
  fEvtWord       = dsph.fEvtWord;

  return *this;
}
void AliPMDDspHeader::SetHeader(Int_t *header)
{
  fTotalLength   = header[0];
  fRawDataLength = header[1];
  fTrWord1       = header[2];
  fTrWord2       = header[3];
  fTrWord3       = header[4];
  fTrWord4       = header[5];
  fDspId         = header[6];
  fEvtWord       = header[7];
}
      
