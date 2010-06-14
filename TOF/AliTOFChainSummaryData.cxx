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
***************************************************************************/

/*
  author: Roberto Preghenella (R+), preghenella@bo.infn.it
*/


//////////////////////////////////////////////////////////////////////
//                                                                  //
//                                                                  //
//        This class provides a summary for TRM chain data.         //
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include "AliTOFChainSummaryData.h"

ClassImp(AliTOFChainSummaryData)

AliTOFChainSummaryData::AliTOFChainSummaryData() :
  TObject(),
  fHeader(kFALSE),
  fTrailer(kFALSE),
  fChain(0),
  fBunchID(0),
  fPB24Temp(0),
  fPB24ID(0),
  fTSBit(0),
  fStatus(0),
  fEventCounter(0),
  fTDCHitBuffer(0x0),
  fTDCPackedHitBuffer(0x0),
  fTDCErrorBuffer(0x0)
{
  /* default constructor */
  fTDCHitBuffer = new AliTOFTDCHitBuffer();
  fTDCPackedHitBuffer = new AliTOFTDCHitBuffer();
  fTDCErrorBuffer = new AliTOFTDCErrorBuffer();
}

//_________________________________________________________________

AliTOFChainSummaryData::AliTOFChainSummaryData(const AliTOFChainSummaryData &source) :
  TObject(),
  fHeader(source.fHeader),
  fTrailer(source.fTrailer),
  fChain(source.fChain),
  fBunchID(source.fBunchID),
  fPB24Temp(source.fPB24Temp),
  fPB24ID(source.fPB24ID),
  fTSBit(source.fTSBit),
  fStatus(source.fStatus),
  fEventCounter(source.fEventCounter),
  fTDCHitBuffer(0x0),
  fTDCPackedHitBuffer(0x0),
  fTDCErrorBuffer(0x0)
{
/* copy constructor */
  fTDCHitBuffer = new AliTOFTDCHitBuffer(*source.fTDCHitBuffer);
  fTDCPackedHitBuffer = new AliTOFTDCHitBuffer(*source.fTDCPackedHitBuffer);
  fTDCErrorBuffer = new AliTOFTDCErrorBuffer(*source.fTDCErrorBuffer);
}

//_________________________________________________________________

AliTOFChainSummaryData &
AliTOFChainSummaryData::operator = (const AliTOFChainSummaryData &source)
{
  /* operator = */
  fHeader = source.fHeader;
  fTrailer = source.fTrailer;
  fChain = source.fChain;
  fBunchID = source.fBunchID;
  fPB24Temp = source.fPB24Temp;
  fPB24ID = source.fPB24ID;
  fTSBit = source.fTSBit;
  fStatus = source.fStatus;
  fEventCounter = source.fEventCounter;
  *fTDCHitBuffer = *source.fTDCHitBuffer;
  *fTDCPackedHitBuffer = *source.fTDCPackedHitBuffer;
  *fTDCErrorBuffer = *source.fTDCErrorBuffer;
  return *this;
}

//_________________________________________________________________

AliTOFChainSummaryData::~AliTOFChainSummaryData()
{
  /* default destructor */
  delete fTDCHitBuffer;
  delete fTDCPackedHitBuffer;
  delete fTDCErrorBuffer;
}

//_________________________________________________________________

void
AliTOFChainSummaryData::Reset()
{
  /* reset function */
  fHeader = kFALSE;
  fTrailer = kFALSE;
  fChain = 0;
  fBunchID = 0;
  fPB24Temp = 0;
  fPB24ID = 0;
  fTSBit = 0;
  fStatus = 0;
  fEventCounter = 0;
  fTDCHitBuffer->Reset();
  fTDCPackedHitBuffer->Reset();
  fTDCErrorBuffer->Reset();
}
