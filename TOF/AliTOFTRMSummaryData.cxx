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
//        This class provides a summary for TRM data.               //
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include "AliTOFTRMSummaryData.h"

ClassImp(AliTOFTRMSummaryData)

AliTOFTRMSummaryData::AliTOFTRMSummaryData() :
  TObject(),
  fHeader(kFALSE),
  fTrailer(kFALSE),
  fSlotID(0),
  fEventWords(0),
  fACQBits(0),
  fLBit(0),
  fEBit(0),
  fEventCRC(0),
  fEventCounter(0),
  fDecoderCRC(0)
{
  /* default constructor */
  for (Int_t iChain = 0; iChain < N_CHAIN; iChain++)
    fChainSummaryData[iChain] = new AliTOFChainSummaryData();
}

//_________________________________________________________________

AliTOFTRMSummaryData::AliTOFTRMSummaryData(const AliTOFTRMSummaryData &source) :
  TObject(),
  fHeader(source.fHeader),
  fTrailer(source.fTrailer),
  fSlotID(source.fSlotID),
  fEventWords(source.fEventWords),
  fACQBits(source.fACQBits),
  fLBit(source.fLBit),
  fEBit(source.fEBit),
  fEventCRC(source.fEventCRC),
  fEventCounter(source.fEventCounter),
  fDecoderCRC(source.fDecoderCRC)
{
  /* copy constructor */
  for (Int_t iChain = 0; iChain < N_CHAIN; iChain++)
    fChainSummaryData[iChain] = new AliTOFChainSummaryData(*source.fChainSummaryData[iChain]);
}

//_________________________________________________________________

AliTOFTRMSummaryData &
AliTOFTRMSummaryData::operator = (const AliTOFTRMSummaryData &source)
{
  /* operator = */
  fHeader = source.fHeader;
  fTrailer = source.fTrailer;
  fSlotID = source.fSlotID;
  fEventWords = source.fEventWords;
  fACQBits = source.fACQBits;
  fLBit = source.fLBit;
  fEBit = source.fEBit;
  fEventCRC = source.fEventCRC;
  fEventCounter = source.fEventCounter;
  fDecoderCRC = source.fDecoderCRC;
  for (Int_t iChain = 0; iChain < N_CHAIN; iChain++)
    *fChainSummaryData[iChain] = *source.fChainSummaryData[iChain];
  return *this;
}

//_________________________________________________________________

AliTOFTRMSummaryData::~AliTOFTRMSummaryData()
{
  /* default destructor */
  for (Int_t iChain = 0; iChain < N_CHAIN; iChain++)
    delete fChainSummaryData[iChain];
}

//_________________________________________________________________

void
AliTOFTRMSummaryData::Reset()
{
  /* reset function */
  fHeader = kFALSE;
  fTrailer = kFALSE;
  fSlotID = 0;
  fEventWords = 0;
  fACQBits = 0;
  fLBit = 0;
  fEBit = 0;
  fEventCRC = 0;
  fEventCounter = 0;
  fDecoderCRC = 0;
  for (Int_t iChain = 0; iChain < N_CHAIN; iChain++)
    fChainSummaryData[iChain]->Reset();
}

