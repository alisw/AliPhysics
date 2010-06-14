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
//        This classes provide decoder summaries for data.          //
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////
                               
#include "AliTOFDecoderSummaryData.h"

ClassImp(AliTOFDecoderSummaryData)

AliTOFDecoderSummaryData::AliTOFDecoderSummaryData() :
  TObject(),
  fRunNumber(0),
  fEventNumber(0),
  fEquipmentID(0),
  fInputWords(0),
  fDecodedWords(0),
  fDecoderStatus(0),
  fErrorDetected(kFALSE),
  fErrorSlotID(0),
  fCurrentDRMID(0),
  fCurrentSlotID(0),
  fCurrentChain(0),
  fV2718Patch(kFALSE),
  fRecoverError(kFALSE),
  fRecoveringError(kFALSE),
  fSpider(kFALSE),
  fDRMSummaryData(0x0)
{
  /* default constructor */
  fDRMSummaryData = new AliTOFDRMSummaryData();
}

//_________________________________________________________________

AliTOFDecoderSummaryData::AliTOFDecoderSummaryData(const AliTOFDecoderSummaryData &source) :
  TObject(),
  fRunNumber(source.fRunNumber),
  fEventNumber(source.fEventNumber),
  fEquipmentID(source.fEquipmentID),
  fInputWords(source.fInputWords),
  fDecodedWords(source.fDecodedWords),
  fDecoderStatus(source.fDecoderStatus),
  fErrorDetected(source.fErrorDetected),
  fErrorSlotID(source.fErrorSlotID),
  fCurrentDRMID(source.fCurrentDRMID),
  fCurrentSlotID(source.fCurrentSlotID),
  fCurrentChain(source.fCurrentChain),
  fV2718Patch(source.fV2718Patch),
  fRecoverError(source.fRecoverError),
  fRecoveringError(source.fRecoveringError),
  fSpider(kFALSE),
  fDRMSummaryData(0x0)
{
  /* copy constructor */
  fDRMSummaryData = new AliTOFDRMSummaryData(*source.fDRMSummaryData);
}

//_________________________________________________________________

AliTOFDecoderSummaryData &
AliTOFDecoderSummaryData::operator = (const AliTOFDecoderSummaryData &source)
{
  /* operator = */
  fRunNumber = source.fRunNumber;
  fEventNumber = source.fEventNumber;
  fEquipmentID = source.fEquipmentID;
  fInputWords = source.fInputWords;
  fDecodedWords = source.fDecodedWords;
  fDecoderStatus = source.fDecoderStatus;
  fErrorDetected = source.fErrorDetected;
  fErrorSlotID = source.fErrorSlotID;
  fCurrentDRMID = source.fCurrentDRMID;
  fCurrentSlotID = source.fCurrentSlotID;
  fCurrentChain = source.fCurrentChain;
  fV2718Patch = source.fV2718Patch;
  fRecoverError = source.fRecoverError;
  fRecoveringError = source.fRecoveringError;
  fSpider = source.fSpider;
  *fDRMSummaryData = *source.fDRMSummaryData;
  return *this;
}

//_________________________________________________________________

AliTOFDecoderSummaryData::~AliTOFDecoderSummaryData()
{
  /* default destructor */
  delete fDRMSummaryData;
}

//_________________________________________________________________

void 
AliTOFDecoderSummaryData::Reset()
{
  /* reset function */
  fRunNumber = 0;
  fEventNumber = 0;
  fEquipmentID = 0;
  fInputWords = 0;
  fDecodedWords = 0;
  fDecoderStatus = 0;
  fErrorDetected = kFALSE;
  fErrorSlotID = 0;
  fCurrentDRMID = 0;
  fCurrentSlotID = 0;
  fCurrentChain = 0;
  fV2718Patch = kFALSE;
  fRecoverError = kFALSE;
  fRecoveringError = kFALSE;
  fSpider = kFALSE;
  fDRMSummaryData->Reset();
}

