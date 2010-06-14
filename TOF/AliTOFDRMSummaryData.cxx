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
//        This class provides a summary for DRM data.               //
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include "AliTOFDRMSummaryData.h"

ClassImp(AliTOFDRMSummaryData)

AliTOFDRMSummaryData::AliTOFDRMSummaryData() :
  TObject(),
  fHeader(kFALSE),
  fTrailer(kFALSE),
  fSlotID(0),
  fEventWords(0),
  fDRMID(0),
  fLocalEventCounter(0),
  fPartecipatingSlotID(0),
  fCBit(0),
  fVersID(0),
  fDRMhSize(0),
  fSlotEnableMask(0),
  fFaultID(0),
  fRTOBit(0),
  fL0BCID(0),
  fRunTimeInfo(0),
  fTemperature(0),
  fACKBit(0),
  fSensAD(0),
  fEventCRC(0),
  fDecoderCRC(0),
  fDecoderSlotEnableMask(0),
  fLTMSummaryData(0x0)
{
  /* default constructor */
  fLTMSummaryData = new AliTOFLTMSummaryData();
  for (Int_t iTRM = 0; iTRM < N_TRM; iTRM++)
    fTRMSummaryData[iTRM] = new AliTOFTRMSummaryData();
}

//_________________________________________________________________

AliTOFDRMSummaryData::AliTOFDRMSummaryData(const AliTOFDRMSummaryData &source) :
  TObject(),
  fHeader(source.fHeader),
  fTrailer(source.fTrailer),
  fSlotID(source.fSlotID),
  fEventWords(source.fEventWords),
  fDRMID(source.fDRMID),
  fLocalEventCounter(source.fLocalEventCounter),
  fPartecipatingSlotID(source.fPartecipatingSlotID),
  fCBit(source.fCBit),
  fVersID(source.fVersID),
  fDRMhSize(source.fDRMhSize),
  fSlotEnableMask(source.fSlotEnableMask),
  fFaultID(source.fFaultID),
  fRTOBit(source.fRTOBit),
  fL0BCID(source.fL0BCID),
  fRunTimeInfo(source.fRunTimeInfo),
  fTemperature(source.fTemperature),
  fACKBit(source.fACKBit),
  fSensAD(source.fSensAD),
  fEventCRC(source.fEventCRC),
  fDecoderCRC(source.fDecoderCRC),
  fDecoderSlotEnableMask(source.fDecoderSlotEnableMask),
  fLTMSummaryData(0x0)
{
  /* copy constructor */
  fLTMSummaryData = new AliTOFLTMSummaryData(*source.fLTMSummaryData);
  for (Int_t iTRM = 0; iTRM < N_TRM; iTRM++)
    fTRMSummaryData[iTRM] = new AliTOFTRMSummaryData(*source.fTRMSummaryData[iTRM]);
}

//_________________________________________________________________

AliTOFDRMSummaryData &
AliTOFDRMSummaryData::operator = (const AliTOFDRMSummaryData &source)
{
  /* operator = */
  fHeader = source.fHeader;
  fTrailer = source.fTrailer;
  fSlotID = source.fSlotID;
  fEventWords = source.fEventWords;
  fDRMID = source.fDRMID;
  fLocalEventCounter = source.fLocalEventCounter;
  fPartecipatingSlotID = source.fPartecipatingSlotID;
  fCBit = source.fCBit;
  fVersID = source.fVersID;
  fDRMhSize = source.fDRMhSize;
  fSlotEnableMask = source.fSlotEnableMask;
  fFaultID = source.fFaultID;
  fRTOBit = source.fRTOBit;
  fL0BCID = source.fL0BCID;
  fRunTimeInfo = source.fRunTimeInfo;
  fTemperature = source.fTemperature;
  fACKBit = source.fACKBit;
  fSensAD = source.fSensAD;
  fEventCRC = source.fEventCRC;
  fDecoderCRC = source.fDecoderCRC;
  fDecoderSlotEnableMask = source.fDecoderSlotEnableMask;
  *fLTMSummaryData = *source.fLTMSummaryData;
  for (Int_t iTRM = 0; iTRM < N_TRM; iTRM++)
    *fTRMSummaryData[iTRM] = *source.fTRMSummaryData[iTRM];
  return *this;
}

//_________________________________________________________________

AliTOFDRMSummaryData::~AliTOFDRMSummaryData()
{
  /* default destructor */
  delete fLTMSummaryData;
  for (Int_t iTRM = 0; iTRM < N_TRM; iTRM++)
    delete fTRMSummaryData[iTRM];
}

//_________________________________________________________________

void
AliTOFDRMSummaryData::Reset()
{
  /* reset function */
  fHeader = kFALSE;
  fTrailer = kFALSE;
  fSlotID = 0;
  fEventWords = 0;
  fDRMID = 0;
  fLocalEventCounter = 0;
  fPartecipatingSlotID = 0;
  fCBit = 0;
  fVersID = 0;
  fDRMhSize = 0;
  fSlotEnableMask = 0;
  fFaultID = 0;
  fRTOBit = 0;
  fL0BCID = 0;
  fRunTimeInfo = 0;
  fTemperature = 0;
  fACKBit = 0;
  fSensAD = 0;
  fEventCRC = 0;
  fDecoderCRC = 0;
  fDecoderSlotEnableMask = 0;
  fLTMSummaryData->Reset();
  for (Int_t iTRM = 0; iTRM < N_TRM; iTRM++)
    fTRMSummaryData[iTRM]->Reset();
}
