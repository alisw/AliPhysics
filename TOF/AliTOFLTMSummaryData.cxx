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
//        This class provides a summary for LTM data.               //
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include "AliTOFLTMSummaryData.h"

ClassImp(AliTOFLTMSummaryData)

AliTOFLTMSummaryData::AliTOFLTMSummaryData() :
  TObject(),
  fHeader(kFALSE),
  fTrailer(kFALSE),
  fSlotID(0),
  fEventWords(0),
  fCBit(0),
  fFault(0),
  fEventCRC(0),
  fEventNumber(0),
  fDecoderCRC(0)
{
  /* default constructor */
  for (Int_t iPDL = 0; iPDL < LTM_N_PDL; iPDL++)
    fPDL[iPDL] = 0;
  for (Int_t iADC = 0; iADC < LTM_N_ADC; iADC++)
    fADC[iADC] = 0;
  for (Int_t iOR = 0; iOR < LTM_N_OR; iOR++)
    fOR[iOR] = 0;
}

//_________________________________________________________________

AliTOFLTMSummaryData::AliTOFLTMSummaryData(const AliTOFLTMSummaryData &source) :
  TObject(),
  fHeader(source.fHeader),
  fTrailer(source.fTrailer),
  fSlotID(source.fSlotID),
  fEventWords(source.fEventWords),
  fCBit(source.fCBit),
  fFault(source.fFault),
  fEventCRC(source.fEventCRC),
  fEventNumber(source.fEventNumber),
  fDecoderCRC(source.fDecoderCRC)
{
  /* copy constructor */
  for (Int_t iPDL = 0; iPDL < LTM_N_PDL; iPDL++)
    fPDL[iPDL] = source.fPDL[iPDL];
  for (Int_t iADC = 0; iADC < LTM_N_ADC; iADC++)
    fADC[iADC] = source.fADC[iADC];
  for (Int_t iOR = 0; iOR < LTM_N_OR; iOR++)
    fOR[iOR] = source.fOR[iOR];
}

//_________________________________________________________________

AliTOFLTMSummaryData &
AliTOFLTMSummaryData::operator = (const AliTOFLTMSummaryData &source)
{
  /* operator = */
  fHeader = source.fHeader;
  fTrailer = source.fTrailer;
  fSlotID = source.fSlotID;
  fEventWords = source.fEventWords;
  fCBit = source.fCBit;
  fFault = source.fFault;
  for (Int_t iPDL = 0; iPDL < LTM_N_PDL; iPDL++)
    fPDL[iPDL] = source.fPDL[iPDL];
  for (Int_t iADC = 0; iADC < LTM_N_ADC; iADC++)
    fADC[iADC] = source.fADC[iADC];
  for (Int_t iOR = 0; iOR < LTM_N_OR; iOR++)
    fOR[iOR] = source.fOR[iOR];
  fEventCRC = source.fEventCRC;
  fEventNumber = source.fEventNumber;
  fDecoderCRC = source.fDecoderCRC;
  return *this;
}

//_________________________________________________________________

AliTOFLTMSummaryData::~AliTOFLTMSummaryData()
{
  /* default destructor */
}

//_________________________________________________________________

void
AliTOFLTMSummaryData::Reset()
{
  /* reset function */
  fHeader = kFALSE;
  fTrailer = kFALSE;
  fSlotID = 0;
  fEventWords = 0;
  fCBit = 0;
  fFault = 0;
  for (Int_t iPDL = 0; iPDL < LTM_N_PDL; iPDL++)
    fPDL[iPDL] = 0;
  for (Int_t iADC = 0; iADC < LTM_N_ADC; iADC++)
    fADC[iADC] = 0;
  for (Int_t iOR = 0; iOR < LTM_N_OR; iOR++)
    fOR[iOR] = 0;
  fEventCRC = 0;
  fEventNumber = 0;
  fDecoderCRC = 0;
}

