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
//        This class provides a definition for TDC hits.            //
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include "AliTOFTDCHit.h"
#define TIME_BIN_WIDTH          24.4e-3//ns
#define TOT_BIN_WIDTH           48.8e-3//ns
#define TIME_TO_TOT_BIN_WIDTH   ( TIME_BIN_WIDTH / TOT_BIN_WIDTH )
#define TOT_TO_TIME_BIN_WIDTH   ( TOT_BIN_WIDTH / TIME_BIN_WIDTH )

ClassImp(AliTOFTDCHit)

AliTOFTDCHit::AliTOFTDCHit() :
  TObject(),
  fHitTime(0),
  fTOTWidth(0),
  fChan(0),
  fTDCID(0),
  fEBit(0),
  fPSBits(0)
{
  /* default constructor */
}

//_________________________________________________________________

AliTOFTDCHit::AliTOFTDCHit(const AliTOFTDCHit &source) :
  TObject(),
  fHitTime(source.fHitTime),
  fTOTWidth(source.fTOTWidth),
  fChan(source.fChan),
  fTDCID(source.fTDCID),
  fEBit(source.fEBit),
  fPSBits(source.fPSBits)
{
  /* copy constructor */
}

//_________________________________________________________________

AliTOFTDCHit &
AliTOFTDCHit::operator = (const AliTOFTDCHit &source)
{
  /* operator = */
  fHitTime = source.fHitTime;
  fTOTWidth = source.fTOTWidth;
  fChan = source.fChan;
  fTDCID = source.fTDCID;
  fEBit = source.fEBit;
  fPSBits = source.fPSBits;
  return *this;
}

#if 0
//_________________________________________________________________

AliTOFTDCHit &
AliTOFTDCHit::operator - (const AliTOFTDCHit &source)
{
  /* operator - */
  fHitTime = fHitTime - source.fHitTime;
  return *this;
}
#endif

//_________________________________________________________________

AliTOFTDCHit &
AliTOFTDCHit::operator -= (const AliTOFTDCHit &source)
{
  /* operator -= */
  fHitTime -= source.fHitTime;
  return *this;
}

//_________________________________________________________________

AliTOFTDCHit &
AliTOFTDCHit::operator << (const AliTOFTDCHit &source)
{
  /* operator << */
  /* build packed hit */
  fTOTWidth = source.fHitTime - fHitTime; /* compute TOT width */
  fTOTWidth = (UShort_t)(fTOTWidth * TIME_TO_TOT_BIN_WIDTH); /* convert into 48.8 ps bins */
  fEBit = fEBit | source.fEBit; /* set E bit as or */
  fPSBits = 0; /* set PB bits as packed hit */
  return *this;
}

//_________________________________________________________________

AliTOFTDCHit::~AliTOFTDCHit()
{
  /* destructor */
}
