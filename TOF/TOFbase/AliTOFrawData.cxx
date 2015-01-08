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

/*
$Log$
Revision 1.2  2007/03/28 10:50:33  decaro
Rounding off problem in rawData coding/decoding: solved

Revision 1.1  2007/02/20 15:57:00  decaro
Raw data update: to read the TOF raw data defined in UNPACKED mode

Revision 0.1  2006/12/15 A.De Caro
   Introuction
*/

//////////////////////////////////////////////////////
//                                                  //
//  This class provides the TOF raw data object     //
//                                                  //
//////////////////////////////////////////////////////

#include "AliLog.h"

#include "AliTOFGeometry.h"
#include "AliTOFrawData.h"

ClassImp(AliTOFrawData)

//_____________________________________________________________________________
AliTOFrawData::AliTOFrawData() :
  TObject(),
  fACQflag(-1),
  fPSbit(-1),
  fTRM(-1),
  fTRMchain(-1),
  fTDC(-1),
  fTDCchannel(-1),
  fLeading(-1),
  fTrailing(-1),
  fToT(-1),
  fTime(-1),
  fError(-1),
  fDeltaBC(-1),
  fL0L1Latency(-1)
{

  // default ctr

}

//_____________________________________________________________________________
AliTOFrawData::AliTOFrawData(Int_t a, Int_t b, Int_t c, Int_t d, Int_t e,
			     Int_t f, Int_t g, Int_t h, Int_t l) :
  TObject(),
  fACQflag(h),
  fPSbit(g),
  fTRM(a),
  fTRMchain(b),
  fTDC(c),
  fTDCchannel(d),
  fLeading(-1),
  fTrailing(-1),
  fToT(f),
  fTime(e),
  fError(l),
  fDeltaBC(0),
  fL0L1Latency(0)
{

// ctr

}

//_____________________________________________________________________________
AliTOFrawData::AliTOFrawData(Int_t a, Int_t b, Int_t c, Int_t d, Int_t e,
			     Int_t f, Int_t ee, Int_t ff, Int_t g, Int_t h, Int_t l, Int_t deltaBC, Int_t l0l1) :
  TObject(),
  fACQflag(h),
  fPSbit(g),
  fTRM(a),
  fTRMchain(b),
  fTDC(c),
  fTDCchannel(d),
  fLeading(ee),
  fTrailing(ff),
  fToT(f),
  fTime(e),
  fError(l),
  fDeltaBC(deltaBC),
  fL0L1Latency(l0l1)
{

  // ctr
  if (fTime==-1) fTime = fLeading;
}

//_____________________________________________________________________________
AliTOFrawData::AliTOFrawData(const AliTOFrawData& r) :
  TObject(),
  fACQflag(-1),
  fPSbit(-1),
  fTRM(-1),
  fTRMchain(-1),
  fTDC(-1),
  fTDCchannel(-1),
  fLeading(-1),
  fTrailing(-1),
  fToT(-1),
  fTime(-1),
  fError(-1),
  fDeltaBC(-1),
  fL0L1Latency(-1)
{

  // dummy copy constructor

  fACQflag    = r.fACQflag;
  fPSbit      = r.fPSbit;
  fTRM        = r.fTRM;
  fTRMchain   = r.fTRMchain;
  fTDC        = r.fTDC;
  fTDCchannel = r.fTDCchannel;
  fLeading    = r.fLeading;
  fTrailing   = r.fTrailing;
  fToT        = r.fToT;
  fTime       = r.fTime;
  fError      = r.fError;
  fDeltaBC      = r.fDeltaBC;
  fL0L1Latency      = r.fL0L1Latency;

}

//_____________________________________________________________________________
AliTOFrawData& AliTOFrawData::operator=(const AliTOFrawData& r)
{

  // Assignment operator

  if (this == &r)
    return *this;

  TObject::operator=(r);
  fACQflag    = r.fACQflag;
  fPSbit      = r.fPSbit;
  fTRM        = r.fTRM;
  fTRMchain   = r.fTRMchain;
  fTDC        = r.fTDC;
  fTDCchannel = r.fTDCchannel;
  fLeading    = r.fLeading;
  fTrailing   = r.fTrailing;
  fToT        = r.fToT;
  fTime       = r.fTime;
  fError      = r.fError;
  fDeltaBC    = r.fDeltaBC;
  fL0L1Latency= r.fL0L1Latency;
  return *this;

}

//_____________________________________________________________________________
void AliTOFrawData::Update(Int_t tof, Int_t tot, Int_t leading, Int_t trailing, Int_t psBit, Int_t acq, Int_t errorFlag)
{
  //
  // To update a raw data object:
  //  if there is just a leading edge measurement,
  //  this method adds the trailing edge measurement
  //  to evaluate the time-of-flight and time-over-threshold measurements
  //

  AliDebug(2,Form(" %d %d %d %d %d %d",tof, tot, leading, psBit, acq, errorFlag));

  if (fLeading!=-1 /*&& fTime==-1*/ && fToT==-1 && trailing!=-1) { // adc

    fTrailing = trailing;
    fTime = fLeading;
    fToT = Int_t((trailing - fLeading)*AliTOFGeometry::TdcBinWidth()/AliTOFGeometry::ToTBinWidth());
    
  }

}

//_____________________________________________________________________________
Int_t AliTOFrawData::GetTOT() const
{
  //
  // Return the TOT measurement
  //

  Int_t dummyToT = 0;
  if (fLeading!=-1 && fToT==-1) dummyToT = 0;
  else dummyToT = fToT;

  return dummyToT;

}
