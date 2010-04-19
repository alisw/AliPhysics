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

/* $Id$ */
/////////////////////////////////////////////////////////////////////////
//  Class AliT0RecPoint for T0 time and ADC signals
//  fTimeA  - A side TOF signal
//  fTimeC  - C side TOF signal
//  fTimeBestA - TOF first particle on the A side
//  TimeBestC - TOF first particle on the C side
//  fTimeAverage = (fTimeBestA + TimeBestLeft ) /2. T0 signal
//  fVertex - vertex position 
//
///////////////////////////////////////////////////////////////////////



 
#include "AliT0RecPoint.h"
#include "AliLog.h"


ClassImp(AliT0RecPoint)

//------------------------------------
  AliT0RecPoint::AliT0RecPoint() : TObject(),
				   fTimeAverage(99999),
				   fTimeOnlineMean(99999),
				   fVertexPosition(999999),
				   fTimeBestA(0),fTimeBestC(0),
                                   fMultC(0),fMultA(0),
                                   fT0clock(9999999),
				   fT0trig(0)

{
  //ctor
  // fTimeAverage=99999;
  fTimeBestA=99999;
  fTimeBestC=99999;
  //  fVertexPosition=99999;
  fMultA=0;
  fMultC=0;
  for (Int_t i=0; i<24; i++) { fTime[i]=0; fADC[i]=0; fADCLED[i]=0;}
}
//_____________________________________________________________________________

AliT0RecPoint::AliT0RecPoint(const AliT0RecPoint &r):TObject(),
 						     fTimeAverage(999999),
						     fTimeOnlineMean(999999),
						     fVertexPosition(999999),
						     fTimeBestA(0),fTimeBestC(0),
						     fMultC(0),fMultA(0),
                                                     fT0clock(9999999),
						     fT0trig(0)
{
  //
  // AliT0RecPoint copy constructor
  //

  ((AliT0RecPoint &) r).Copy(*this);

}
//_____________________________________________________________________________

void AliT0RecPoint::SetT0Trig(Bool_t *tr)
{
  fT0trig=0;
  for (Int_t i=0; i<5; i++) fT0trig=fT0trig<<1|tr[i];
}
//_____________________________________________________________________________

void AliT0RecPoint::PrintTriggerSignals(Int_t trig)
{
  Bool_t tr[5];
  for (Int_t i=0; i<5; i++) tr[i]=trig&(1<<i);

  AliInfo(Form("T0 triggers %d %d %d %d %d",tr[0],tr[1],tr[2],tr[3],tr[4]));
}
