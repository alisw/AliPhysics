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
				   fTimeBestA(99999),
				   fTimeBestC(99999),
                                   fMultC(0),fMultA(0),
                                   fT0clock(9999999),
				   fT0trig(0),
				   fPileup(kFALSE),
				   fSattelite(kFALSE),
				   fTime1stA(99999),
				   fTime1stC(99999)
    

{
  //ctor
  for (Int_t i=0; i<24; i++) { fTime[i]=0; fADC[i]=0; fADCLED[i]=0;
    for(Int_t iHit=0; iHit<5; iHit++) {
      fTimeFull[i][iHit] = 0;   
      fOrA[iHit] = 0; 
      fOrC[iHit] = 0;  
      fTVDC[iHit] = 0; 
    }
 }
}
//_____________________________________________________________________________

AliT0RecPoint::AliT0RecPoint(const AliT0RecPoint &r):TObject(),
 						     fTimeAverage(r.fTimeAverage),
						     fTimeOnlineMean(r.fTimeOnlineMean),
						     fVertexPosition(r.fVertexPosition),
						     fTimeBestA(r.fTimeBestA),
						     fTimeBestC(r.fTimeBestC),
						     fMultC(r.fMultC),
						     fMultA(r.fMultA),
                                                     fT0clock(r.fT0clock),
						     fT0trig(r.fT0trig),
						     fPileup(r.fPileup),
						     fSattelite(r.fSattelite),
						     fTime1stA(r.fTime1stA),
						     fTime1stC(r.fTime1stC)

{
  //
  // AliT0RecPoint copy constructor
  //
  for (Int_t i=0; i<24; i++) {
    fTime[i] = r. fTime[i];
    fADC[i] = r.fADC[i]; 
    fADCLED[i] = r. fADCLED[i];
    for(Int_t iHit=0; iHit<5; iHit++) {
      fTimeFull[i][iHit] = r.fTimeFull[i][iHit];   
      fOrA[iHit] = r.fOrA[iHit]; 
      fOrC[iHit] = r.fOrC[iHit];  
      fTVDC[iHit] = r.fTVDC[iHit]; 
    }
  }
  //  ((AliT0RecPoint &) r).Copy(*this);

}
//_____________________________________________________________________________

/*
//_____________________________________________________________________________

AliT0RecPoint& AliT0RecPoint:: operator=(const AliT0RecPoint &r)
{
  //
  // assign. operator
  //

  if (this == &r)
    return *this;
  
  fTimeAverage = r.fTimeAverage;
  fTimeOnlineMean = r.fTimeOnlineMean;
  fVertexPosition = r.fVertexPosition;
  fTimeBestA =  r.fTimeBestA;
  fTimeBestC = r.fTimeBestC;
  fMultC = r.fMultC;
  fMultA = r.fMultA;
  fT0clock = r.fT0clock;
  fT0trig = r.fT0trig;
  fPileup = r.fPileup;
  fSattelite = r.fSattelite;
  fTime1stA = r.fTime1stA;
  fTime1stC = r.fTime1stC;
  for (Int_t i=0; i<24; i++) {
    fTime[i] = r. fTime[i];
    fADC[i] = r.fADC[i]; 
    fADCLED[i] = r. fADCLED[i];
    for(Int_t iHit=0; iHit<5; iHit++) {
      fTimeFull[i][iHit] = r.fTimeFull[i][iHit];   
      fOrA[iHit] = r.fOrA[iHit]; 
      fOrC[iHit] = r.fOrC[iHit];  
      fTVDC[iHit] = r.fTVDC[iHit]; 
    }
  }
  
  return *this;
}
*/
//_____________________________________________________________________________
void AliT0RecPoint::SetT0Trig(Bool_t *tr)
{
  fT0trig=0;
  for (Int_t i=0; i<5; i++) fT0trig = fT0trig | (tr[i]?(1<<i):0);
}
//_____________________________________________________________________________

void AliT0RecPoint::PrintTriggerSignals(Int_t trig)
{
  Bool_t tr[5];
  for (Int_t i=0; i<5; i++) tr[i] = (trig&(1<<i))!=0;

  AliInfo(Form("T0 triggers tvdc %d orA %d orC %d centr %d semicentral %d",tr[0],tr[1],tr[2],tr[3],tr[4]));
}
