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

/* $Id: AliFITRecPointUp.cxx 52436 2011-10-31 14:29:49Z alla $ */
/////////////////////////////////////////////////////////////////////////
//  Class AliFITRecPointUp for FIT time and ADC signals
//  fTimeA  - A side TOF signal
//  fTimeC  - C side TOF signal
//  fTimeBestA - TOF first particle on the A side
//  TimeBestC - TOF first particle on the C side
//  fTimeAverage = (fTimeBestA + TimeBestLeft ) /2. FIT signal
//  fVertex - vertex position 
//
///////////////////////////////////////////////////////////////////////



 
#include "AliFITRecPoint.h"
#include "AliLog.h"


ClassImp(AliFITRecPoint)

//------------------------------------
  AliFITRecPoint::AliFITRecPoint() : TObject()

  //ctor
  for (Int_t i=0; i<160; i++) { fTime[i]=0; fADCQTC[i]=0; }
}
//_____________________________________________________________________________

AliFITRecPoint::AliFITRecPointUp(const AliFITRecPointUp &r):TObject(){
  //
  // AliFITRecPoint copy constructor
  //
  for (Int_t i=0; i<160; i++) {
    fTime[i] = r. fTime[i];
    fADCQTC[i] = r.fADCQTC[i]; 
  }
  //  ((AliFITRecPointUp &) r).Copy(*this);

}
//_____________________________________________________________________________

/*
//_____________________________________________________________________________

AliFITRecPointUp& AliFITRecPointUp:: operator=(const AliFITRecPointUp &r)
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
  fFITclock = r.fFITclock;
  fFITtrig = r.fFITtrig;
  fPileup = r.fPileup;
  fSattelite = r.fSattelite;
  fTime1stA = r.fTime1stA;
  fTime1stC = r.fTime1stC;
  for (Int_t i=0; i<160 i++) {
    fTime[i] = r. fTime[i];
    fADC[i] = r.fADC[i]; 
  }
  
  return *this;
}

//_____________________________________________________________________________
void AliFITRecPoint::SetFITTrig(Bool_t *tr)
{
  fFITtrig=0;
  for (Int_t i=0; i<5; i++) fFITtrig = fFITtrig | (tr[i]?(1<<i):0);
}
//_____________________________________________________________________________

void AliFITRecPoint::PrintTriggerSignals(Int_t trig)
{
  Bool_t tr[5];
  for (Int_t i=0; i<5; i++) tr[i] = (trig&(1<<i))!=0;

  AliInfo(Form("FIT triggers tvdc %d orA %d orC %d centr %d semicentral %d",tr[0],tr[1],tr[2],tr[3],tr[4]));
}
*/
