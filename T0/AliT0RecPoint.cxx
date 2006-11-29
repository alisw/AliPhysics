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
//  fTimeRight  - right side TOF signal
//  fTimeLeft  - left side TOF signal
//  fTimeBestRight - TOF first particle on the right side
//  TimeBestLeft - TOF first particle on the left side
//  fTimeAverage = (fTimeBestRight + TimeBestLeft ) /2. T0 signal
//  fVertex - vertex position 
//
///////////////////////////////////////////////////////////////////////



 
#include "AliT0RecPoint.h"
#include <Riostream.h>

ClassImp(AliT0RecPoint)

//------------------------------------
  AliT0RecPoint::AliT0RecPoint() : TObject(), fTimeAverage(0),fVertexPosition(0),fTimeBestRight(0),fTimeBestLeft(0),fMultC(0),fMultA(0)
{
  //ctor
  fTimeAverage=99999;
  fTimeBestRight=99999;
  fTimeBestLeft=99999;
  fVertexPosition=99999;
  fMultA=0;
  fMultC=0;
  for (Int_t i=0; i<24; i++) { fTime[i]=0; fADC[i]=0; fADCLED[i]=0;}
}
