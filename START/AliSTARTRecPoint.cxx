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
//  Class AliSTARTRecPoint for START time and ADC signals
//  fTimeRight  - right side TOF signal
//  fTimeLeft  - left side TOF signal
//  fTimeBestRight - TOF first particle on the right side
//  TimeBestLeft - TOF first particle on the left side
//  fTimeAverage = (fTimeBestRight + TimeBestLeft ) /2. START signal
//  fVertex - vertex position 
//
///////////////////////////////////////////////////////////////////////



 
#include <TArrayI.h>
#include "AliSTARTRecPoint.h"
#include <Riostream.h>

ClassImp(AliSTARTRecPoint)

//------------------------------------
 AliSTARTRecPoint::AliSTARTRecPoint() : TObject()
{
  //ctor
  fTimeAverage=9999;
  fTimeBestRight=9999;
  fTimeBestLeft=9999;

  fTime = new TArrayI(24);  
  fADC  = new TArrayI(24);  
}
//-----------------------------------
AliSTARTRecPoint::~AliSTARTRecPoint() {
  // destructor
  delete fTime;
  delete fADC;
}
//-----------------------------------
void AliSTARTRecPoint::SetTime (TArrayI &o)
{
  ////////////////////////////////////////

  Int_t i;
  for (i=0; i<24; i++)
    {
      Int_t buf=o.At(i);
      fTime->AddAt(buf,i);
    }
}
//--------------------------------------------
void AliSTARTRecPoint::GetTime (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fTime->At(i);
    }
}
//--------------------------------------------
void AliSTARTRecPoint::GetADC (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fADC->At(i);
    }
}
//--------------------------------------------
void AliSTARTRecPoint::SetADC (TArrayI &o)
{
  //
  Int_t i;
  //  Float_t fProcessKoef=1; // for pb 0.001
  for (i=0; i<24; i++)
    {
      Int_t buf=(o.At(i));
      fADC->AddAt(buf,i);
    }
}
