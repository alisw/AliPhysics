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
//  Class AliT0digit for T0 digits
//  fTimeRight  - right side TOF signal
//  fTimeLeft  - left side TOF signal
//  fTimeBestRight - TOF first particle on the right side
//  TimeBestLeft - TOF first particle on the left side
//  fTimeAverage = (fTimeBestRight + TimeBestLeft ) /2. T0 signal
//  fTimeDiff = fTimeBestRight - TimeBestLeft  
//
///////////////////////////////////////////////////////////////////////

#include "AliT0digit.h"
#include "TArrayI.h"

ClassImp(AliT0digit)

//-----------------------------------------------
  AliT0digit::AliT0digit() :TObject()
{

  fTimeAverage   = 99999;
  fTimeDiff      = 99999;
  fBestTimeRight = 99999;
  fBestTimeLeft  = 99999;

  fTime = new TArrayI(24);
  fADC  = new TArrayI(24);
  fTimeAmp = new TArrayI(24);
  fADCAmp  = new TArrayI(24);
}

//-----------------------------------
AliT0digit::~AliT0digit() {
  // destructor
  delete fTime;
  delete fADC;
  delete fTimeAmp;
  delete fADCAmp;
}
//-----------------------------------
void AliT0digit::SetTime (TArrayI &o)
{
  ////////////////////////////////////////
  fTime = new TArrayI(24);

  Int_t i;
  for (i=0; i<24; i++)
    {
      Int_t buf=o.At(i);
      fTime->AddAt(buf,i);
    }
}
//--------------------------------------------
void AliT0digit::GetTime (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fTime->At(i);
    }
}
//--------------------------------------------
void AliT0digit::GetADC (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fADC->At(i);
    }
}
//--------------------------------------------
void AliT0digit::SetADC (TArrayI &o)
{
  //
  fADC  = new TArrayI(24);
  Int_t i;
  for (i=0; i<24; i++)
    {
      Int_t buf=(o.At(i));
      fADC->AddAt(buf,i);
    }
}
//-----------------------------------
void AliT0digit::SetTimeAmp (TArrayI &o)
{
  ////////////////////////////////////////
  fTimeAmp = new TArrayI(24);

  Int_t i;
  for (i=0; i<24; i++)
    {
      Int_t buf=o.At(i);
      fTimeAmp->AddAt(buf,i);
    }
}
//--------------------------------------------
void AliT0digit::GetTimeAmp (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fTimeAmp->At(i);
    }
}
//--------------------------------------------
void AliT0digit::GetADCAmp (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fADCAmp->At(i);
    }
}
//--------------------------------------------
void AliT0digit::SetADCAmp (TArrayI &o)
{
  //
  fADCAmp  = new TArrayI(24);
  Int_t i;
  for (i=0; i<24; i++)
    {
      Int_t buf=(o.At(i));
      fADCAmp->AddAt(buf,i);
    }
}
