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
//  Class AliSTARTdigit for START digits
//  fTimeRight  - right side TOF signal
//  fTimeLeft  - left side TOF signal
//  fTimeBestRight - TOF first particle on the right side
//  TimeBestLeft - TOF first particle on the left side
//  fTimeAverage = (fTimeBestRight + TimeBestLeft ) /2. START signal
//  fTimeDiff = fTimeBestRight - TimeBestLeft  
//
///////////////////////////////////////////////////////////////////////

#include "AliSTARTdigit.h"
#include "TArrayI.h"

ClassImp(AliSTARTdigit)

//-----------------------------------------------
  AliSTARTdigit::AliSTARTdigit() :TObject()
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
AliSTARTdigit::~AliSTARTdigit() {
  // destructor
  delete fTime;
  delete fADC;
  delete fTimeAmp;
  delete fADCAmp;
}
//-----------------------------------
void AliSTARTdigit::SetTime (TArrayI &o)
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
void AliSTARTdigit::GetTime (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fTime->At(i);
    }
}
//--------------------------------------------
void AliSTARTdigit::GetADC (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fADC->At(i);
    }
}
//--------------------------------------------
void AliSTARTdigit::SetADC (TArrayI &o)
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
void AliSTARTdigit::SetTimeAmp (TArrayI &o)
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
void AliSTARTdigit::GetTimeAmp (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fTimeAmp->At(i);
    }
}
//--------------------------------------------
void AliSTARTdigit::GetADCAmp (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fADCAmp->At(i);
    }
}
//--------------------------------------------
void AliSTARTdigit::SetADCAmp (TArrayI &o)
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
