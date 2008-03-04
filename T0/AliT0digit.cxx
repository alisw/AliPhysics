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
//  fTimeA  - right side TOF signal
//  fTimeC  - left side TOF signal
//  fTimeBestA - TOF first particle on the right side
//  TimeBestC - TOF first particle on the left side
//  fTimeAverage = (fTimeBestA + TimeBestC ) /2. T0 signal
//  fTimeDiff = fTimeBestA - TimeBestC  
//
///////////////////////////////////////////////////////////////////////

#include "AliT0digit.h"
#include "TArrayI.h"

ClassImp(AliT0digit)

//-----------------------------------------------
  AliT0digit::AliT0digit() :TObject(),
			    fTimeCFD(new TArrayI(24)),    
			    fQT0( new TArrayI(24)),     
			    fTimeLED( new TArrayI(24)), 
			    fQT1( new TArrayI(24)),     
			    fTimeAverage(99999),
			    fTimeDiff(99999),
			    fBestTimeA(99999),
			    fBestTimeC (99999),
			    fSumMult(0),
			    fRefPoint(0)

{
  //
}

//_____________________________________________________________________________

AliT0digit::~AliT0digit() {
  // destructor
  delete fTimeCFD;
  delete fQT0;
  delete fTimeLED;
  delete fQT1;
}
//-----------------------------------
void AliT0digit::SetTimeCFD (TArrayI &o)
{
  ////////////////////////////////////////
  if(fTimeCFD)delete  fTimeCFD;
  fTimeCFD = new TArrayI(24);

  Int_t i;
  for (i=0; i<24; i++)
    {
      Int_t buf=o.At(i);
      fTimeCFD->AddAt(buf,i);
    }
}
//--------------------------------------------
void AliT0digit::GetTimeCFD (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fTimeCFD->At(i);
    }
}
//--------------------------------------------
void AliT0digit::GetQT0 (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fQT0->At(i);
    }
}
//--------------------------------------------
void AliT0digit::SetQT0 (TArrayI &o)
{
  //
  if(fQT0)delete fQT0;
  fQT0  = new TArrayI(24);
  Int_t i;
  for (i=0; i<24; i++)
    {
      Int_t buf=(o.At(i));
      fQT0->AddAt(buf,i);
    }
}
//-----------------------------------
void AliT0digit::SetTimeLED (TArrayI &o)
{
  ////////////////////////////////////////
  if(fTimeLED)delete fTimeLED;
  fTimeLED = new TArrayI(24);

  Int_t i;
  for (i=0; i<24; i++)
    {
      Int_t buf=o.At(i);
      fTimeLED->AddAt(buf,i);
    }
}
//--------------------------------------------
void AliT0digit::GetTimeLED (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fTimeLED->At(i);
    }
}
//--------------------------------------------
void AliT0digit::GetQT1 (TArrayI &o)
{
  //
  Int_t i;
  for (i=0; i<24; i++)
    {
      o[i]=fQT1->At(i);
    }
}
//--------------------------------------------
void AliT0digit::SetQT1 (TArrayI &o)
{
  //
  if(fQT1)delete fQT1;
  fQT1  = new TArrayI(24);
  Int_t i;
  for (i=0; i<24; i++)
    {
      Int_t buf=(o.At(i));
      fQT1->AddAt(buf,i);
    }
}
