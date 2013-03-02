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

/* $Id: AliTPCcluster.cxx,v 1.7 2007/03/15 17:24:53 matyja Exp $ */

//-----------------------------------------------------------------
//           Implementation of the TPC cluster class
//
//    Origin: Adam Matyja, INP PAN, adam.matyja@ifj.edu.pl
//-----------------------------------------------------------------

#include "AliPadMax.h"
#include "AliTPCvtpr.h"

ClassImp(AliPadMax)
AliPadMax::AliPadMax()
:AliTPCvtpr(),
  fBegin(0),
  fEnd(0),
  fSumAdc(0)
{
//
// constructor
//
}


AliPadMax::AliPadMax(AliTPCvtpr vtpr,
		     Short_t beg,Short_t end,Short_t sum)
:AliTPCvtpr(),
  fBegin(0),
  fEnd(0),
  fSumAdc(0)
{
  fAdc  =  vtpr.GetAdc();
  fTime =  vtpr.GetTime();
  fPad  =  vtpr.GetPad();
  fRow  =  vtpr.GetRow();    
  fX  =  vtpr.GetX();  
  fY  =  vtpr.GetY();  
  fT  =  vtpr.GetT();  

  fBegin=beg;
  fEnd=end;
  fSumAdc=sum;
}

AliPadMax::AliPadMax(Short_t max,Short_t nt,Short_t np,Short_t nr,
		     Double_t x,Double_t y,Double_t t,
		     Short_t beg,Short_t end,Short_t sum)
:AliTPCvtpr(),
  fBegin(0),
  fEnd(0),
  fSumAdc(0)
{
  fAdc  =  max;
  fTime =  nt;
  fPad  =  np;
  fRow  =  nr; 
  fX  = x;
  fY  = y;
  fT  = t;

  fBegin=beg;
  fEnd=end;
  fSumAdc=sum;
}

AliPadMax::~AliPadMax()
{
  //
  // destructor
  //
}
