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

/* $Id: AliTPCcluster.cxx,v 1.7 2008/01/22 17:24:53 matyja Exp $ */

//-----------------------------------------------------------------
//           TPC cordinate Class
//
//  Origin: Adam Matyja, INP PAN, adam.matyja@ifj.edu.pl
//-----------------------------------------------------------------

#include "AliTPCvtpr.h"

ClassImp(AliTPCvtpr)

AliTPCvtpr::AliTPCvtpr()
:fAdc(0),
  fTime(0),
  fPad(0),
  fRow(0)
{
//
// constructor
//
}

AliTPCvtpr::AliTPCvtpr(Short_t max,Short_t nt,Short_t np,Short_t nr)
:fAdc(0),
  fTime(0),
  fPad(0),
  fRow(0)
{
//
// another constructor
//
  fAdc=max;
  fTime=nt;
  fPad=np;
  fRow=nr;
}

AliTPCvtpr::AliTPCvtpr(const AliTPCvtpr & param)
:fAdc(0),
  fTime(0),
  fPad(0),
  fRow(0)
{
//
// copy constructor
//
  fAdc  = param.fAdc;
  fTime = param.fTime;
  fPad  = param.fPad;
  fRow  = param.fRow;
}

AliTPCvtpr::~AliTPCvtpr()
{
  //
  // destructor
  //
}

AliTPCvtpr & AliTPCvtpr::operator = (const AliTPCvtpr & param)
{
  fAdc  = param.fAdc;
  fTime = param.fTime;
  fPad  = param.fPad;
  fRow  = param.fRow;
  return (*this);
}
