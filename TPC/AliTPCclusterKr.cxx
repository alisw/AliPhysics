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

/* $Id: AliTPCclusterKr.cxx,v 1.7 2008/01/22 17:24:53 matyja Exp $ */

//-----------------------------------------------------------------
//           Implementation of the TPC Kr cluster class
//
// Origin: Adam Matyja, INP PAN, adam.matyja@ifj.edu.pl
//-----------------------------------------------------------------

#include "AliTPCclusterKr.h"
#include "TObject.h"
#include <vector>
#include "AliTPCvtpr.h"

ClassImp(AliTPCclusterKr)


AliTPCclusterKr::AliTPCclusterKr()
:fMax(),
  fADCcluster(0),
  fNsec(0),
  fNpads(0),
  fSize(0)
{
//
// default constructor
//
}

AliTPCclusterKr::AliTPCclusterKr(const AliTPCclusterKr &param)
:fMax(),
  fADCcluster(0),
  fNsec(0),
  fNpads(0),
  fSize(0)
{
//
// copy constructor
//
  fADCcluster = param.fADCcluster;
  fNsec  = param.fNsec ;
  fNpads = param.fNpads;
  fMax = param.fMax;
  fCluster=param.fCluster;
  fSize=param.fSize;
} 

AliTPCclusterKr &AliTPCclusterKr::operator = (const AliTPCclusterKr & param)
{
  fADCcluster = param.fADCcluster;
  fNsec  = param.fNsec ;
  fNpads = param.fNpads;
  fMax = param.fMax;
  fCluster=param.fCluster;
  fSize=param.fSize;
  return (*this);
}

AliTPCclusterKr::~AliTPCclusterKr()
{
  //
  // destructor
  //
}
