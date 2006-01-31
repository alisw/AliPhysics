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

// $Id$

#include "AliMUONCalibParam.h"

#include "Riostream.h"

ClassImp(AliMUONCalibParam)

//_____________________________________________________________________________
AliMUONCalibParam::AliMUONCalibParam(Float_t mean, Float_t sigma) 
: TObject()
{
  Set(mean,sigma);
}

//_____________________________________________________________________________
AliMUONCalibParam::~AliMUONCalibParam()
{
}

//_____________________________________________________________________________
Float_t 
AliMUONCalibParam::Mean() const
{
  return fMean;
}

//_____________________________________________________________________________
void
AliMUONCalibParam::Print(Option_t*) const
{
  cout << "Mean=" << Mean() << " Sigma=" << Sigma() << endl;
}

//_____________________________________________________________________________
void 
AliMUONCalibParam::Set(Float_t mean, Float_t sigma)
{
  fMean = mean;
  fSigma = sigma;
}

//_____________________________________________________________________________
Float_t
AliMUONCalibParam::Sigma() const
{
  return fSigma;
}



