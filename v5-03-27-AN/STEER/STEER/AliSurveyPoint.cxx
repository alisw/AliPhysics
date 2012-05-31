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

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliSurveyPoint					   //
//  Retrieve and Convert survey data into ROOT Objects		   //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliSurveyPoint.h"

ClassImp(AliSurveyPoint)
  
//_____________________________________________________________________________
AliSurveyPoint::AliSurveyPoint():
  fPointName(""),
  fX(0.0),
  fY(0.0),
  fZ(0.0),
  fPrecisionX(0.0),
  fPrecisionY(0.0),
  fPrecisionZ(0.0),
  fType(' '),
  fTargetUsed(kTRUE)
{
  // Constructor
  
}

//_____________________________________________________________________________
AliSurveyPoint::AliSurveyPoint(TString name, Float_t x, Float_t y,
                               Float_t z, Float_t precX, Float_t precY,
                               Float_t precZ, Char_t type, Bool_t Target):
  fPointName(name),
  fX(x),
  fY(y),
  fZ(z),
  fPrecisionX(precX),
  fPrecisionY(precY),
  fPrecisionZ(precZ),
  fType(type),
  fTargetUsed(Target)
{
  // Constructor
  
}

//_____________________________________________________________________________
AliSurveyPoint::~AliSurveyPoint() {
  //destructor
  
}

//_____________________________________________________________________________
void AliSurveyPoint::PrintPoint() {
  // Prints X, Y and Z coordinates of the point
  printf("Point Coordinates \"%s\": %f %f %f\n", (const char*) fPointName, fX, fY, fZ);
  return;
}

