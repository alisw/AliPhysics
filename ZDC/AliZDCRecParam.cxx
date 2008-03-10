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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for ZDC calibration                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliZDCRecParam.h"
#include "TH1.h"
#include "AliCDBEntry.h"

ClassImp(AliZDCRecParam)

//________________________________________________________________
AliZDCRecParam::AliZDCRecParam():
TNamed(),
fZEMEndValue(0),
fZEMCutFraction(0),
fDZEMSup(0),
fDZEMInf(0),
fEZN1MaxValue(0),
fEZP1MaxValue(0),
fEZDC1MaxValue(0),
fEZN2MaxValue(0),
fEZP2MaxValue(0),
fEZDC2MaxValue(0)
{
  Reset();
}

//________________________________________________________________
AliZDCRecParam::AliZDCRecParam(const char* name):
TNamed(),
fZEMEndValue(0),
fZEMCutFraction(0),
fDZEMSup(0),
fDZEMInf(0),
fEZN1MaxValue(0),
fEZP1MaxValue(0),
fEZDC1MaxValue(0),
fEZN2MaxValue(0),
fEZP2MaxValue(0),
fEZDC2MaxValue(0)
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
}

//________________________________________________________________
AliZDCRecParam::AliZDCRecParam(const AliZDCRecParam& calibda) :
TNamed(calibda),
fZEMEndValue(calibda.GetZEMEndValue()),  
fZEMCutFraction(calibda.GetZEMCutFraction()),
fDZEMSup(calibda.GetDZEMSup()),
fDZEMInf(calibda.GetDZEMInf()),
fEZN1MaxValue(0),
fEZP1MaxValue(0),
fEZDC1MaxValue(0),
fEZN2MaxValue(0),
fEZP2MaxValue(0),
fEZDC2MaxValue(0)
{
  // Copy constructor
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
}

//________________________________________________________________
AliZDCRecParam &AliZDCRecParam::operator =(const AliZDCRecParam& calibda)
{
// assignment operator
  SetName(calibda.GetName());
  SetTitle(calibda.GetName());
  Reset();
  fZEMEndValue    = calibda.GetZEMEndValue();
  fZEMCutFraction = calibda.GetZEMCutFraction();

  return *this;
}

//________________________________________________________________
AliZDCRecParam::~AliZDCRecParam()
{
}

//________________________________________________________________
void AliZDCRecParam::Reset()
{
  // Reset
}                                                                                       


//________________________________________________________________
void  AliZDCRecParam::Print(Option_t *) const
{
   // Printing calibration object
   printf("\n\n ####### Parameters from EZDC vs. ZEM correlation #######	\n");
   printf("  ZEMEndPoint = %1.2f, ZEMCutFraction = %1.2f \n"
     "  DZEMInf = %1.2f, DZEMSup = %1.2f\n",
     fZEMEndValue, fZEMCutFraction, fDZEMInf, fDZEMSup);
 
   printf("\n\n ####### Parameters from EZDC vs. Nspec correlation #######	\n");
   printf("  EZN1MaxValue = %1.2f, EZP1MaxValue = %1.2f, EZDC1MaxValue = %1.2f \n"
     "  EZN2MaxValue = %1.2f, EZP2MaxValue = %1.2f, EZDC2MaxValue = %1.2f \n\n",
     fEZN1MaxValue, fEZP1MaxValue, fEZDC1MaxValue,
     fEZN2MaxValue, fEZP2MaxValue, fEZDC2MaxValue);

} 

