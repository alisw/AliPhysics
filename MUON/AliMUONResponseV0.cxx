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

#include <TMath.h>
#include <TRandom.h>

#include "AliMUONResponseV0.h"
#include "AliSegmentation.h"
#include "AliMUONGeometrySegmentation.h"

ClassImp(AliMUONResponseV0)
	
//__________________________________________________________________________
AliMUONResponseV0::AliMUONResponseV0()
  : AliMUONResponse()
{
// Default constructor

  fMathieson = new AliMUONMathieson();
  fChargeCorrel = 0;
}
   //__________________________________________________________________________
AliMUONResponseV0::~AliMUONResponseV0()
{
  delete fMathieson;
}
  //__________________________________________________________________________
void AliMUONResponseV0::SetSqrtKx3AndDeriveKx2Kx4(Float_t SqrtKx3)
{
  // Set to "SqrtKx3" the Mathieson parameter K3 ("fSqrtKx3")
  // in the X direction, perpendicular to the wires,
  // and derive the Mathieson parameters K2 ("fKx2") and K4 ("fKx4")
  // in the same direction
  fMathieson->SetSqrtKx3AndDeriveKx2Kx4(SqrtKx3);
}
	
  //__________________________________________________________________________
void AliMUONResponseV0::SetSqrtKy3AndDeriveKy2Ky4(Float_t SqrtKy3)
{
  // Set to "SqrtKy3" the Mathieson parameter K3 ("fSqrtKy3")
  // in the Y direction, along the wires,
  // and derive the Mathieson parameters K2 ("fKy2") and K4 ("fKy4")
  // in the same direction
  fMathieson->SetSqrtKy3AndDeriveKy2Ky4(SqrtKy3);
}
  //__________________________________________________________________________
Float_t AliMUONResponseV0::IntPH(Float_t eloss)
{
  // Calculate charge from given ionization energy loss
  Int_t nel;
  nel= Int_t(eloss*1.e9/27.4);
  Float_t charge=0;
  if (nel == 0) nel=1;
  for (Int_t i=1;i<=nel;i++) {
      Float_t arg=0.;
      while(!arg) arg = gRandom->Rndm();
      charge -= fChargeSlope*TMath::Log(arg);    
  }
  return charge;
}
  //-------------------------------------------
Float_t AliMUONResponseV0::IntXY(AliSegmentation * segmentation)
{
  // Calculate charge on current pad according to Mathieson distribution

  return fMathieson->IntXY(segmentation);

}
  //-------------------------------------------
Float_t AliMUONResponseV0::IntXY(Int_t idDE, AliMUONGeometrySegmentation* segmentation)
{
 // Calculate charge on current pad according to Mathieson distribution

  return fMathieson->IntXY(idDE, segmentation);
}
  //-------------------------------------------
Int_t  AliMUONResponseV0::DigitResponse(Int_t digit, AliMUONTransientDigit* /*where*/)
{
    // add white noise and do zero-suppression and signal truncation
//     Float_t meanNoise = gRandom->Gaus(1, 0.2);
    // correct noise for slat chambers;
    // one more field to add to AliMUONResponseV0 to allow different noises ????
    Float_t meanNoise = gRandom->Gaus(1., 0.2);
    Float_t noise     = gRandom->Gaus(0., meanNoise);
    digit += TMath::Nint(noise); 
    if ( digit <= ZeroSuppression()) digit = 0;
    // if ( digit >  MaxAdc())          digit=MaxAdc();
    if ( digit >  Saturation())          digit=Saturation();

    return digit;
}









