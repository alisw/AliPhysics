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

/*
$Log$
Revision 1.5  2000/11/21 13:47:55  gosset
All Mathieson parameters (Sqrt(K3), K2 and K4) set in one function,
SetSqrtKx3AndDeriveKx2Kx4 or SetSqrtKx3AndDeriveKx2Kx4,
for each cathode plane

Revision 1.4  2000/10/25 10:41:52  morsch
IntPH(..): Protec Log against random numbers equal to 0.

Revision 1.3  2000/07/03 11:54:57  morsch
AliMUONSegmentation and AliMUONHitMap have been replaced by AliSegmentation and AliHitMap in STEER
The methods GetPadIxy and GetPadXxy of AliMUONSegmentation have changed name to GetPadI and GetPadC.

Revision 1.2  2000/06/15 07:58:48  morsch
Code from MUON-dev joined

Revision 1.1.2.1  2000/06/09 21:33:35  morsch
AliMUONResponse code  from  AliMUONSegResV0.cxx

*/

#include "AliMUONResponseV0.h"
#include "AliSegmentation.h"
#include <TMath.h>
#include <TRandom.h>


ClassImp(AliMUONResponseV0)
	
  //__________________________________________________________________________
void AliMUONResponseV0::SetSqrtKx3AndDeriveKx2Kx4(Float_t SqrtKx3)
{
  // Set to "SqrtKx3" the Mathieson parameter K3 ("fSqrtKx3")
  // in the X direction, perpendicular to the wires,
  // and derive the Mathieson parameters K2 ("fKx2") and K4 ("fKx4")
  // in the same direction
  fSqrtKx3 = SqrtKx3;
  fKx2 = TMath::Pi() / 2. * (1. - 0.5 * fSqrtKx3);
  Float_t cx1 = fKx2 * fSqrtKx3 / 4. / TMath::ATan(Double_t(fSqrtKx3));
  fKx4 = cx1 / fKx2 / fSqrtKx3;
}
	
  //__________________________________________________________________________
void AliMUONResponseV0::SetSqrtKy3AndDeriveKy2Ky4(Float_t SqrtKy3)
{
  // Set to "SqrtKy3" the Mathieson parameter K3 ("fSqrtKy3")
  // in the Y direction, along the wires,
  // and derive the Mathieson parameters K2 ("fKy2") and K4 ("fKy4")
  // in the same direction
  fSqrtKy3 = SqrtKy3;
  fKy2 = TMath::Pi() / 2. * (1. - 0.5 * fSqrtKy3);
  Float_t cy1 = fKy2 * fSqrtKy3 / 4. / TMath::ATan(Double_t(fSqrtKy3));
  fKy4 = cy1 / fKy2 / fSqrtKy3;
}

Float_t AliMUONResponseV0::IntPH(Float_t eloss)
{
  // Calculate charge from given ionization energy loss
  Int_t nel;
  nel= Int_t(eloss*1.e9/32.);
  Float_t charge=0;
  if (nel == 0) nel=1;
  for (Int_t i=1;i<=nel;i++) {
      Float_t arg=0.;
      while(!arg) arg = gRandom->Rndm();
      charge -= fChargeSlope*TMath::Log(arg);    
  }
  return charge;
}
// -------------------------------------------

Float_t AliMUONResponseV0::IntXY(AliSegmentation * segmentation)
{
// Calculate charge on current pad according to Mathieson distribution
// 
    const Float_t kInversePitch = 1/fPitch;
//
//  Integration limits defined by segmentation model
//  
    Float_t xi1, xi2, yi1, yi2;
    segmentation->IntegrationLimits(xi1,xi2,yi1,yi2);
    xi1=xi1*kInversePitch;
    xi2=xi2*kInversePitch;
    yi1=yi1*kInversePitch;
    yi2=yi2*kInversePitch;
//
// The Mathieson function 
    Double_t ux1=fSqrtKx3*TMath::TanH(fKx2*xi1);
    Double_t ux2=fSqrtKx3*TMath::TanH(fKx2*xi2);

    Double_t uy1=fSqrtKy3*TMath::TanH(fKy2*yi1);
    Double_t uy2=fSqrtKy3*TMath::TanH(fKy2*yi2);

    
    return Float_t(4.*fKx4*(TMath::ATan(ux2)-TMath::ATan(ux1))*
		      fKy4*(TMath::ATan(uy2)-TMath::ATan(uy1)));
}

Int_t  AliMUONResponseV0::DigitResponse(Int_t digit)
{
    // add white noise and do zero-suppression and signal truncation
//     Float_t meanNoise = gRandom->Gaus(1, 0.2);
    // correct noise for slat chambers;
    // one more field to add to AliMUONResponseV0 to allow different noises ????
    Float_t meanNoise = gRandom->Gaus(1.5, 0.2);
    Float_t noise     = gRandom->Gaus(0, meanNoise);
    digit+=(Int_t)noise; 
    if ( digit <= ZeroSuppression()) digit = 0;
    if ( digit >  MaxAdc())          digit=MaxAdc();
    return digit;
}









