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
Float_t AliMUONResponseV0::IntPH(Float_t eloss)
{
  // Calculate charge from given ionization energy loss
  Int_t nel;
  nel= Int_t(eloss*1.e9/32.);
  Float_t charge=0;
  if (nel == 0) nel=1;
  for (Int_t i=1;i<=nel;i++) {
    charge -= fChargeSlope*TMath::Log(gRandom->Rndm());    
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
    Float_t meanNoise = gRandom->Gaus(1, 0.2);
    Float_t noise     = gRandom->Gaus(0, meanNoise);
    digit+=(Int_t)noise; 
    if ( digit <= ZeroSuppression()) digit = 0;
    if ( digit >  MaxAdc())          digit=MaxAdc();
    return digit;
}









