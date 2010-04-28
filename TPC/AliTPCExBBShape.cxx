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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliExBBShape class                                                     //
// The class calculates the space point distortions due to the B field    //
// shape imperfections using a second order technique based on integrals  //
// over Bz (e.g. int By/Bz) obtained via the AliMagF class                //
// The class allows "effective Omega Tau" corrections.                    //
//                                                                        //
// date: 27/04/2010                                                       //
// Authors: Magnus Mager, Jim Thomas, Stefan Rossegger                    //
//                                                                        //
// Example usage:                                                         //
//  AliMagF mag("mag","mag");                                             //
//  AliTPCExBBShape exb;                                                  //
//  exb.SetBField(&mag);             // use Bfield from AliMagF           //
//  exb.SetOmegaTauT1T2(0.32,1.,1.); // values ideally from OCDB          //
//  // plot dRPhi distortions ...                                         //
//  exb.CreateHistoDRPhiinZR(0,100,100)->Draw("surf2");                   //
////////////////////////////////////////////////////////////////////////////

#include <AliMagF.h>

#include "AliTPCExBBShape.h"

AliTPCExBBShape::AliTPCExBBShape()
  : AliTPCCorrection("exb_bshape","ExB B-shape"),
    fC1(0.),fC2(0.),
    fBField(0)
{
  //
  // default constructor
  //
}

AliTPCExBBShape::~AliTPCExBBShape() {
  //
  // virtual destructor
  //
}

void AliTPCExBBShape::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the space point corrections of the B field inperfections (B field shape) 
  //

  if (!fBField) {
    for (Int_t j=0;j<3;++j) dx[j]=0.;
    return;
  }

  const Double_t xStart[3]={ x[0], x[1], x[2] };
  const Double_t xEnd[3]={ x[0],  x[1],  roc%36<18?fgkTPC_Z0:-fgkTPC_Z0 };

  Double_t intBStart[3];
  Double_t intBEnd[3];

  fBField->GetTPCRatInt(xStart,intBStart);
  fBField->GetTPCRatInt(xEnd,  intBEnd  );

  const Float_t intBxOverBz=(intBEnd[0]-intBStart[0]);
  const Float_t intByOverBz=(intBEnd[1]-intBStart[1]);
  
  dx[0]=fC2*intBxOverBz-fC1*intByOverBz;
  dx[1]=fC1*intBxOverBz+fC2*intByOverBz;
  dx[2]=0.;


}

void AliTPCExBBShape::GetBxAndByOverBz(const Float_t x[],const Short_t roc,Float_t BxByOverBz[]) {
  //
  // This function is purely for calibration purposes
  // Returns the via AliMagF obtaind B field integrals  
  // 

  if (!fBField) {
    for (Int_t j=0;j<3;++j) BxByOverBz[j]=0.;
    return;
  }

  const Double_t xStart[3]={ x[0], x[1], x[2] };
  const Double_t xEnd[3]={ x[0],  x[1],  roc%36<18?fgkTPC_Z0:-fgkTPC_Z0 };

  Double_t intBStart[3];
  Double_t intBEnd[3];

  fBField->GetTPCRatInt(xStart,intBStart);
  fBField->GetTPCRatInt(xEnd,  intBEnd  );

  const Float_t intBxOverBz=(intBEnd[0]-intBStart[0]);
  const Float_t intByOverBz=(intBEnd[1]-intBStart[1]);
  
  BxByOverBz[0]=intBxOverBz;
  BxByOverBz[1]=intByOverBz;

}

void AliTPCExBBShape::Print(Option_t* option) const {
  //
  // Print function to check the settings (e.g. voltage offsets)
  // option=="a" prints details of the B field settings and the 
  // C0 and C1 coefficents (for calibration purposes)
  //
  TString opt = option; opt.ToLower();
  printf("%s\n - B field settings:\n",GetTitle());
  fBField->Print(option);
  //  printf(" - B field: X-Twist: %1.5lf rad, Y-Twist: %1.5lf rad \n",fBField->Print(option));
  if (opt.Contains("a")) { // Print all details
    printf(" - C1: %1.4f, C2: %1.4f \n",fC1,fC2);
  }    
 
}
