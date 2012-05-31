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
// AliTPCCorrectionDrift class  
// linear drift corrections                                               //
//  
////////////////////////////////////////////////////////////////////////////
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliLog.h"

#include "TMath.h"
#include "AliTPCROC.h"
#include "AliTPCCorrectionDrift.h"
ClassImp(AliTPCCorrectionDrift)

AliTPCCorrectionDrift::AliTPCCorrectionDrift()
  : AliTPCCorrection("CorrectionDrift","CorrectionDrift") ,
    fZ0Aside(0),     // z- t0*vdrift shift A side
    fZ0Cside(0),     // z- t0*vdrift shift C side
    fVScale0(0),     // drift velocity scaling - constant
    fVScaleR(0),     // drift velocity scaling - radial
    fVScaleX(0),     // drift velocity scaling - global x
    fVScaleY(0),     // drift velocity scaling - global y
    fIROCZ0(0),
    fOROCDZ(0)
{
  //
  // default constructor
  //
}

AliTPCCorrectionDrift::~AliTPCCorrectionDrift() {
  //
  // default destructor
  //
}



void AliTPCCorrectionDrift::Init() {
  //
  // Initialization funtion
  //
  


}

void AliTPCCorrectionDrift::Update(const TTimeStamp &/*timeStamp*/) {
  //
  // Update function 
  //

}



void AliTPCCorrectionDrift::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the correction due conical shape
  //   
  AliTPCROC * calROC = AliTPCROC::Instance();
  //const Double_t kRTPC0  =calROC->GetPadRowRadii(0,0);
  const Double_t kRTPC1  =calROC->GetPadRowRadii(36,calROC->GetNRows(36)-1);
  const Double_t kRIO    =0.5*(calROC->GetPadRowRadii(0,calROC->GetNRows(0)-1)+calROC->GetPadRowRadii(36,0));
  //  Float_t rmiddle=(kRTPC0+kRTPC1)/2.;
  //
  //Double_t phi      = TMath::ATan2(x[1],x[0]);
  Double_t r        = TMath::Sqrt(x[1]*x[1]+x[0]*x[0]);
  Double_t driftN   = 1.-TMath::Abs(x[2])/calROC->GetZLength(0);  // drift from 0 to 1
  //
  Double_t dz0   =(roc%36<18) ? fZ0Aside:-fZ0Cside;
  Double_t dscale= (fVScale0+(fVScaleR*r+fVScaleX*x[0]+fVScaleY*x[1])/kRTPC1);
  Double_t ddrift=(roc%36<18) ? driftN*dscale*calROC->GetZLength(0):-driftN*dscale*calROC->GetZLength(0);
  if (r<kRIO) dz0+=(roc%36<18) ? fIROCZ0:-fIROCZ0;
  if (r>kRIO) dz0+=(roc%36<18) ? fOROCDZ*(r-kRIO):-fOROCDZ*(r-kRIO);
  // Calculate correction in cartesian coordinates
  dx[0] = 0;
  dx[1] = 0;
  dx[2] = dz0+ddrift; // z distortion not implemented (1st order distortions)

}





void AliTPCCorrectionDrift::Print(const Option_t* option) const {
  //
  // Print function to check the settings (e.g. the twist in the X direction)
  // 
  //

  TString opt = option; opt.ToLower();
  printf("%s\t%s\n",GetName(),GetTitle());
  
  if (opt.Contains("a")) { // Print all details
    printf(" - T0A: %1.4f, T0C: %1.4f (cm)\n",fZ0Aside,fZ0Cside);
    printf(" - Scale0: %1.4f, ScaleR: %1.4f \n",fVScale0,fVScaleR);
    printf(" - ScaleX: %1.4f, ScaleY: %1.4f \n",fVScaleX,fVScaleY);

  }    
 
 
}
