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
// AliTPCExBBShape class                                                  //
////////////////////////////////////////////////////////////////////////////

#include <AliMagF.h>
#include "TGeoGlobalMagField.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliLog.h"

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

void AliTPCExBBShape::Init() {
  //
  // Initialization funtion
  //
  
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magneticd field - not initialized");
  Double_t bzField = magF->SolenoidField()/10.; //field in T
  SetBField(magF);
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 
  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  SetOmegaTauT1T2(wt,fT1,fT2);


}

void AliTPCExBBShape::Update(const TTimeStamp &/*timeStamp*/) {
  //
  // Update function 
  //
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magneticd field - not initialized");
  Double_t bzField = magF->SolenoidField()/10.; //field in T
  SetBField(magF);
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 
  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  SetOmegaTauT1T2(wt,fT1,fT2);


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
  const Double_t xEnd[3]={ x[0],  x[1],  roc%36<18?fgkTPCZ0:-fgkTPCZ0 };

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
  const Double_t xEnd[3]={ x[0],  x[1],  roc%36<18?fgkTPCZ0:-fgkTPCZ0 };

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
  printf("%s\t%s\n - B field settings:\n",GetTitle(),GetName());
  fBField->Print(option);
  //  printf(" - B field: X-Twist: %1.5lf rad, Y-Twist: %1.5lf rad \n",fBField->Print(option));
  if (opt.Contains("a")) { // Print all details  
    printf(" - T1: %1.4f, T2: %1.4f \n",fT1,fT2);
    printf(" - C1: %1.4f, C2: %1.4f \n",fC1,fC2);
  }    
}

Double_t AliTPCExBBShape::GetBFieldXYZ(Double_t gx, Double_t gy, Double_t gz, Int_t axisType){
  //
  // return B field at given x,y,z
  // 
  AliMagF* field = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!field) return 0;
  Double_t xyz[3]={gx,gy,gz};
  Double_t bxyz[3]={0};
  field->Field(xyz,bxyz);
  //
  Double_t r=TMath::Sqrt(gx*gx+gy*gy);
  //  Double_t b=TMath::Sqrt(bxyz[0]*bxyz[0]+bxyz[1]*bxyz[1]);
  if (axisType==0) {
    return (xyz[0]*bxyz[1]-xyz[1]*bxyz[0])/(bxyz[2]*r);
  }
  if (axisType==1){
    return (xyz[0]*bxyz[0]+xyz[1]*bxyz[1])/(bxyz[2]*r);
  }
  if (axisType==2) return bxyz[2];
  if (axisType==3) return bxyz[0];
  if (axisType==4) return bxyz[1];
  if (axisType==5) return bxyz[2];
  return bxyz[2];
}
