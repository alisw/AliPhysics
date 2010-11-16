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
// AliTPCExBTwist class                                                   //
////////////////////////////////////////////////////////////////////////////


#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliLog.h"

#include "AliTPCExBTwist.h"

AliTPCExBTwist::AliTPCExBTwist()
  : AliTPCCorrection("exb_twist","ExB twist"),
    fC1(0.),fC2(0.),
    fXTwist(0.),fYTwist(0.)
{
  //
  // default constructor
  //
}

AliTPCExBTwist::~AliTPCExBTwist() {
  //
  // default destructor
  //
}



void AliTPCExBTwist::Init() {
  //
  // Initialization funtion
  //
  
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magneticd field - not initialized");
  Double_t bzField = magF->SolenoidField()/10.; //field in T
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 
  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  SetOmegaTauT1T2(wt,fT1,fT2);


}

void AliTPCExBTwist::Update(const TTimeStamp &/*timeStamp*/) {
  //
  // Update function 
  //
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magneticd field - not initialized");
  Double_t bzField = magF->SolenoidField()/10.; //field in T
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 
  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  SetOmegaTauT1T2(wt,fT1,fT2);


}



void AliTPCExBTwist::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the correction of a mismatch between the E and B field axis
  // 
  
  const Float_t zstart=x[2];
  const Float_t zend  =(roc%36<18?fgkTPCZ0:-fgkTPCZ0);
  const Float_t zdrift=zstart-zend;
  
  dx[0]=(fC2*fXTwist-fC1*fYTwist)*zdrift;
  dx[1]=(fC1*fXTwist+fC2*fYTwist)*zdrift;
  dx[2]=0.;
}

void AliTPCExBTwist::Print(const Option_t* option) const {
  //
  // Print function to check the settings (e.g. the twist in the X direction)
  // option=="a" prints the C0 and C1 coefficents for calibration purposes
  //

  TString opt = option; opt.ToLower();
  printf("%s\n",GetTitle());
  
  printf(" - Twist settings: X-Twist: %1.5f rad, Y-Twist: %1.5f rad \n",fXTwist,fYTwist);
  if (opt.Contains("a")) { // Print all details
    printf(" - T1: %1.4f, T2: %1.4f \n",fT1,fT2);
    printf(" - C1: %1.4f, C2: %1.4f \n",fC1,fC2);
  }    
 
 
}
