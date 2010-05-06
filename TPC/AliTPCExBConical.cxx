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
// AliTPCExBConical class                                                   //
// The class calculates the space point distortions due to the conical shape 
// of ALICE TPC:
//
// Becasue of mechanical deformation ALICE TPC, chambers are misaligned in z direction
// TPC has roughly conical shape
//
// For the moment ONLY efective correction used - NOT EDGE EFFECT calcualted //  
//                                                           //
//                     //
// The class allows "effective Omega Tau" corrections.                    // 
// 
//                                                                        //
// date: 02/05/2010                                                       //
// Authors: Marian Ivanov, Jim Thomas, Magnus Mager, Stefan Rossegger     //
//                                                                        //
// Example usage:                                                         //
//  AliTPCExBConical conical;                                                //
//  conical.SetOmegaTauT1T2(0.32,1.,1.); // values ideally from OCDB         //
//  conical.SetXConical(0.001);   // set conical in X direction (in rad)     //
//  // plot dRPhi distortions ...                                            //
//  conical.CreateHistoDRPhiinZR(1.,100,100)->Draw("surf2");                //
////////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "AliTPCROC.h"
#include "AliTPCExBConical.h"
ClassImp(AliTPCExBConical)

AliTPCExBConical::AliTPCExBConical()
  : AliTPCCorrection("exb_conical","ExB conical"),
    fC1(0.),fC2(0.),fConicalFactor(0)
{
  //
  // default constructor
  //
  for (Int_t i=0; i<3; i++){
    fConicalA[i]= 0;  
    fConicalC[i]= 0;
  }
}

AliTPCExBConical::~AliTPCExBConical() {
  //
  // default destructor
  //
}

void AliTPCExBConical::Init() {
  //
  // Initialization funtion (not used at the moment)
  //
  
  // Set default parameters
  // FIXME: Ask the database for these entries
  
  Double_t vdrift = 2.6; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t bzField = -0.5; // [Tesla] // From dataBase: to be updated: per run

  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 

  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  Double_t t1 = 0.9;   // ideally from database
  Double_t t2 = 1.5;   // ideally from database

  SetOmegaTauT1T2(wt,t1,t2);

}

void AliTPCExBConical::Update(const TTimeStamp &/*timeStamp*/) {
  //
  // Update function 
  //

  Double_t vdrift = 2.6; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t bzField = -0.5; // [Tesla] // From dataBase: to be updated: per run

  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 

  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  Double_t t1 = 0.9;   // ideally from database
  Double_t t2 = 1.5;   // ideally from database

 SetOmegaTauT1T2(wt,t1,t2); 
}



void AliTPCExBConical::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the correction due conical shape
  //   
  AliTPCROC * calROC = AliTPCROC::Instance();
  const Double_t kRTPC0  =calROC->GetPadRowRadii(0,0);
  const Double_t kRTPC1  =calROC->GetPadRowRadii(36,calROC->GetNRows(36)-1);
  Float_t rmiddle=(kRTPC0+kRTPC1)/2.;
  //
  Double_t phi =TMath::ATan2(x[1],x[0]);
  Double_t r =TMath::Sqrt(x[1]*x[1]+x[0]*x[0]);
  Double_t dTheta=0;
  if (roc%36<18)  dTheta = fConicalA[0]+TMath::Cos(phi)*fConicalA[1]+TMath::Sin(phi)*fConicalA[2];
  if (roc%36>=18) {
    dTheta = fConicalC[0]+TMath::Cos(phi)*fConicalC[1]+TMath::Sin(phi)*fConicalC[2];
  }
  Double_t corr=dTheta*fConicalFactor;
  if (roc%36>=18) corr*=-1.;
  Double_t drphi=fC1*corr;
  Double_t dr   =fC2*corr;
  dx[0]= TMath::Cos(phi)*dr-TMath::Sin(phi)*drphi;
  dx[1]= TMath::Sin(phi)*dr+TMath::Cos(phi)*drphi;
  dx[2]= -0.001*dTheta*(r-rmiddle); // dtheta in mrad

}

void AliTPCExBConical::Print(Option_t* option) const {
  //
  // Print function to check the settings (e.g. the conical in the X direction)
  // option=="a" prints the C0 and C1 coefficents for calibration purposes
  //

  TString opt = option; opt.ToLower();
  printf("%s\n",GetTitle());
  
  printf("Conical settings: Empirical: %1.5f mm/mrad\n",fConicalFactor);
  printf("Conical settings: A-Conical: %1.5f mrad:%1.5f mrad:%1.5f mrad \n",fConicalA[0],fConicalA[1],fConicalA[2]);
  printf("Conical settings: C-Conical: %1.5f mrad:%1.5f mrad:%1.5f mrad \n",fConicalC[0],fConicalC[1],fConicalC[2]);

}


void AliTPCExBConical::SetConicalA(Float_t conicalA[3]){
  //
  // set paramters of conical shape - A side - obtained from alignment
  //
  fConicalA[0]= conicalA[0];  // constant
  fConicalA[1]= conicalA[1];  // cos 
  fConicalA[2]= conicalA[2];  // sin
}
void AliTPCExBConical::SetConicalC(Float_t conicalC[3]){
  //
  // set paramters of conical shape -C side obtained form the alignemnt
  // 
  fConicalC[0]= conicalC[0];  // constant
  fConicalC[1]= conicalC[1];  // cos 
  fConicalC[2]= conicalC[2];  // sin
}
