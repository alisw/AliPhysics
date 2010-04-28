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
// AliTPCExBTwist class                                                   //
// The class calculates the space point distortions due to a mismatch     //
// of the E and B field axis (original code from STAR)                    //
// The class allows "effective Omega Tau" corrections.                    // 
//                                                                        //
// date: 27/04/2010                                                       //
// Authors: Jim Thomas, Magnus Mager, Stefan Rossegger                    //
//                                                                        //
// Example usage:                                                         //
//  AliTPCExBTwist twist;                                                 //
//  twist.SetOmegaTauT1T2(0.32,1.,1.); // values ideally from OCDB        //
//  twist.SetXTwist(0.001);   // set twist in X direction (in rad)        //
//  // plot dRPhi distortions ...                                         //
//  twist.CreateHistoDRPhiinZR(1.,100,100)->Draw("surf2");                //
////////////////////////////////////////////////////////////////////////////

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

void AliTPCExBTwist::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the correction of a mismatch between the E and B field axis
  // 
  
  const Float_t zstart=x[2];
  const Float_t zend  =(roc%36<18?fgkTPC_Z0:-fgkTPC_Z0);
  const Float_t zdrift=zstart-zend;
  
  dx[0]=(fC2*fXTwist-fC1*fYTwist)*zdrift;
  dx[1]=(fC1*fXTwist+fC2*fYTwist)*zdrift;
  dx[2]=0.;
}

void AliTPCExBTwist::Print(Option_t* option) const {
  //
  // Print function to check the settings (e.g. the twist in the X direction)
  // option=="a" prints the C0 and C1 coefficents for calibration purposes
  //

  TString opt = option; opt.ToLower();
  printf("%s\n",GetTitle());
  
  printf(" - Twist settings: X-Twist: %1.5f rad, Y-Twist: %1.5f rad \n",fXTwist,fYTwist);
  if (opt.Contains("a")) { // Print all details
    printf(" - C1: %1.4f, C2: %1.4f \n",fC1,fC2);
  }    
 
 
}
