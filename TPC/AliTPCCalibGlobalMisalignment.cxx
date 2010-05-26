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
// AliTPCCalibGlobalMisalignment class                                        //
// The class calculates the space point distortions due to simple         // 
// misalignments like shifts in caresian coordinates or a rotation        //
// of the TPC read out planes (A and C side)                              //
//                                                                        //
// date: 06/05/2010                                                       //
// Authors: Stefan Rossegger, Jim Thomas, Magnus Mager                    //
////////////////////////////////////////////////////////////////////////////

#include "AliTPCCalibGlobalMisalignment.h"
#include "TMath.h"

AliTPCCalibGlobalMisalignment::AliTPCCalibGlobalMisalignment()
  : AliTPCCorrection("mialign","Misalignment"),
    fXShift(0.),fYShift(0.),fZShift(0.),
    fRotPhiA(0.),fRotPhiC(0.),
    fdRPhiOffsetA(0.), fdRPhiOffsetC(0.)
{
  //
  // default constructor
  //
}

AliTPCCalibGlobalMisalignment::~AliTPCCalibGlobalMisalignment() {
  //
  // default destructor
  //
}



//void AliTPCCalibGlobalMisalignment::Init() {
//  //
// // Initialization funtion
//  //

//  // nothing to be initialized, results of this calibration class will go to the global aligment structure

//}

//void AliTPCCalibGlobalMisalignment::Update(const TTimeStamp &/*timeStamp*/) {
//  //
//  // Update function 
//  //
//
//  // nothing to be updated, results of this calibration class will go to the global aligment structure
//
//}



void AliTPCCalibGlobalMisalignment::GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // Calculates the simple correction due to a shift (in x,y,z) or an rotation of the TPC (around z)
  // 
 
  Double_t r, phi;
  r   = TMath::Sqrt( x[0]*x[0] + x[1]*x[1] );
  phi = TMath::ATan2(x[1],x[0]);


  // rotation of the read-out planes
  if  (roc%36<18) // A side
    phi += fRotPhiA;
  else         // C side
    phi += fRotPhiC;

  // Simply adding a constant dRPHi residual. PURELY FOR CALIBRATION PURPOSES
  if  (roc%36<18) // A side
    phi += fdRPhiOffsetA/r;
  else         // C side
    phi += fdRPhiOffsetC/r;

  dx[0] = r * TMath::Cos(phi) - x[0];
  dx[1] = r * TMath::Sin(phi) - x[1]; 
  dx[2] = 0.; 

  // Simple shifts
  dx[0] -= fXShift;
  dx[1] -= fYShift;
  dx[2] -= fZShift;

}

void AliTPCCalibGlobalMisalignment::Print(Option_t* /*option*/ ) const {
  //
  // Print function to check the settings 
  //
  printf("%s",GetTitle());  
  printf(" - Trivial Misalignments for calibration purposes: \n");
  printf(" - X-Shift: %1.3f cm, Y-Shift: %1.3f cm, Z-Shift: %1.3f cm \n",fXShift,fYShift,fZShift);
  printf(" - Phi-Rotations: A side: %1.5f rad, C side: %1.5f rad\n",fRotPhiA,fRotPhiC);
  printf(" - dRPhi offsets: A side: %1.5f cm, C side: %1.5f cm\n",fdRPhiOffsetA,fdRPhiOffsetC);
 
 
}
