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

*/

#include "AliMUONResponseTriggerV1.h"
#include "AliSegmentation.h"
#include <TMath.h>
#include <TRandom.h>
#include <iostream.h> 

ClassImp(AliMUONResponseTriggerV1)

//------------------------------------------------------------------   
AliMUONResponseTriggerV1::AliMUONResponseTriggerV1(){
// default constructor 
  Float_t hv=9.2;
  SetParameters(hv);
}

//------------------------------------------------------------------   
AliMUONResponseTriggerV1::AliMUONResponseTriggerV1(Float_t hv){
// Constructor 
  SetParameters(hv);
}

//------------------------------------------------------------------   
void AliMUONResponseTriggerV1::SetParameters(Float_t hv){
// initialize parameters accoring to HV
// (see V.Barret B.Espagnon and P.Rosnet Alice/note xxx)
  fA = 6.089 * hv - 52.70;
  fB = 2.966;
  fC = 4.3e-4 * hv - 3.5e-3;
}

//------------------------------------------------------------------   
Int_t AliMUONResponseTriggerV1::SetGenerCluster(){
// Set the GenerCluster parameter and return 1
  fGenerCluster = gRandom->Rndm();
  return 1;
} 

//------------------------------------------------------------------   
Float_t AliMUONResponseTriggerV1::IntXY(AliSegmentation * segmentation){
// Returns 1 or 0 if the current strip is fired or not 
// get the "parameters" needed to evaluate the strip response
// x1 : hit x(y) position
// x2 : x(y) coordinate of the main strip
// x3 : current strip real x(y) coordinate  
// x4 : dist. between x(y) hit pos. and the closest border of the current strip

  Float_t x1,x2,x3,x4;  
  segmentation->IntegrationLimits(x1,x2,x3,x4);    
  Float_t theta = 0.; // incident angle to be implemented

  return (fGenerCluster < FireStripProb(x4,theta)) ? 1:0; 
}

//------------------------------------------------------------------   
Float_t AliMUONResponseTriggerV1::FireStripProb(Float_t x4, Float_t theta){
// parametrisation of the probability that a strip neighbour of the main 
// strip is fired (V.Barret B.Espagnon and P.Rosnet Alice/note xxx)
// WARNING : need to convert x4 from cm to mm

  return 
    (TMath::Cos(theta)*fA/(fA+TMath::Cos(theta)*TMath::Power(x4*10.,fB))+fC)/
    (TMath::Cos(theta)+fC);
}

