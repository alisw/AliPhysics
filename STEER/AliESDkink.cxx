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

/* $Id$ */

//-------------------------------------------------------------------------
//    Origin: Marian Ivanov marian.ivanov@cern.ch
//-------------------------------------------------------------------------

#include <Riostream.h>
#include <TMath.h>
#include <TPDGCode.h>
#include "AliESDkink.h"
#include "AliHelix.h"


ClassImp(AliESDkink)

AliESDkink::AliESDkink(){
  //
  //Dafault constructor
  //
  fID = -1;
  fDist1  = -1;
  fDist2  = -1;
  fPdr[0] = 0;
  fPdr[1] = 0;
  fPdr[2] = 0;
  fXr[0] = 0;
  fXr[1] = 0;
  fXr[2] = 0;
  fPm[0] = 0;
  fPm[1] = 0;
  fPm[2] = 0;
  fAngle[0] = 0;
  fAngle[1] = 0;
  fAngle[2] = 0;
  fRr     = -1;
  fLab[0] = -1;
  fLab[1] = -1;
  fIndex[0] = -1;
  fIndex[1] = -1;
  fStatus = 0;
  fTPCdensity[0][0]=-1;
  fTPCdensity[0][1]=-1;
  fTPCdensity[1][0]=-1;
  fTPCdensity[1][1]=-1;
  fTPCdensity2[0][0]=-1;
  fTPCdensity2[0][1]=-1;
  fTPCdensity2[1][0]=-1;
  fTPCdensity2[1][1]=-1;
  fShapeFactor = 0;
  fRow0   =-1;
  fMultiple[0] = 0;
  fMultiple[1] = 0;
  fZm[0] = 0;
  fZm[1] = 0;
  fFi[0] = 0;
  fFi[1] = 0;
}

void AliESDkink::SetMother(const AliExternalTrackParam & pmother)  {
  //
  // set mother
  //
  fParamMother   = pmother;
}

void AliESDkink::SetDaughter(const AliExternalTrackParam & pdaughter){
  //
  //set daughter
  //
  fParamDaughter = pdaughter;

}
  
void  AliESDkink::Update()
{
  //
  // updates Kink Info
  //
  Float_t distance1,distance2;
  //
  AliHelix dhelix1(fParamDaughter);
  AliHelix mhelix(fParamMother);    
  //
  //find intersection linear
  //
  Double_t phase[2][2],radius[2];
  Int_t  points = dhelix1.GetRPHIintersections(mhelix, phase, radius,200);
  Double_t delta1=10000,delta2=10000;  
  
  if (points>0){
    dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    dhelix1.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
  }
  if (points==2){    
    dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    dhelix1.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
  }
  distance1 = TMath::Min(delta1,delta2);
  //
  //find intersection parabolic
  //
  points = dhelix1.GetRPHIintersections(mhelix, phase, radius);
  delta1=10000,delta2=10000;  
  Double_t d1=1000.,d2=10000.;
  if (points>0){
    dhelix1.ParabolicDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    dhelix1.ParabolicDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    Double_t xd[3],xm[3];
    dhelix1.Evaluate(phase[0][0],xd);
    mhelix.Evaluate(phase[0][1],xm);
    d1 = (xd[0]-xm[0])*(xd[0]-xm[0])+(xd[1]-xm[1])*(xd[1]-xm[1])+(xd[2]-xm[2])*(xd[2]-xm[2]);
  }
  if (points==2){    
    dhelix1.ParabolicDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    dhelix1.ParabolicDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    Double_t xd[3],xm[3];
    dhelix1.Evaluate(phase[1][0],xd);
    mhelix.Evaluate(phase[1][1],xm);
    d2 = (xd[0]-xm[0])*(xd[0]-xm[0])+(xd[1]-xm[1])*(xd[1]-xm[1])+(xd[2]-xm[2])*(xd[2]-xm[2]);
  }
  //
  distance2 = TMath::Min(delta1,delta2);
  if (delta1<delta2){
    //get V0 info
    //    dhelix1.Evaluate(phase[0][0],fXr);
    Double_t xd[3],xm[3];
    dhelix1.Evaluate(phase[0][0],xd);
    mhelix.Evaluate(phase[0][1], xm);
    fXr[0] = 0.5*(xd[0]+xm[0]);
    fXr[1] = 0.5*(xd[1]+xm[1]);
    fXr[2] = 0.5*(xd[2]+xm[2]);
    //
    dhelix1.GetMomentum(phase[0][0],fPdr);
    mhelix.GetMomentum(phase[0][1],fPm);
    dhelix1.GetAngle(phase[0][0],mhelix,phase[0][1],fAngle);
    //fRr = TMath::Sqrt(radius[0]);
    fRr = TMath::Sqrt(fXr[0]*fXr[0]+fXr[1]*fXr[1]);
  }
  else{
    //dhelix1.Evaluate(phase[1][0],fXr);
    Double_t xd[3],xm[3];
    dhelix1.Evaluate(phase[1][0],xd);
    mhelix.Evaluate(phase[1][1], xm);
    fXr[0] = 0.5*(xd[0]+xm[0]);
    fXr[1] = 0.5*(xd[1]+xm[1]);
    fXr[2] = 0.5*(xd[2]+xm[2]);
    //
    dhelix1.GetMomentum(phase[1][0], fPdr);
    mhelix.GetMomentum(phase[1][1], fPm);
    dhelix1.GetAngle(phase[1][0],mhelix,phase[1][1],fAngle);
    //    fRr = TMath::Sqrt(radius[1]); 
    fRr = TMath::Sqrt(fXr[0]*fXr[0]+fXr[1]*fXr[1]);
  }
  fDist1 = TMath::Sqrt(TMath::Min(d1,d2));
  fDist2 = TMath::Sqrt(distance2);      
  //            
  //

}

 
Float_t AliESDkink::GetTPCDensityFactor() const
{
  //
  //
  return fTPCdensity[0][0]+fTPCdensity[1][1]-TMath::Max(fTPCdensity[0][1],Float_t(0.0))-TMath::Max(fTPCdensity[1][0],Float_t(0.0)); 
}

Float_t AliESDkink::GetQt() const
{
  Float_t dmomentum = TMath::Sqrt(fPdr[0]*fPdr[0]+fPdr[1]*fPdr[1]+fPdr[2]*fPdr[2]);
  return TMath::Sin(fAngle[2])*dmomentum;
}
