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

//-------------------------------------------------------------------------
//
//    Implementation of the V0 vertex class
//    Numerical part - AliHelix functionality used             
//
//    Origin: Marian Ivanov marian.ivanov@cern.ch
//-------------------------------------------------------------------------
#include <TMath.h>

#include "AliV0.h"
#include "AliHelix.h"


ClassImp(AliV0)

void  AliV0::Update(Float_t vertex[3])
{
  //
  // updates Kink Info
  //
  //  Float_t distance1,distance2;
  Float_t distance2;
  //
  AliHelix phelix(fParamP);
  AliHelix mhelix(fParamN);    
  //
  //find intersection linear
  //
  Double_t phase[2][2],radius[2];
  Int_t  points = phelix.GetRPHIintersections(mhelix, phase, radius,200);
  Double_t delta1=10000,delta2=10000;  
  /*
  if (points<=0) return;
  if (points>0){
    phelix.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    phelix.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    phelix.LinearDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
  }
  if (points==2){    
    phelix.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    phelix.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    phelix.LinearDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
  }
  distance1 = TMath::Min(delta1,delta2);
  */
  //
  //find intersection parabolic
  //
  points = phelix.GetRPHIintersections(mhelix, phase, radius);
  delta1=10000,delta2=10000;  
  Double_t d1=1000.,d2=10000.;
  Double_t err[3],angles[3];
  if (points<=0) return;
  if (points>0){
    phelix.ParabolicDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    phelix.ParabolicDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    if (TMath::Abs(fParamP.GetX()-TMath::Sqrt(radius[0])<3) && TMath::Abs(fParamN.GetX()-TMath::Sqrt(radius[0])<3)){
      // if we are close to vertex use error parama
      //
      err[1] = fParamP.GetCovariance()[0]+fParamN.GetCovariance()[0]+0.05*0.05
	+0.3*(fParamP.GetCovariance()[2]+fParamN.GetCovariance()[2]);
      err[2] = fParamP.GetCovariance()[2]+fParamN.GetCovariance()[2]+0.05*0.05
	+0.3*(fParamP.GetCovariance()[0]+fParamN.GetCovariance()[0]);
      
      phelix.GetAngle(phase[0][0],mhelix,phase[0][1],angles);
      Double_t tfi  = TMath::Abs(TMath::Tan(angles[0]));
      Double_t tlam = TMath::Abs(TMath::Tan(angles[1]));
      err[0] = err[1]/((0.2+tfi)*(0.2+tfi))+err[2]/((0.2+tlam)*(0.2+tlam));
      err[0] = ((err[1]*err[2]/((0.2+tfi)*(0.2+tfi)*(0.2+tlam)*(0.2+tlam))))/err[0];
      phelix.ParabolicDCA2(mhelix,phase[0][0],phase[0][1],radius[0],delta1,err);
    }
    Double_t xd[3],xm[3];
    phelix.Evaluate(phase[0][0],xd);
    mhelix.Evaluate(phase[0][1],xm);
    d1 = (xd[0]-xm[0])*(xd[0]-xm[0])+(xd[1]-xm[1])*(xd[1]-xm[1])+(xd[2]-xm[2])*(xd[2]-xm[2]);
  }
  if (points==2){    
    phelix.ParabolicDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    phelix.ParabolicDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    if (TMath::Abs(fParamP.GetX()-TMath::Sqrt(radius[1])<3) && TMath::Abs(fParamN.GetX()-TMath::Sqrt(radius[1])<3)){
      // if we are close to vertex use error paramatrization
      //
      err[1] = fParamP.GetCovariance()[0]+fParamN.GetCovariance()[0]+0.05*0.05
	+0.3*(fParamP.GetCovariance()[2]+fParamN.GetCovariance()[2]);
      err[2] = fParamP.GetCovariance()[2]+fParamN.GetCovariance()[2]+0.05*0.05
	+0.3*(fParamP.GetCovariance()[0]+fParamN.GetCovariance()[0]);
      
      phelix.GetAngle(phase[1][0],mhelix,phase[1][1],angles);
      Double_t tfi  = TMath::Abs(TMath::Tan(angles[0]));
      Double_t tlam = TMath::Abs(TMath::Tan(angles[1]));
      err[0] = err[1]/((0.2+tfi)*(0.2+tfi))+err[2]/((0.2+tlam)*(0.2+tlam));     
      err[0] = ((err[1]*err[2]/((0.2+tfi)*(0.2+tfi)*(0.2+tlam)*(0.2+tlam))))/err[0];
      phelix.ParabolicDCA2(mhelix,phase[1][0],phase[1][1],radius[1],delta2,err);
    }
    Double_t xd[3],xm[3];
    phelix.Evaluate(phase[1][0],xd);
    mhelix.Evaluate(phase[1][1],xm);
    d2 = (xd[0]-xm[0])*(xd[0]-xm[0])+(xd[1]-xm[1])*(xd[1]-xm[1])+(xd[2]-xm[2])*(xd[2]-xm[2]);
  }
  //
  distance2 = TMath::Min(delta1,delta2);
  if (delta1<delta2){
    //get V0 info
    Double_t xd[3],xm[3];
    phelix.Evaluate(phase[0][0],xd);
    mhelix.Evaluate(phase[0][1], xm);
    fPos[0] = 0.5*(xd[0]+xm[0]);
    fPos[1] = 0.5*(xd[1]+xm[1]);
    fPos[2] = 0.5*(xd[2]+xm[2]);

    Float_t wy = fParamP.GetCovariance()[0]/(fParamP.GetCovariance()[0]+fParamN.GetCovariance()[0]);
    Float_t wz = fParamP.GetCovariance()[2]/(fParamP.GetCovariance()[2]+fParamN.GetCovariance()[2]);
    fPos[0] = 0.5*( (1.-wy)*xd[0]+ wy*xm[0] + (1.-wz)*xd[0]+ wz*xm[0] );
    fPos[1] = (1.-wy)*xd[1]+ wy*xm[1];
    fPos[2] = (1.-wz)*xd[2]+ wz*xm[2];
    //
    phelix.GetMomentum(phase[0][0],fPmom);
    mhelix.GetMomentum(phase[0][1],fNmom);
    phelix.GetAngle(phase[0][0],mhelix,phase[0][1],fAngle);
    fRr = TMath::Sqrt(fPos[0]*fPos[0]+fPos[1]*fPos[1]);
  }
  else{
    Double_t xd[3],xm[3];
    phelix.Evaluate(phase[1][0],xd);
    mhelix.Evaluate(phase[1][1], xm);
    fPos[0] = 0.5*(xd[0]+xm[0]);
    fPos[1] = 0.5*(xd[1]+xm[1]);
    fPos[2] = 0.5*(xd[2]+xm[2]);
    Float_t wy = fParamP.GetCovariance()[0]/(fParamP.GetCovariance()[0]+fParamN.GetCovariance()[0]);
    Float_t wz = fParamP.GetCovariance()[2]/(fParamP.GetCovariance()[2]+fParamN.GetCovariance()[2]);
    fPos[0] = 0.5*( (1.-wy)*xd[0]+ wy*xm[0] + (1.-wz)*xd[0]+ wz*xm[0] );
    fPos[1] = (1.-wy)*xd[1]+ wy*xm[1];
    fPos[2] = (1.-wz)*xd[2]+ wz*xm[2];
    //
    phelix.GetMomentum(phase[1][0], fPmom);
    mhelix.GetMomentum(phase[1][1], fNmom);
    phelix.GetAngle(phase[1][0],mhelix,phase[1][1],fAngle);
    fRr = TMath::Sqrt(fPos[0]*fPos[0]+fPos[1]*fPos[1]);
  }
  //Bo:  fDist1 = TMath::Sqrt(TMath::Min(d1,d2));
  //Bo:  fDist2 = TMath::Sqrt(distance2);
  fDcaV0Daughters = TMath::Sqrt(distance2);//Bo:
  //            
  //
  Double_t v[3] = {fPos[0]-vertex[0],fPos[1]-vertex[1],fPos[2]-vertex[2]};
  Double_t p[3] = {fPmom[0]+fNmom[0], fPmom[1]+fNmom[1],fPmom[2]+fNmom[2]};
  Double_t vnorm2 = v[0]*v[0]+v[1]*v[1];
  if (TMath::Abs(v[2])>100000) return;
  Double_t vnorm3 = TMath::Sqrt(TMath::Abs(v[2]*v[2]+vnorm2));
  vnorm2 = TMath::Sqrt(vnorm2);
  Double_t pnorm2 = p[0]*p[0]+p[1]*p[1];
  Double_t pnorm3 = TMath::Sqrt(p[2]*p[2]+pnorm2);
  pnorm2 = TMath::Sqrt(pnorm2);  
  fPointAngleFi = (v[0]*p[0]+v[1]*p[1])/(vnorm2*pnorm2);
  fPointAngleTh = (v[2]*p[2]+vnorm2*pnorm2)/(vnorm3*pnorm3);  
  fPointAngle   = (v[0]*p[0]+v[1]*p[1]+v[2]*p[2])/(vnorm3*pnorm3);
  //
}

