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
#include "AliESDV0MI.h"
#include "AliHelix.h"


ClassImp(AliESDV0MI)

AliESDV0MI::AliESDV0MI(){
  //
  //Dafault constructor
  //
  fID =0;
  fDist1  =-1;
  fDist2  =-1;
  fRr     =-1;
  fStatus = 0;
  fRow0   =-1;  
}

void AliESDV0MI::SetP(const AliExternalTrackParam & paramp)  {
  //
  // set mother
  //
  fParamP   = paramp;
}

void AliESDV0MI::SetM(const AliExternalTrackParam & paramm){
  //
  //set daughter
  //
  fParamM = paramm;

}
  
void  AliESDV0MI::UpdatePID(Double_t pidp[5], Double_t pidm[5])
{
  //
  // set PID hypothesy
  //
  // norm PID to 1
  Float_t sump =0;
  Float_t summ =0;
  for (Int_t i=0;i<5;i++){
    fRP[i]=pidp[i];
    sump+=fRP[i];
    fRM[i]=pidm[i];
    summ+=fRM[i];
  }
  for (Int_t i=0;i<5;i++){
    fRP[i]/=sump;
    fRM[i]/=summ;
  }
}

Float_t AliESDV0MI::GetProb(UInt_t p1, UInt_t p2){
  //
  //
  //
  //
  return TMath::Max(fRP[p1]+fRM[p2], fRP[p2]+fRM[p1]);
}

Float_t AliESDV0MI::GetEffMass(UInt_t p1, UInt_t p2){
  //
  // calculate effective mass
  //
  const Float_t pmass[5] = {5.10000000000000037e-04,1.05660000000000004e-01,1.39570000000000000e-01,
		      4.93599999999999983e-01, 9.38270000000000048e-01};
  if (p1>4) return -1;
  if (p2>4) return -1;
  Float_t mass1 = pmass[p1]; 
  Float_t mass2 = pmass[p2];   
  Double_t *m1 = fPP;
  Double_t *m2 = fPM;
  //
  if (fRP[p1]+fRM[p2]<fRP[p2]+fRM[p1]){
    m1 = fPM;
    m2 = fPP;
  }
  //
  Float_t e1    = TMath::Sqrt(mass1*mass1+
                              m1[0]*m1[0]+
                              m1[1]*m1[1]+
                              m1[2]*m1[2]);
  Float_t e2    = TMath::Sqrt(mass2*mass2+
                              m2[0]*m2[0]+
                              m2[1]*m2[1]+
                              m2[2]*m2[2]);  
  Float_t mass =  
    (m2[0]+m1[0])*(m2[0]+m1[0])+
    (m2[1]+m1[1])*(m2[1]+m1[1])+
    (m2[2]+m1[2])*(m2[2]+m1[2]);
  
  mass = TMath::Sqrt((e1+e2)*(e1+e2)-mass);
  return mass;
}

void  AliESDV0MI::Update(Float_t vertex[3])
{
  //
  // updates Kink Info
  //
  Float_t distance1,distance2;
  //
  AliHelix phelix(fParamP);
  AliHelix mhelix(fParamM);    
  //
  //find intersection linear
  //
  Double_t phase[2][2],radius[2];
  Int_t  points = phelix.GetRPHIintersections(mhelix, phase, radius,200);
  Double_t delta1=10000,delta2=10000;  
  
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
  //
  //find intersection parabolic
  //
  points = phelix.GetRPHIintersections(mhelix, phase, radius);
  delta1=10000,delta2=10000;  
  Double_t d1=1000.,d2=10000.;
  if (points>0){
    phelix.ParabolicDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    phelix.ParabolicDCA(mhelix,phase[0][0],phase[0][1],radius[0],delta1);
    Double_t xd[3],xm[3];
    phelix.Evaluate(phase[0][0],xd);
    mhelix.Evaluate(phase[0][1],xm);
    d1 = (xd[0]-xm[0])*(xd[0]-xm[0])+(xd[1]-xm[1])*(xd[1]-xm[1])+(xd[2]-xm[2])*(xd[2]-xm[2]);
  }
  if (points==2){    
    phelix.ParabolicDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
    phelix.ParabolicDCA(mhelix,phase[1][0],phase[1][1],radius[1],delta2);
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
    fXr[0] = 0.5*(xd[0]+xm[0]);
    fXr[1] = 0.5*(xd[1]+xm[1]);
    fXr[2] = 0.5*(xd[2]+xm[2]);
    //
    phelix.GetMomentum(phase[0][0],fPP);
    mhelix.GetMomentum(phase[0][1],fPM);
    phelix.GetAngle(phase[0][0],mhelix,phase[0][1],fAngle);
    fRr = TMath::Sqrt(fXr[0]*fXr[0]+fXr[1]*fXr[1]);
  }
  else{
    Double_t xd[3],xm[3];
    phelix.Evaluate(phase[1][0],xd);
    mhelix.Evaluate(phase[1][1], xm);
    fXr[0] = 0.5*(xd[0]+xm[0]);
    fXr[1] = 0.5*(xd[1]+xm[1]);
    fXr[2] = 0.5*(xd[2]+xm[2]);
    //
    phelix.GetMomentum(phase[1][0], fPP);
    mhelix.GetMomentum(phase[1][1], fPM);
    phelix.GetAngle(phase[1][0],mhelix,phase[1][1],fAngle);
    fRr = TMath::Sqrt(fXr[0]*fXr[0]+fXr[1]*fXr[1]);
  }
  fDist1 = TMath::Sqrt(TMath::Min(d1,d2));
  fDist2 = TMath::Sqrt(distance2);      
  //            
  //
  Float_t v[3] = {fXr[0]-vertex[0],fXr[1]-vertex[1],fXr[2]-vertex[2]};
  Float_t p[3] = {fPP[0]+fPM[0], fPP[1]+fPM[1],fPP[2]+fPM[2]};
  Float_t vnorm2 = v[0]*v[0]+v[1]*v[1];
  Float_t vnorm3 = TMath::Sqrt(v[2]*v[2]+vnorm2);
  vnorm2 = TMath::Sqrt(vnorm2);
  Float_t pnorm2 = p[0]*p[0]+p[1]*p[1];
  Float_t pnorm3 = TMath::Sqrt(p[2]*p[2]+pnorm2);
  pnorm2 = TMath::Sqrt(pnorm2);  
  fPointAngleFi = (v[0]*p[0]+v[1]*p[1])/(vnorm2*pnorm2);
  fPointAngleTh = (v[2]*p[2]+vnorm2*pnorm2)/(vnorm3*pnorm3);  
  fPointAngle   = (v[0]*p[0]+v[1]*p[1]+v[2]*p[2])/(vnorm3*pnorm3);
  //
}

