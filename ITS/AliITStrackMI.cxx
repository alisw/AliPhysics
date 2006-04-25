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
//                Implementation of the ITS track class
//
//          Origin: Marian Ivanov, CERN, Marian.Ivanov@cern.ch
//     dEdx analysis by: Boris Batyunya, JINR, Boris.Batiounia@cern.ch
//-------------------------------------------------------------------------

#include <TMatrixD.h>

#include <TMath.h>

#include "AliCluster.h"
#include "AliESDtrack.h"
#include "AliITStrackMI.h"

ClassImp(AliITStrackMI)

const Int_t kWARN=5;

//____________________________________________________________________________
AliITStrackMI::AliITStrackMI():AliITStrackV2(),
  fNUsed(0),
  fNSkipped(0),
  fNDeadZone(0),
  fDeadZoneProbability(0),			       
  fReconstructed(kFALSE),
  fConstrain(kFALSE)
{
    for(Int_t i=0; i<kMaxLayer; i++) fClIndex[i]=-1;
    for(Int_t i=0; i<6; i++) { fNy[i]=0; fNz[i]=0; fNormQ[i]=0; fNormChi2[i]=1000;}
    for(Int_t i=0; i<12; i++) {fDy[i]=0; fDz[i]=0; fSigmaY[i]=0; fSigmaZ[i]=0; fChi2MIP[i]=0;}
    fD[0]=0; fD[1]=0;
    fExpQ=40;
    fdEdxMismatch=0;
    fChi22=0;
    fGoldV0 = kFALSE;
}

//____________________________________________________________________________
AliITStrackMI::AliITStrackMI(AliESDtrack& t,Bool_t c) throw (const Char_t *) :
AliITStrackV2(t,c) {
  //------------------------------------------------------------------
  // Conversion ESD track -> ITS track.
  // If c==kTRUE, create the ITS track out of the constrained params.
  //------------------------------------------------------------------
  fNUsed = 0;
  fReconstructed = kFALSE;
  fNSkipped =0; 
  fNDeadZone = 0;
  fDeadZoneProbability = 0;
  for(Int_t i=0; i<6; i++) {fClIndex[i]=-1; fNy[i]=0; fNz[i]=0; fNormQ[i]=0; fNormChi2[i]=1000;}
  for(Int_t i=0; i<12; i++) {fDy[i]=0; fDz[i]=0; fSigmaY[i]=0; fSigmaZ[i]=0;fChi2MIP[i]=0;}
  fD[0]=0; fD[1]=0;
  fExpQ=40;
  fConstrain = kFALSE;
  fdEdxMismatch=0;
  fChi22 =0;
  fGoldV0 = kFALSE;
  //if (!Invariant()) throw "AliITStrackV2: conversion failed !\n";

}

void AliITStrackMI::UpdateESDtrack(ULong_t flags) {
  fESDtrack->UpdateTrackParams(this,flags);
  //if (flags == AliESDtrack::kITSin) fESDtrack->SetITSChi2MIP(fChi2MIP);
}

//____________________________________________________________________________
AliITStrackMI::AliITStrackMI(const AliITStrackMI& t) : AliITStrackV2(t) {
  //------------------------------------------------------------------
  //Copy constructor
  //------------------------------------------------------------------
  fNUsed = t.fNUsed;
  fReconstructed = t.fReconstructed;
  fNSkipped = t.fNSkipped;
  fNDeadZone = t.fNDeadZone;
  fDeadZoneProbability = t.fDeadZoneProbability;
  fLab = t.fLab;
  fFakeRatio = t.fFakeRatio;
  fdEdxMismatch = t.fdEdxMismatch;
  fChi22 = t.fChi22;
  fGoldV0 = t.fGoldV0;;

  fD[0]=t.fD[0]; fD[1]=t.fD[1];
  fDnorm[0] = t.fDnorm[0]; fDnorm[1]=t.fDnorm[1];
  fExpQ= t.fExpQ;
  for(Int_t i=0; i<6; i++) {
    fClIndex[i]= t.fClIndex[i]; fNy[i]=t.fNy[i]; fNz[i]=t.fNz[i]; fNormQ[i]=t.fNormQ[i]; fNormChi2[i] = t.fNormChi2[i];
  }
  for(Int_t i=0; i<12; i++) {fDy[i]=t.fDy[i]; fDz[i]=t.fDz[i]; 
    fSigmaY[i]=t.fSigmaY[i]; fSigmaZ[i]=t.fSigmaZ[i];fChi2MIP[i]=t.fChi2MIP[i];}
  fConstrain = t.fConstrain;
  //memcpy(fDy,t.fDy,6*sizeof(Float_t));
  //memcpy(fDz,t.fDz,6*sizeof(Float_t));
  //memcpy(fSigmaY,t.fSigmaY,6*sizeof(Float_t));
  //memcpy(fSigmaZ,t.fSigmaZ,6*sizeof(Float_t));
  //memcpy(fChi2MIP,t.fChi2MIP,12*sizeof(Float_t));  
}

//_____________________________________________________________________________
Int_t AliITStrackMI::Compare(const TObject *o) const {
  //-----------------------------------------------------------------
  // This function compares tracks according to the their curvature
  //-----------------------------------------------------------------
  AliITStrackMI *t=(AliITStrackMI*)o;
  //Double_t co=TMath::Abs(t->Get1Pt());
  //Double_t c =TMath::Abs(Get1Pt());
  Double_t co=t->GetSigmaY2()*t->GetSigmaZ2()*(0.5+TMath::Sqrt(0.5*t->fD[0]*t->fD[0]+t->fD[1]*t->fD[1]));
  Double_t c =GetSigmaY2()*GetSigmaZ2()*(0.5+TMath::Sqrt(0.5*fD[0]*fD[0]+fD[1]*fD[1]));
  if (c>co) return 1;
  else if (c<co) return -1;
  return 0;
}


Double_t AliITStrackMI::GetPredictedChi2MI(Double_t cy, Double_t cz, Double_t cerry, Double_t cerrz) const
{
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  Double_t r00=cerry*cerry, r01=0., r11=cerrz*cerrz;
  r00+=fC00; r01+=fC10; r11+=fC11;
  //
  Double_t det=r00*r11 - r01*r01;
  if (TMath::Abs(det) < 1.e-30) {
    Int_t n=GetNumberOfClusters();
    if (n>kWARN) 
      Warning("GetPredictedChi2","Singular matrix (%d) !\n",n);
    return 1e10;
  }
  Double_t tmp=r00; r00=r11; r11=tmp; r01=-r01;

  Double_t dy=cy - fP0, dz=cz - fP1;

  return (dy*r00*dy + 2*r01*dy*dz + dz*r11*dz)/det;
}

//____________________________________________________________________________
Int_t AliITStrackMI::CorrectForMaterial(Double_t d, Double_t x0) {
  //------------------------------------------------------------------
  //This function corrects the track parameters for crossed material
  //------------------------------------------------------------------
  //  Double_t p2=(1.+ GetTgl()*GetTgl())/(Get1Pt()*Get1Pt());
  Double_t p2=(1.+ fP3*fP3)/(Get1Pt()*Get1Pt());
  Double_t et   = p2 + GetMass()*GetMass();
  Double_t beta2=p2/et;
  et = sqrt(et);  
  d*=TMath::Sqrt((1.+ fP3*fP3)/(1.- fP2*fP2));
  //d*=TMath::Sqrt(1.+ fP3*fP3 +fP2*fP2/(1.- fP2*fP2));

  //Multiple scattering******************
  if (d!=0) {
     Double_t theta2=14.1*14.1/(beta2*p2*1e6)*TMath::Abs(d);
     //Double_t theta2=1.0259e-6*14*14/28/(beta2*p2)*TMath::Abs(d)*9.36*2.33;
     fC22 += theta2*(1.- fP2*fP2)*(1. + fP3*fP3);
     fC33 += theta2*(1. + fP3*fP3)*(1. + fP3*fP3);
     fC43 += theta2*fP3*fP4*(1. + fP3*fP3);
     fC44 += theta2*fP3*fP4*fP3*fP4;
  }

  //Energy losses************************
  if (x0!=0.) {
     d*=x0;
     //     Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2)) - beta2)*d;
     //Double_t dE=0; 
     Double_t dE = 0.265*0.153e-3*(39.2-55.6*beta2+28.7*beta2*beta2+27.41/beta2)*d;
     /*
     if (beta2/(1-beta2)>3.5*3.5){
       dE=0.153e-3/beta2*(log(3.5*5940)+0.5*log(beta2/(1-beta2)) - beta2)*d;
     }
     else{
       dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2)) - beta2)*d;
       dE+=0.06e-3/(beta2*beta2)*d;
     }
     */
     fP4*=(1.- et/p2*dE);
     Double_t delta44 = (dE*fP4*et/p2);
     delta44*=delta44;
     fC44+= delta44/400.;
  }

  if (!Invariant()) return 0;

  return 1;
}

//____________________________________________________________________________
Int_t AliITStrackMI::UpdateMI(Double_t cy, Double_t cz, Double_t cerry, Double_t cerrz, Double_t chi2,UInt_t index) {
  //------------------------------------------------------------------
  //This function updates track parameters
  //------------------------------------------------------------------
  Double_t p0=fP0,p1=fP1,p2=fP2,p3=fP3,p4=fP4;
  Double_t c00=fC00;
  Double_t c10=fC10, c11=fC11;
  Double_t c20=fC20, c21=fC21, c22=fC22;
  Double_t c30=fC30, c31=fC31, c32=fC32, c33=fC33;
  Double_t c40=fC40, c41=fC41, c42=fC42, c43=fC43, c44=fC44;


  Double_t r00=cerry*cerry, r01=0., r11=cerrz*cerrz;
  r00+=fC00; r01+=fC10; r11+=fC11;
  Double_t det=r00*r11 - r01*r01;
  Double_t tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;

 
  Double_t k00=fC00*r00+fC10*r01, k01=fC00*r01+fC10*r11;
  Double_t k10=fC10*r00+fC11*r01, k11=fC10*r01+fC11*r11;
  Double_t k20=fC20*r00+fC21*r01, k21=fC20*r01+fC21*r11;
  Double_t k30=fC30*r00+fC31*r01, k31=fC30*r01+fC31*r11;
  Double_t k40=fC40*r00+fC41*r01, k41=fC40*r01+fC41*r11;

  Double_t dy=cy - fP0, dz=cz - fP1;
  Int_t layer = (index & 0xf0000000) >> 28;
  fDy[layer] = dy;
  fDz[layer] = dz;
  fSigmaY[layer] = TMath::Sqrt(cerry*cerry+fC00);
  fSigmaZ[layer] = TMath::Sqrt(cerrz*cerrz+fC11);

  Double_t sf=fP2 + k20*dy + k21*dz;
  
  fP0 += k00*dy + k01*dz;
  fP1 += k10*dy + k11*dz;
  fP2  = sf;
  fP3 += k30*dy + k31*dz;
  fP4 += k40*dy + k41*dz;
  
  Double_t c01=fC10, c02=fC20, c03=fC30, c04=fC40;
  Double_t c12=fC21, c13=fC31, c14=fC41;

  fC00-=k00*fC00+k01*fC10; fC10-=k00*c01+k01*fC11;
  fC20-=k00*c02+k01*c12;   fC30-=k00*c03+k01*c13;
  fC40-=k00*c04+k01*c14; 

  fC11-=k10*c01+k11*fC11;
  fC21-=k10*c02+k11*c12;   fC31-=k10*c03+k11*c13;
  fC41-=k10*c04+k11*c14; 

  fC22-=k20*c02+k21*c12;   fC32-=k20*c03+k21*c13;
  fC42-=k20*c04+k21*c14; 

  fC33-=k30*c03+k31*c13;
  fC43-=k30*c04+k31*c14; 

  fC44-=k40*c04+k41*c14; 

  if (!Invariant()) {
     fP0=p0; fP1=p1; fP2=p2; fP3=p3; fP4=p4;
     fC00=c00;
     fC10=c10; fC11=c11;
     fC20=c20; fC21=c21; fC22=c22;
     fC30=c30; fC31=c31; fC32=c32; fC33=c33;
     fC40=c40; fC41=c41; fC42=c42; fC43=c43; fC44=c44;
     return 0;
  }

  if (chi2<0) return 1;
  Int_t n=GetNumberOfClusters();
  fIndex[n]=index;
  SetNumberOfClusters(n+1);
  SetChi2(GetChi2()+chi2);

  return 1;
}

Int_t AliITStrackMI::GetProlongationFast(Double_t alp, Double_t xk,Double_t &y, Double_t &z)
{
  //-----------------------------------------------------------------------------
  //get fast prolongation 
  //-----------------------------------------------------------------------------
  Double_t ca=TMath::Cos(alp-fAlpha), sa=TMath::Sin(alp-fAlpha);
  Double_t cf=TMath::Sqrt(1.- fP2*fP2);  
  // **** rotation **********************  
  y= -fX*sa + fP0*ca;
  // **** translation ******************  
  Double_t dx = xk- fX*ca - fP0*sa;
  Double_t f1=fP2*ca - cf*sa, f2=f1 + fP4*dx;
  if (TMath::Abs(f2) >= 0.9999) {
    return 0;
  }
  Double_t r1=TMath::Sqrt(1.- f1*f1), r2=TMath::Sqrt(1.- f2*f2);  
  y += dx*(f1+f2)/(r1+r2);
  z  = fP1+dx*(f1+f2)/(f1*r2 + f2*r1)*fP3;  
  return 1;
}


Bool_t AliITStrackMI::IsGoldPrimary()
{
  //
  // Indicates gold pimary track
  //
  Bool_t isGold=kTRUE;
  if (!fConstrain) return kFALSE;                // 
  if (fNDeadZone+fNDeadZone<5.5) isGold =  kFALSE; // short track
  //
  if (fChi2/Float_t(fN)>2.){
    if (fChi2MIP[0]+fNUsed>3.5) isGold = kFALSE;    
  }
  if (fChi2MIP[2]>4.5) isGold = kFALSE;         //back propagation chi2
  //
  if (fDnorm[0]>0&&fDnorm[1]>0){
    const Float_t distcut2 =2.5*2.5;  //normalize distance  cut 
    Float_t dist2 = fD[0]*fD[0]/(fDnorm[0]*fDnorm[0])+fD[1]*fD[1]/(fDnorm[1]*fDnorm[1]);  //normalize distance to the vertex (pools)
    if (dist2>distcut2) isGold = kFALSE;
  }
  return isGold;
}
