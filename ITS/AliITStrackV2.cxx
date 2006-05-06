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

///////////////////////////////////////////////////////////////////////////
//                Implementation of the ITS track class
//
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//     dEdx analysis by: Boris Batyunya, JINR, Boris.Batiounia@cern.ch
///////////////////////////////////////////////////////////////////////////
#include <TMath.h>

#include "AliCluster.h"
#include "AliESDtrack.h"
#include "AliITStrackV2.h"
#include "AliStrLine.h"

ClassImp(AliITStrackV2)

const Int_t kWARN=5;

//____________________________________________________________________________
AliITStrackV2::AliITStrackV2():AliKalmanTrack(),
  fX(0),
  fAlpha(0),
  fdEdx(0),
  fP0(0),
  fP1(0),
  fP2(0),
  fP3(0),
  fP4(0),
  fC00(0),
  fC10(0),
  fC11(0),
  fC20(0),
  fC21(0),
  fC22(0),
  fC30(0),
  fC31(0),
  fC32(0),
  fC33(0),
  fC40(0),
  fC41(0),
  fC42(0),
  fC43(0),
  fC44(0),
  fESDtrack(0)
{
    for(Int_t i=0; i<2*kMaxLayer; i++) fIndex[i]=-1;
    for(Int_t i=0; i<4; i++) fdEdxSample[i]=0;
}


//____________________________________________________________________________
AliITStrackV2::AliITStrackV2(AliESDtrack& t,Bool_t c) throw (const Char_t *) :
AliKalmanTrack() {
  //------------------------------------------------------------------
  // Conversion ESD track -> ITS track.
  // If c==kTRUE, create the ITS track out of the constrained params.
  //------------------------------------------------------------------
  SetNumberOfClusters(t.GetITSclusters(fIndex));
  SetLabel(t.GetLabel());
  SetMass(t.GetMass());
  //
  //

  fdEdx=t.GetITSsignal();
  fAlpha = t.GetAlpha();
  if      (fAlpha < -TMath::Pi()) fAlpha += 2*TMath::Pi();
  else if (fAlpha >= TMath::Pi()) fAlpha -= 2*TMath::Pi();

  //Conversion of the track parameters
  Double_t x,p[5]; 
  if (c) t.GetConstrainedExternalParameters(fAlpha,x,p);
  else t.GetExternalParameters(x,p);
  fX=x;   
  fP0=p[0]; 
  fP1=p[1];   SaveLocalConvConst(); 
  fP2=p[2];
  fP3=p[3];   x=GetLocalConvConst();
  fP4=p[4]/x; 

  //Conversion of the covariance matrix
  Double_t cv[15]; 
  if (c) t.GetConstrainedExternalCovariance(cv);
  else t.GetExternalCovariance(cv);
  fC00=cv[0 ];
  fC10=cv[1 ];   fC11=cv[2 ];
  fC20=cv[3 ];   fC21=cv[4 ];   fC22=cv[5 ];
  fC30=cv[6 ];   fC31=cv[7 ];   fC32=cv[8 ];   fC33=cv[9 ];
  fC40=cv[10]/x; fC41=cv[11]/x; fC42=cv[12]/x; fC43=cv[13]/x; fC44=cv[14]/x/x;

  if (t.GetStatus()&AliESDtrack::kTIME) {
    StartTimeIntegral();
    Double_t times[10]; t.GetIntegratedTimes(times); SetIntegratedTimes(times);
    SetIntegratedLength(t.GetIntegratedLength());
  }
  fESDtrack=&t;

  //  if (!Invariant()) throw "AliITStrackV2: conversion failed !\n";
  for(Int_t i=0; i<4; i++) fdEdxSample[i]=0;
}

void AliITStrackV2::UpdateESDtrack(ULong_t flags) const {
  fESDtrack->UpdateTrackParams(this,flags);
}

//____________________________________________________________________________
AliITStrackV2::AliITStrackV2(const AliITStrackV2& t) : AliKalmanTrack(t) {
  //------------------------------------------------------------------
  //Copy constructor
  //------------------------------------------------------------------
  fX=t.fX;
  fAlpha=t.fAlpha;
  fdEdx=t.fdEdx;

  fP0=t.fP0; fP1=t.fP1; fP2=t.fP2; fP3=t.fP3; fP4=t.fP4;

  fC00=t.fC00;
  fC10=t.fC10;  fC11=t.fC11;
  fC20=t.fC20;  fC21=t.fC21;  fC22=t.fC22;
  fC30=t.fC30;  fC31=t.fC31;  fC32=t.fC32;  fC33=t.fC33;
  fC40=t.fC40;  fC41=t.fC41;  fC42=t.fC42;  fC43=t.fC43;  fC44=t.fC44;

  Int_t i;
  for (i=0; i<2*kMaxLayer; i++) fIndex[i]=t.fIndex[i];
  for (i=0; i<4; i++) fdEdxSample[i]=t.fdEdxSample[i];

  fESDtrack=t.fESDtrack;
}

//_____________________________________________________________________________
Int_t AliITStrackV2::Compare(const TObject *o) const {
  //-----------------------------------------------------------------
  // This function compares tracks according to the their curvature
  //-----------------------------------------------------------------
  AliITStrackV2 *t=(AliITStrackV2*)o;
  //Double_t co=TMath::Abs(t->Get1Pt());
  //Double_t c =TMath::Abs(Get1Pt());
  Double_t co=t->GetSigmaY2()*t->GetSigmaZ2();
  Double_t c =GetSigmaY2()*GetSigmaZ2();
  if (c>co) return 1;
  else if (c<co) return -1;
  return 0;
}

//_____________________________________________________________________________
void AliITStrackV2::GetExternalCovariance(Double_t cc[15]) const {
  //-------------------------------------------------------------------------
  // This function returns an external representation of the covriance matrix.
  //   (See comments in AliTPCtrack.h about external track representation)
  //-------------------------------------------------------------------------
  Double_t a=GetLocalConvConst();

  cc[0 ]=fC00;
  cc[1 ]=fC10;   cc[2 ]=fC11;
  cc[3 ]=fC20;   cc[4 ]=fC21;   cc[5 ]=fC22;
  cc[6 ]=fC30;   cc[7 ]=fC31;   cc[8 ]=fC32;   cc[9 ]=fC33;
  cc[10]=fC40*a; cc[11]=fC41*a; cc[12]=fC42*a; cc[13]=fC43*a; cc[14]=fC44*a*a;
}

//____________________________________________________________________________
Int_t AliITStrackV2::PropagateToVertex(Double_t d,Double_t x0) {
  //------------------------------------------------------------------
  //This function propagates a track to the minimal distance from the origin
  //------------------------------------------------------------------
  //Double_t xv=fP2*(fX*fP2 - fP0*TMath::Sqrt(1.- fP2*fP2)); //linear approxim.
  Double_t tgf=-(fP4*fX - fP2)/(fP4*fP0 + TMath::Sqrt(1 - fP2*fP2));
  Double_t snf=tgf/TMath::Sqrt(1.+ tgf*tgf);
  Double_t xv=(snf - fP2)/fP4 + fX;
  return PropagateTo(xv,d,x0);
}

//____________________________________________________________________________
Int_t AliITStrackV2::
GetGlobalXYZat(Double_t xk, Double_t &x, Double_t &y, Double_t &z) const {
  //------------------------------------------------------------------
  //This function returns a track position in the global system
  //------------------------------------------------------------------
  Double_t dx=xk-fX;
  Double_t f1=fP2, f2=f1 + fP4*dx;
  if (TMath::Abs(f2) >= 0.9999) {
    Int_t n=GetNumberOfClusters();
    if (n>kWARN) 
      Warning("GetGlobalXYZat","Propagation failed (%d) !\n",n);
    return 0;
  }

  Double_t r1=sqrt(1.- f1*f1), r2=sqrt(1.- f2*f2);
  
  Double_t yk = fP0 + dx*(f1+f2)/(r1+r2);
  Double_t zk = fP1 + dx*(f1+f2)/(f1*r2 + f2*r1)*fP3;

  Double_t cs=TMath::Cos(fAlpha), sn=TMath::Sin(fAlpha);
  x = xk*cs - yk*sn;
  y = xk*sn + yk*cs;
  z = zk;

  return 1;
}

//_____________________________________________________________________________
void AliITStrackV2::ApproximateHelixWithLine(Double_t xk, AliStrLine *line)
{
  //------------------------------------------------------------
  // Approximate the track (helix) with a straight line tangent to the
  // helix in the point defined by r (F. Prino, prino@to.infn.it)
  //------------------------------------------------------------
  Double_t mom[3];
  Double_t azim = TMath::ASin(fP2)+fAlpha;
  Double_t theta = TMath::Pi()/2. - TMath::ATan(fP3);
  mom[0] = TMath::Sin(theta)*TMath::Cos(azim);
  mom[1] = TMath::Sin(theta)*TMath::Sin(azim);
  mom[2] = TMath::Cos(theta);
  Double_t pos[3];
  GetGlobalXYZat(xk,pos[0],pos[1],pos[2]);
  line->SetP0(pos);
  line->SetCd(mom);
}
//_____________________________________________________________________________
Double_t AliITStrackV2::GetPredictedChi2(const AliCluster *c) const 
{
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  Double_t r00=c->GetSigmaY2(), r01=0., r11=c->GetSigmaZ2();
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

  Double_t dy=c->GetY() - fP0, dz=c->GetZ() - fP1;

  return (dy*r00*dy + 2*r01*dy*dz + dz*r11*dz)/det;
}

//____________________________________________________________________________
Int_t AliITStrackV2::CorrectForMaterial(Double_t d, Double_t x0) {
  //------------------------------------------------------------------
  //This function corrects the track parameters for crossed material
  //------------------------------------------------------------------
  Double_t p2=(1.+ fP3*fP3)/(Get1Pt()*Get1Pt());
  Double_t beta2=p2/(p2 + GetMass()*GetMass());
  d*=TMath::Sqrt((1.+ fP3*fP3)/(1.- fP2*fP2));

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
     Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2)) - beta2)*d;
     if (beta2/(1-beta2)>3.5*3.5)
       dE=0.153e-3/beta2*(log(3.5*5940)+0.5*log(beta2/(1-beta2)) - beta2)*d;

     fP4*=(1.- TMath::Sqrt(p2+GetMass()*GetMass())/p2*dE);
  }

  if (!Invariant()) return 0;

  return 1;
}

//____________________________________________________________________________
Int_t AliITStrackV2::PropagateTo(Double_t xk, Double_t d, Double_t x0) {
  //------------------------------------------------------------------
  //This function propagates a track
  //------------------------------------------------------------------
  Double_t x1=fX, x2=xk, dx=x2-x1;
  Double_t f1=fP2, f2=f1 + fP4*dx;
  if (TMath::Abs(f2) >= 0.98) {
    // MI change  - don't propagate highly inclined tracks
    //              covariance matrix distorted
    //Int_t n=GetNumberOfClusters();
    //if (n>kWARN) 
    //   Warning("PropagateTo","Propagation failed !\n",n);
    return 0;
  }
  Double_t lcc=GetLocalConvConst();  

  // old position [SR, GSI, 17.02.2003]
  Double_t oldX = fX, oldY = fP0, oldZ = fP1;

  Double_t r1=sqrt(1.- f1*f1), r2=sqrt(1.- f2*f2);
  
  fP0 += dx*(f1+f2)/(r1+r2);
  fP1 += dx*(f1+f2)/(f1*r2 + f2*r1)*fP3;
  fP2 += dx*fP4;

  //f = F - 1
  
  Double_t f02=    dx/(r1*r1*r1);
  Double_t f04=0.5*dx*dx/(r1*r1*r1);
  Double_t f12=    dx*fP3*f1/(r1*r1*r1);
  Double_t f14=0.5*dx*dx*fP3*f1/(r1*r1*r1);
  Double_t f13=    dx/r1;
  Double_t f24=    dx; 
  
  //b = C*ft
  Double_t b00=f02*fC20 + f04*fC40, b01=f12*fC20 + f14*fC40 + f13*fC30;
  Double_t b02=f24*fC40;
  Double_t b10=f02*fC21 + f04*fC41, b11=f12*fC21 + f14*fC41 + f13*fC31;
  Double_t b12=f24*fC41;
  Double_t b20=f02*fC22 + f04*fC42, b21=f12*fC22 + f14*fC42 + f13*fC32;
  Double_t b22=f24*fC42;
  Double_t b40=f02*fC42 + f04*fC44, b41=f12*fC42 + f14*fC44 + f13*fC43;
  Double_t b42=f24*fC44;
  Double_t b30=f02*fC32 + f04*fC43, b31=f12*fC32 + f14*fC43 + f13*fC33;
  Double_t b32=f24*fC43;
  
  //a = f*b = f*C*ft
  Double_t a00=f02*b20+f04*b40,a01=f02*b21+f04*b41,a02=f02*b22+f04*b42;
  Double_t a11=f12*b21+f14*b41+f13*b31,a12=f12*b22+f14*b42+f13*b32;
  Double_t a22=f24*b42;

  //F*C*Ft = C + (b + bt + a)
  fC00 += b00 + b00 + a00;
  fC10 += b10 + b01 + a01; 
  fC20 += b20 + b02 + a02;
  fC30 += b30;
  fC40 += b40;
  fC11 += b11 + b11 + a11;
  fC21 += b21 + b12 + a12;
  fC31 += b31; 
  fC41 += b41;
  fC22 += b22 + b22 + a22;
  fC32 += b32;
  fC42 += b42;

  fX=x2;

  //Change of the magnetic field *************
  SaveLocalConvConst();
  fP4*=lcc/GetLocalConvConst();

  if (!CorrectForMaterial(d,x0)) return 0;

  // Integrated Time [SR, GSI, 17.02.2003]
  if (IsStartedTimeIntegral() && fX>oldX) {
    Double_t l2 = (fX-oldX)*(fX-oldX)+(fP0-oldY)*(fP0-oldY)+
                  (fP1-oldZ)*(fP1-oldZ);
    AddTimeStep(TMath::Sqrt(l2));
  }
  //

  return 1;
}

//____________________________________________________________________________
Int_t AliITStrackV2::Update(const AliCluster* c, Double_t chi2, UInt_t index) {
  //------------------------------------------------------------------
  //This function updates track parameters
  //------------------------------------------------------------------
  Double_t p0=fP0,p1=fP1,p2=fP2,p3=fP3,p4=fP4;
  Double_t c00=fC00;
  Double_t c10=fC10, c11=fC11;
  Double_t c20=fC20, c21=fC21, c22=fC22;
  Double_t c30=fC30, c31=fC31, c32=fC32, c33=fC33;
  Double_t c40=fC40, c41=fC41, c42=fC42, c43=fC43, c44=fC44;


  Double_t r00=c->GetSigmaY2(), r01=0., r11=c->GetSigmaZ2();
  r00+=fC00; r01+=fC10; r11+=fC11;
  Double_t det=r00*r11 - r01*r01;
  Double_t tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;

 
  Double_t k00=fC00*r00+fC10*r01, k01=fC00*r01+fC10*r11;
  Double_t k10=fC10*r00+fC11*r01, k11=fC10*r01+fC11*r11;
  Double_t k20=fC20*r00+fC21*r01, k21=fC20*r01+fC21*r11;
  Double_t k30=fC30*r00+fC31*r01, k31=fC30*r01+fC31*r11;
  Double_t k40=fC40*r00+fC41*r01, k41=fC40*r01+fC41*r11;

  Double_t dy=c->GetY() - fP0, dz=c->GetZ() - fP1;
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

Int_t AliITStrackV2::Invariant() const {
  //------------------------------------------------------------------
  // This function is for debugging purpose only
  //------------------------------------------------------------------
  Int_t n=GetNumberOfClusters();
  
  if (TMath::Abs(fP2)>=0.9999){
     if (n>kWARN) Warning("Invariant","fP2=%f\n",fP2);
     return 0;
  }
  if (fC00<=0 || fC00>9.) {
     if (n>kWARN) Warning("Invariant","fC00=%f\n",fC00); 
     return 0;
  }
  if (fC11<=0 || fC11>9.) {
     if (n>kWARN) Warning("Invariant","fC11=%f\n",fC11); 
     return 0;
  }
  if (fC22<=0 || fC22>1.) {
     if (n>kWARN) Warning("Invariant","fC22=%f\n",fC22); 
     return 0;
  }
  if (fC33<=0 || fC33>1.) {
     if (n>kWARN) Warning("Invariant","fC33=%f\n",fC33); 
     return 0;
  }
  if (fC44<=0 || fC44>6e-5) {
     if (n>kWARN) Warning("Invariant","fC44=%f\n",fC44);
     return 0;
  }
  return 1;
}

//____________________________________________________________________________
Int_t AliITStrackV2::Propagate(Double_t alp,Double_t xk) {
  //------------------------------------------------------------------
  //This function propagates a track
  //------------------------------------------------------------------
  Double_t alpha=fAlpha, x=fX;
  Double_t p0=fP0,p1=fP1,p2=fP2,p3=fP3,p4=fP4;
  Double_t c00=fC00;
  Double_t c10=fC10, c11=fC11;
  Double_t c20=fC20, c21=fC21, c22=fC22;
  Double_t c30=fC30, c31=fC31, c32=fC32, c33=fC33;
  Double_t c40=fC40, c41=fC41, c42=fC42, c43=fC43, c44=fC44;

  if      (alp < -TMath::Pi()) alp += 2*TMath::Pi();
  else if (alp >= TMath::Pi()) alp -= 2*TMath::Pi();
  Double_t ca=TMath::Cos(alp-fAlpha), sa=TMath::Sin(alp-fAlpha);
  Double_t sf=fP2, cf=TMath::Sqrt(1.- fP2*fP2);

  // **** rotation **********************
  {
  fAlpha = alp;
  fX =  x*ca + p0*sa;
  fP0= -x*sa + p0*ca;
  fP2=  sf*ca - cf*sa;

  Double_t rr=(ca+sf/cf*sa);  

  fC00 *= (ca*ca);
  fC10 *= ca; 
  fC20 *= ca*rr;
  fC30 *= ca;
  fC40 *= ca;
  //fC11 = fC11;
  fC21 *= rr;
  //fC31 = fC31; 
  //fC41 = fC41;
  fC22 *= rr*rr;
  fC32 *= rr;
  fC42 *= rr;
  //fC33=fC33;
  //fC43=fC43;
  //fC44=fC44;
 
  }

  // **** translation ******************
  {
  Double_t dx=xk-fX;
  Double_t f1=fP2, f2=f1 + fP4*dx;
  if (TMath::Abs(f2) >= 0.98) {
    // don't propagate highly inclined tracks MI
    return 0;
  }
  //    Int_t n=GetNumberOfClusters();
  //  if (n>kWARN) 
  //     Warning("Propagate","Propagation failed (%d) !\n",n);
  //  return 0;
  //}
  Double_t lcc=GetLocalConvConst();  

  Double_t r1=TMath::Sqrt(1.- f1*f1), r2=TMath::Sqrt(1.- f2*f2);

  fX=xk;
  fP0 += dx*(f1+f2)/(r1+r2);
  fP1 += dx*(f1+f2)/(f1*r2 + f2*r1)*fP3;
  fP2 += dx*fP4;

  //Change of the magnetic field *************
  SaveLocalConvConst();
  fP4*=lcc/GetLocalConvConst();

  //f = F - 1
  
  Double_t f02=    dx/(r1*r1*r1);
  Double_t f04=0.5*dx*dx/(r1*r1*r1);
  Double_t f12=    dx*fP3*f1/(r1*r1*r1);
  Double_t f14=0.5*dx*dx*fP3*f1/(r1*r1*r1);
  Double_t f13=    dx/r1;
  Double_t f24=    dx; 
  
  //b = C*ft
  Double_t b00=f02*fC20 + f04*fC40, b01=f12*fC20 + f14*fC40 + f13*fC30;
  Double_t b02=f24*fC40;
  Double_t b10=f02*fC21 + f04*fC41, b11=f12*fC21 + f14*fC41 + f13*fC31;
  Double_t b12=f24*fC41;
  Double_t b20=f02*fC22 + f04*fC42, b21=f12*fC22 + f14*fC42 + f13*fC32;
  Double_t b22=f24*fC42;
  Double_t b40=f02*fC42 + f04*fC44, b41=f12*fC42 + f14*fC44 + f13*fC43;
  Double_t b42=f24*fC44;
  Double_t b30=f02*fC32 + f04*fC43, b31=f12*fC32 + f14*fC43 + f13*fC33;
  Double_t b32=f24*fC43;
  
  //a = f*b = f*C*ft
  Double_t a00=f02*b20+f04*b40,a01=f02*b21+f04*b41,a02=f02*b22+f04*b42;
  Double_t a11=f12*b21+f14*b41+f13*b31,a12=f12*b22+f14*b42+f13*b32;
  Double_t a22=f24*b42;

  //F*C*Ft = C + (b + bt + a)
  fC00 += b00 + b00 + a00;
  fC10 += b10 + b01 + a01; 
  fC20 += b20 + b02 + a02;
  fC30 += b30;
  fC40 += b40;
  fC11 += b11 + b11 + a11;
  fC21 += b21 + b12 + a12;
  fC31 += b31; 
  fC41 += b41;
  fC22 += b22 + b22 + a22;
  fC32 += b32;
  fC42 += b42;

  if (!Invariant()) {
     fAlpha=alpha; 
     fX=x; 
     fP0=p0; fP1=p1; fP2=p2; fP3=p3; fP4=p4;
     fC00=c00;
     fC10=c10; fC11=c11;
     fC20=c20; fC21=c21; fC22=c22;
     fC30=c30; fC31=c31; fC32=c32; fC33=c33;
     fC40=c40; fC41=c41; fC42=c42; fC43=c43; fC44=c44;
     return 0;
  }
  }

  return 1;
}


Double_t AliITStrackV2::GetD(Double_t x, Double_t y) const {
  //------------------------------------------------------------------
  // This function calculates the transverse impact parameter
  // with respect to a point with global coordinates (x,y)
  //------------------------------------------------------------------
  Double_t xt=fX, yt=fP0;

  Double_t sn=TMath::Sin(fAlpha), cs=TMath::Cos(fAlpha);
  Double_t a = x*cs + y*sn;
  y = -x*sn + y*cs; x=a;
  xt-=x; yt-=y;

  sn=fP4*xt - fP2; cs=fP4*yt + TMath::Sqrt(1.- fP2*fP2);
  a=2*(xt*fP2 - yt*TMath::Sqrt(1.- fP2*fP2))-fP4*(xt*xt + yt*yt);
  if (fP4<0) a=-a;
  return a/(1 + TMath::Sqrt(sn*sn + cs*cs));
}

Double_t AliITStrackV2::GetZat(Double_t x) const {
  //------------------------------------------------------------------
  // This function calculates the z at given x point - in current coordinate system
  //------------------------------------------------------------------
  Double_t x1=fX, x2=x, dx=x2-x1;
  //
  Double_t f1=fP2, f2=f1 + fP4*dx;
  if (TMath::Abs(f2) >= 0.9999) {
    return 10000000;
  }
  Double_t r1=sqrt(1.- f1*f1), r2=sqrt(1.- f2*f2);
  Double_t z =  fP1 + dx*(f1+f2)/(f1*r2 + f2*r1)*fP3;
  return z;
}




Int_t AliITStrackV2::Improve(Double_t x0,Double_t xyz[3],Double_t ers[3]) {
  //------------------------------------------------------------------
  //This function improves angular track parameters  
  //------------------------------------------------------------------
  Double_t cs=TMath::Cos(fAlpha), sn=TMath::Sin(fAlpha);
  //Double_t xv = xyz[0]*cs + xyz[1]*sn; // vertex
    Double_t yv =-xyz[0]*sn + xyz[1]*cs; // in the
    Double_t zv = xyz[2];                // local frame
  Double_t dy=fP0-yv, dz=fP1-zv;
  Double_t r2=fX*fX+dy*dy;
  Double_t p2=(1.+ GetTgl()*GetTgl())/(Get1Pt()*Get1Pt());
  Double_t beta2=p2/(p2 + GetMass()*GetMass());
  x0*=TMath::Sqrt((1.+ GetTgl()*GetTgl())/(1.- GetSnp()*GetSnp()));
  Double_t theta2=14.1*14.1/(beta2*p2*1e6)*x0;
  //Double_t theta2=1.0259e-6*14*14/28/(beta2*p2)*x0*9.36*2.33;
  {
  Double_t dummy=4/r2-fP4*fP4;
  if (dummy < 0) return 0;
  Double_t parp=0.5*(fP4*fX + dy*TMath::Sqrt(dummy));
  Double_t sigma2p = theta2*(1.- GetSnp()*GetSnp())*(1. + GetTgl()*GetTgl());
  sigma2p += fC00/r2*(1.- dy*dy/r2)*(1.- dy*dy/r2);
  sigma2p += ers[1]*ers[1]/r2;
  sigma2p += 0.25*fC44*fX*fX;
  Double_t eps2p=sigma2p/(fC22+sigma2p);
  fP0 += fC20/(fC22+sigma2p)*(parp-fP2);
  fP2 = eps2p*fP2 + (1-eps2p)*parp;
  fC22 *= eps2p;
  fC20 *= eps2p;
  }
  {
  Double_t parl=0.5*fP4*dz/TMath::ASin(0.5*fP4*TMath::Sqrt(r2));
  Double_t sigma2l=theta2;
  sigma2l += fC11/r2+fC00*dy*dy*dz*dz/(r2*r2*r2);
  sigma2l += ers[2]*ers[2]/r2;
  Double_t eps2l=sigma2l/(fC33+sigma2l);
  fP1 += fC31/(fC33+sigma2l)*(parl-fP3);
  fP4 += fC43/(fC33+sigma2l)*(parl-fP3);
  fP3 = eps2l*fP3 + (1-eps2l)*parl;
  fC33 *= eps2l; fC43 *= eps2l; 
  fC31 *= eps2l; 
  }
  if (!Invariant()) return 0;
  return 1;
} 

void AliITStrackV2::ResetCovariance() {
  //------------------------------------------------------------------
  //This function makes a track forget its history :)  
  //------------------------------------------------------------------

  fC00*=10.;
  fC10=0.;  fC11*=10.;
  fC20=0.;  fC21=0.;  fC22*=10.;
  fC30=0.;  fC31=0.;  fC32=0.;  fC33*=10.;
  fC40=0.;  fC41=0.;  fC42=0.;  fC43=0.;  fC44*=10.;

}

void AliITStrackV2::CookdEdx(Double_t low, Double_t up) {
  //-----------------------------------------------------------------
  // This function calculates dE/dX within the "low" and "up" cuts.
  // Origin: Boris Batyunya, JINR, Boris.Batiounia@cern.ch 
  //-----------------------------------------------------------------
  // The clusters order is: SSD-2, SSD-1, SDD-2, SDD-1, SPD-2, SPD-1

  Int_t i;
  Int_t nc=0;
  for (i=0; i<GetNumberOfClusters(); i++) {
    Int_t idx=GetClusterIndex(i);
    idx=(idx&0xf0000000)>>28;
    if (idx>1) nc++; // Take only SSD and SDD
  }

  Int_t swap;//stupid sorting
  do {
    swap=0;
    for (i=0; i<nc-1; i++) {
      if (fdEdxSample[i]<=fdEdxSample[i+1]) continue;
      Float_t tmp=fdEdxSample[i];
      fdEdxSample[i]=fdEdxSample[i+1]; fdEdxSample[i+1]=tmp;
      swap++;
    }
  } while (swap);

  Int_t nl=Int_t(low*nc), nu=Int_t(up*nc); //b.b. to take two lowest dEdX
                                           // values from four ones choose
                                           // nu=2
  Float_t dedx=0;
  for (i=nl; i<nu; i++) dedx += fdEdxSample[i];
  if (nu-nl>0) dedx /= (nu-nl);

  SetdEdx(dedx);
}

Double_t AliITStrackV2::
PropagateToDCA(AliKalmanTrack *p, Double_t d, Double_t x0) {
  //--------------------------------------------------------------
  // Propagates this track and the argument track to the position of the
  // distance of closest approach. 
  // Returns the (weighed !) distance of closest approach.
  //--------------------------------------------------------------
  Double_t xthis, xp, dca;
  {
  //Temporary solution
  Double_t b=1./GetLocalConvConst()/kB2C;
  AliExternalTrackParam dummy1(*this), dummy2(*p); 
  dca=dummy1.GetDCA(&dummy2,b,xthis,xp);
  }
  if (!PropagateTo(xthis,d,x0)) {
    //AliWarning(" propagation failed !");
    return 1e+33;
  }  

  if (!p->PropagateTo(xp,d,x0)) {
    //AliWarning(" propagation failed !";
    return 1e+33;
  }  

  return dca;
} 

Double_t AliITStrackV2::Get1Pt() const {
  //--------------------------------------------------------------
  // Returns the inverse Pt (1/GeV/c)
  // (or 1/"most probable pt", if the field is too weak)
  //--------------------------------------------------------------
  if (TMath::Abs(GetLocalConvConst()) > kVeryBigConvConst)
      return 1./kMostProbableMomentum/TMath::Sqrt(1.+ GetTgl()*GetTgl());
  return (TMath::Sign(1e-9,fP4) + fP4)*GetLocalConvConst();
}
