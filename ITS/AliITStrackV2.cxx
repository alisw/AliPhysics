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
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------

#include <TMatrixD.h>

#include <TMath.h>
#include <iostream.h>

#include "AliCluster.h"
#include "../TPC/AliTPCtrack.h"
#include "AliITStrackV2.h"

ClassImp(AliITStrackV2)

const Int_t kWARN=1;

//____________________________________________________________________________
AliITStrackV2::AliITStrackV2(const AliTPCtrack& t) throw (const Char_t *) {
  //------------------------------------------------------------------
  //Convertion TPC track -> ITS track
  //------------------------------------------------------------------
  SetLabel(t.GetLabel());
  SetChi2(0.);
  SetNumberOfClusters(0);
  fdEdx  = 0.;
  fAlpha = t.GetAlpha();
  if      (fAlpha < -TMath::Pi()) fAlpha += 2*TMath::Pi();
  else if (fAlpha >= TMath::Pi()) fAlpha -= 2*TMath::Pi();

  //Convertion of the track parameters
  Double_t x,p[5]; t.GetExternalParameters(x,p);
  fX=x;    x=kConversionConstant;
  fP0=p[0]; 
  fP1=p[1]; 
  fP2=p[2];
  fP3=p[3];
  fP4=p[4]/x; 

  //Convertion of the covariance matrix
  Double_t c[15]; t.GetExternalCovariance(c);

  fC00=c[0 ];
  fC10=c[1 ];   fC11=c[2 ];
  fC20=c[3 ];   fC21=c[4 ];   fC22=c[5 ];
  fC30=c[6 ];   fC31=c[7 ];   fC32=c[8 ];   fC33=c[9 ];
  fC40=c[10]/x; fC41=c[11]/x; fC42=c[12]/x; fC43=c[13]/x; fC44=c[14]/x/x;

  if (!Invariant()) throw "AliITStrackV2: conversion failed !\n";
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

  Int_t n=GetNumberOfClusters();
  for (Int_t i=0; i<n; i++) fIndex[i]=t.fIndex[i];
}

//_____________________________________________________________________________
Int_t AliITStrackV2::Compare(const TObject *o) const {
  //-----------------------------------------------------------------
  // This function compares tracks according to the their curvature
  //-----------------------------------------------------------------
  AliTPCtrack *t=(AliTPCtrack*)o;
  Double_t co=TMath::Abs(t->Get1Pt());
  Double_t c =TMath::Abs(Get1Pt());
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
  Double_t a=kConvConst;

  cc[0 ]=fC00;
  cc[1 ]=fC10;   cc[2 ]=fC11;
  cc[3 ]=fC20;   cc[4 ]=fC21;   cc[5 ]=fC22;
  cc[6 ]=fC30;   cc[7 ]=fC31;   cc[8 ]=fC32;   cc[9 ]=fC33;
  cc[10]=fC40*a; cc[11]=fC41*a; cc[12]=fC42*a; cc[13]=fC43*a; cc[14]=fC44*a*a;
}

//____________________________________________________________________________
Int_t AliITStrackV2::PropagateToVertex(Double_t x0,Double_t rho,Double_t pm) {
  //------------------------------------------------------------------
  //This function propagates a track to the minimal distance from the origin
  //------------------------------------------------------------------
  Double_t xv=fP2*(fX*fP2 - fP0*TMath::Sqrt(1.- fP2*fP2)); //linear approxim.
  Propagate(fAlpha,xv,0.,0.,pm);   
  return 0;
}

//____________________________________________________________________________
Int_t AliITStrackV2::
GetGlobalXYZat(Double_t xk, Double_t &x, Double_t &y, Double_t &z) const {
  //------------------------------------------------------------------
  //This function returns a track position in the global system
  //------------------------------------------------------------------
  Double_t dx=xk-fX;
  Double_t f1=fP2, f2=f1 + fP4*dx;
  if (TMath::Abs(f2) >= 0.99999) {
     Int_t n=GetNumberOfClusters();
    if (n>kWARN) cerr<<n<<" AliITStrackV2 warning: Propagation failed !\n";
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
Double_t AliITStrackV2::GetPredictedChi2(const AliCluster *c) const 
{
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  Double_t r00=c->GetSigmaY2(), r01=0., r11=c->GetSigmaZ2();
  r00+=fC00; r01+=fC10; r11+=fC11;

  Double_t det=r00*r11 - r01*r01;
  if (TMath::Abs(det) < 1.e-10) {
    Int_t n=GetNumberOfClusters();
    if (n>4) cerr<<n<<" AliKalmanTrack warning: Singular matrix !\n";
    return 1e10;
  }
  Double_t tmp=r00; r00=r11; r11=tmp; r01=-r01;
  
  Double_t dy=c->GetY() - fP0, dz=c->GetZ() - fP1;
  
  return (dy*r00*dy + 2*r01*dy*dz + dz*r11*dz)/det;
}

//_____________________________________________________________________________
Double_t AliITStrackV2::GetPredictedChi2(const AliCluster *c,Double_t *m,
Double_t x0, Double_t pm) const {
  //-----------------------------------------------------------------
  // This function calculates a chi2 increment with a vertex contraint 
  //-----------------------------------------------------------------
  TVectorD x(5); x(0)=fP0; x(1)=fP1; x(2)=fP2; x(3)=fP3; x(4)=fP4;
  TMatrixD C(5,5);
  C(0,0)=fC00; 
  C(1,0)=fC10; C(1,1)=fC11; 
  C(2,0)=fC20; C(2,1)=fC21; C(2,2)=fC22;
  C(3,0)=fC30; C(3,1)=fC31; C(3,2)=fC32; C(3,3)=fC33;
  C(4,0)=fC40; C(4,1)=fC41; C(4,2)=fC42; C(4,3)=fC43; C(4,4)=fC44;

  C(0,1)=C(1,0);
  C(0,2)=C(2,0); C(1,2)=C(2,1);
  C(0,3)=C(3,0); C(1,3)=C(3,1); C(2,3)=C(3,2);
  C(0,4)=C(4,0); C(1,4)=C(4,1); C(2,4)=C(4,2); C(3,4)=C(4,3);

  TMatrixD H(4,5); H.UnitMatrix();
  Double_t dy=(c->GetY() - m[0]), dz=(c->GetZ() - m[1]);

  Double_t dr=TMath::Sqrt(fX*fX + dy*dy);
  Double_t r =TMath::Sqrt(4/dr/dr - fP4*fP4);
  Double_t sn=0.5*(fP4*fX + dy*r);
  Double_t tg=0.5*fP4*dz/TMath::ASin(0.5*fP4*dr);
  TVectorD mm(4); 
  mm(0)=m[0]=c->GetY(); mm(1)=m[1]=c->GetZ(); mm(2)=m[2]=sn; mm(3)=m[3]=tg;

  Double_t v22=0.,v33=0.;
  //x0=0.;
  if (x0!=0.) {
     Double_t pp2=(1.+ GetTgl()*GetTgl())/(Get1Pt()*Get1Pt());
     Double_t beta2=pp2/(pp2 + pm*pm);
     x0*=TMath::Sqrt((1.+ GetTgl()*GetTgl())/(1.- GetSnp()*GetSnp()));
     Double_t theta2=14.1*14.1/(beta2*pp2*1e6)*x0;
     v22 = theta2*(1.- GetSnp()*GetSnp())*(1. + GetTgl()*GetTgl());
     v33 = theta2*(1.+ GetTgl()*GetTgl())*(1. + GetTgl()*GetTgl());
  }
  Double_t sy2=c->GetSigmaY2(), sz2=c->GetSigmaZ2();
  v22+=kSigmaYV*kSigmaYV/dr/dr;
  v22+=sy2/dr/dr;
  Double_t v20=sy2/dr;

  v33+=kSigmaZV*kSigmaZV/dr/dr;
  v33+=sz2/dr/dr;
  Double_t v31=sz2/dr;

  TMatrixD V(4,4); 
  V(0,0)=m[4 ]=sy2; V(0,1)=m[5 ]=0.;  V(0,2)=m[6 ]=v20; V(0,3)=m[7 ]=0.;
  V(1,0)=m[8 ]=0.;  V(1,1)=m[9 ]=sz2; V(1,2)=m[10]=0.;  V(1,3)=m[11]=v31;
  V(2,0)=m[12]=v20; V(2,1)=m[13]=0.;  V(2,2)=m[14]=v22; V(2,3)=m[15]=0.;
  V(3,0)=m[16]=0.;  V(3,1)=m[17]=v31; V(3,2)=m[18]=0.;  V(3,3)=m[19]=v33;

  TVectorD res=x;  res*=H; res-=mm; //res*=-1; 
  TMatrixD tmp(H,TMatrixD::kMult,C);
  TMatrixD R(tmp,TMatrixD::kMult,TMatrixD(TMatrixD::kTransposed,H)); R+=V;
  
  Double_t det=R.Determinant();
  if (TMath::Abs(det) < 1.e-25) {
    Int_t n=GetNumberOfClusters();
    if (n>kWARN) cerr<<n<<" AliITStrackV2 warning: Singular matrix !\n";
    return 1e10;
  }

  R.Invert();

  TVectorD rs=res;
  res*=R;
  return rs*res;
}

//____________________________________________________________________________
Int_t 
AliITStrackV2::PropagateTo(Double_t xk,Double_t x0,Double_t rho,Double_t pm) {
  //------------------------------------------------------------------
  //This function propagates a track
  //------------------------------------------------------------------
  Double_t x1=fX, x2=xk, dx=x2-x1;
  Double_t f1=fP2, f2=f1 + fP4*dx;
  if (TMath::Abs(f2) >= 0.99999) {
    Int_t n=GetNumberOfClusters();
    if (n>kWARN) cerr<<n<<" AliITStrackV2 warning: Propagation failed !\n";
    return 0;
  }

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

  Double_t p2=(1.+ GetTgl()*GetTgl())/(Get1Pt()*Get1Pt());
  Double_t beta2=p2/(p2 + pm*pm);

  //Multiple scattering******************
  //x0=0.;
  if (x0!=0) {
     x0*=TMath::Sqrt((1.+ fP3*fP3)/(1.- fP2*fP2));
     Double_t theta2=14.1*14.1/(beta2*p2*1e6)*x0;
     fC22 += theta2*(1.- fP2*fP2)*(1. + fP3*fP3);
     fC33 += theta2*(1. + fP3*fP3)*(1. + fP3*fP3);
     fC43 += theta2*fP3*fP4*(1. + fP3*fP3);
     fC44 += theta2*fP3*fP4*fP3*fP4;
  }

  //Energy losses************************
  if (rho!=0.) {
     rho*=TMath::Sqrt((1.+ fP3*fP3)/(1.- fP2*fP2));
     Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2)) - beta2)*rho;
     if (x1 < x2) dE=-dE;
     fP4*=(1.- sqrt(p2+pm*pm)/p2*dE);
  }

  if (!Invariant()) {cout<<"Propagate !\n"; return 0;}

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
  /*
  if (TMath::Abs(sf) >= 0.99999) {
    Int_t n=GetNumberOfClusters();
    if (n>kWARN) cerr<<n<<" AliITStrackV2 warning: Filtering failed !\n";
    return 0;
  }
  */
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

  Int_t n=GetNumberOfClusters();
  fIndex[n]=index;
  SetNumberOfClusters(n+1);
  SetChi2(GetChi2()+chi2);

  return 1;
}


//____________________________________________________________________________
Int_t AliITStrackV2::Update(const Double_t* m, Double_t chi2, UInt_t index) {
  //------------------------------------------------------------------
  //This function updates track parameters with a vertex constraint
  //------------------------------------------------------------------
  Double_t p0=fP0,p1=fP1,p2=fP2,p3=fP3,p4=fP4;
  Double_t c00=fC00;
  Double_t c10=fC10, c11=fC11;
  Double_t c20=fC20, c21=fC21, c22=fC22;
  Double_t c30=fC30, c31=fC31, c32=fC32, c33=fC33;
  Double_t c40=fC40, c41=fC41, c42=fC42, c43=fC43, c44=fC44;


  TVectorD x(5); x(0)=fP0; x(1)=fP1; x(2)=fP2; x(3)=fP3; x(4)=fP4;
  TMatrixD C(5,5);
  C(0,0)=fC00; 
  C(1,0)=fC10; C(1,1)=fC11; 
  C(2,0)=fC20; C(2,1)=fC21; C(2,2)=fC22;
  C(3,0)=fC30; C(3,1)=fC31; C(3,2)=fC32; C(3,3)=fC33;
  C(4,0)=fC40; C(4,1)=fC41; C(4,2)=fC42; C(4,3)=fC43; C(4,4)=fC44;

  C(0,1)=C(1,0);
  C(0,2)=C(2,0); C(1,2)=C(2,1);
  C(0,3)=C(3,0); C(1,3)=C(3,1); C(2,3)=C(3,2);
  C(0,4)=C(4,0); C(1,4)=C(4,1); C(2,4)=C(4,2); C(3,4)=C(4,3);

  TMatrixD H(4,5); H.UnitMatrix();
  TMatrixD Ht(TMatrixD::kTransposed,H);
  TVectorD mm(4); mm(0)=m[0]; mm(1)=m[1]; mm(2)=m[2]; mm(3)=m[3];
  TMatrixD V(4,4); 
  V(0,0)=m[4 ]; V(0,1)=m[5 ]; V(0,2)=m[6 ]; V(0,3)=m[7 ];
  V(1,0)=m[8 ]; V(1,1)=m[9 ]; V(1,2)=m[10]; V(1,3)=m[11];
  V(2,0)=m[12]; V(2,1)=m[13]; V(2,2)=m[14]; V(2,3)=m[15];
  V(3,0)=m[16]; V(3,1)=m[17]; V(3,2)=m[18]; V(3,3)=m[19];

  TMatrixD tmp(H,TMatrixD::kMult,C);
  TMatrixD R(tmp,TMatrixD::kMult,Ht); R+=V;

  R.Invert();
  
  TMatrixD K(C,TMatrixD::kMult,Ht); K*=R;
  
  TVectorD savex=x;
  x*=H; x-=mm; x*=-1; x*=K; x+=savex;

  TMatrixD saveC=C;
  C.Mult(K,tmp); C-=saveC; C*=-1;

  fP0=x(0); fP1=x(1); fP2=x(2); fP3=x(3); fP4=x(4);
  fC00=C(0,0); 
  fC10=C(1,0); fC11=C(1,1); 
  fC20=C(2,0); fC21=C(2,1); fC22=C(2,2);
  fC30=C(3,0); fC31=C(3,1); fC32=C(3,2); fC33=C(3,3);
  fC40=C(4,0); fC41=C(4,1); fC42=C(4,2); fC43=C(4,3); fC44=C(4,4);


  if (!Invariant()) {
     fP0=p0; fP1=p1; fP2=p2; fP3=p3; fP4=p4;
     fC00=c00;
     fC10=c10; fC11=c11;
     fC20=c20; fC21=c21; fC22=c22;
     fC30=c30; fC31=c31; fC32=c32; fC33=c33;
     fC40=c40; fC41=c41; fC42=c42; fC43=c43; fC44=c44;
     return 0;
  }

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
  //if (TMath::Abs(fP1)>11.5)
  //if (fP1*fP4<0) {cout<<"fP1*fP4="<<fP1*fP4<<' '<<fP1<<endl; return 0;}
  if (TMath::Abs(fP2)>=1) {cout<<"fP2="<<fP2<<endl; return 0;}

  if (fC00<=0) {cout<<"fC00="<<fC00<<endl; return 0;}
  if (fC11<=0) {cout<<"fC11="<<fC11<<endl; return 0;}
  if (fC22<=0) {cout<<"fC22="<<fC22<<endl; return 0;}
  if (fC33<=0) {cout<<"fC33="<<fC33<<endl; return 0;}
  if (fC44<=0) {cout<<"fC44="<<fC44<<endl; return 0;}

   TMatrixD m(5,5);
   m(0,0)=fC00; 
   m(1,0)=fC10; m(1,1)=fC11; 
   m(2,0)=fC20; m(2,1)=fC21; m(2,2)=fC22;
   m(3,0)=fC30; m(3,1)=fC31; m(3,2)=fC32; m(3,3)=fC33;
   m(4,0)=fC40; m(4,1)=fC41; m(4,2)=fC42; m(4,3)=fC43; m(4,4)=fC44;

   m(0,1)=m(1,0);
   m(0,2)=m(2,0); m(1,2)=m(2,1);
   m(0,3)=m(3,0); m(1,3)=m(3,1); m(2,3)=m(3,2);
   m(0,4)=m(4,0); m(1,4)=m(4,1); m(2,4)=m(4,2); m(3,4)=m(4,3);
   /*
   Double_t det=m.Determinant(); 

   if (det <= 0) {
       cout<<" bad determinant "<<det<<endl;
       m.Print(); 
       return 0;
   }
   */
   return 1;
}

//____________________________________________________________________________
Int_t AliITStrackV2::Propagate(Double_t alp, Double_t xk,
Double_t x0,Double_t rho,Double_t pm) {
  //------------------------------------------------------------------
  //This function propagates a track
  //------------------------------------------------------------------
  Double_t p0=fP0,p1=fP1,p2=fP2,p3=fP3,p4=fP4;
  Double_t c00=fC00;
  Double_t c10=fC10, c11=fC11;
  Double_t c20=fC20, c21=fC21, c22=fC22;
  Double_t c30=fC30, c31=fC31, c32=fC32, c33=fC33;
  Double_t c40=fC40, c41=fC41, c42=fC42, c43=fC43, c44=fC44;


  Double_t dalp=alp-fAlpha;

  Double_t ca=TMath::Cos(dalp), sa=TMath::Sin(dalp);
  Double_t sf=fP2, cf=TMath::Sqrt(1.- fP2*fP2);  

  Double_t pp2=fP2*ca - cf*sa;
  if (TMath::Abs(pp2) >= 0.99999) {
     Int_t n=GetNumberOfClusters();
     if (n>kWARN) cerr<<n<<" AliITStrackV2 warning: Rotation failed !\n";
     return 0;
  }

  fAlpha = alp;
  if      (fAlpha < -TMath::Pi()) fAlpha += 2*TMath::Pi();
  else if (fAlpha >= TMath::Pi()) fAlpha -= 2*TMath::Pi();
  
  Double_t x1=fX, y1=fP0;

  fX = x1*ca + y1*sa;
  fP0=-x1*sa + y1*ca;
  fP2 = pp2;

  cf=ca + sf*sa/cf;

  if (!Invariant()) {cout<<dalp<<" Rotate !\n"; return 0;}

  x1=fX; Double_t x2=xk, dx=x2-x1;
  Double_t f1=fP2, f2=f1 + fP4*dx;
  if (TMath::Abs(f2) >= 0.99999) {
    Int_t n=GetNumberOfClusters();
    if (n>kWARN) cerr<<n<<" AliITStrackV2 warning: Propagation failed !\n";
    return 0;
  }

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
  /*
  //b = C*ft
  Double_t b00=f02*fC20 + f03*fC30, b01=f12*fC20 + f13*fC30 + f14*fC40;
  Double_t b02=f23*fC30;
  Double_t b10=f02*fC21 + f03*fC31, b11=f12*fC21 + f13*fC31 + f14*fC41;
  Double_t b12=f23*fC31;
  Double_t b20=f02*fC22 + f03*fC32, b21=f12*fC22 + f13*fC32 + f14*fC42;
  Double_t b22=f23*fC32;
  Double_t b30=f02*fC32 + f03*fC33, b31=f12*fC32 + f13*fC33 + f14*fC43;
  Double_t b32=f23*fC33;
  Double_t b40=f02*fC42 + f03*fC43, b41=f12*fC42 + f13*fC43 + f14*fC44;
  Double_t b42=f23*fC43;
  
  //a = f*b = f*C*ft
  Double_t a00=f02*b20+f03*b30,a01=f02*b21+f03*b31,a02=f02*b22+f03*b32;
  Double_t a11=f12*b21+f13*b31+f14*b41,a12=f12*b22+f13*b32+f14*b42;
  Double_t a22=f23*b32;

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
*/

 TMatrixD F(5,5); F.UnitMatrix();
 F(0,0)=-(f1+f2)/(r1+r2)*sa + ca; F(0,2)=f02*cf; F(0,4)=f04;
 F(1,0)=-(f1+f2)/(f1*r2 + f2*r1)*fP3*sa; F(1,2)=f12*cf; F(1,4)=f14; F(1,3)=f13;
 F(2,0)=-fP4*sa; F(2,2)=cf; F(2,4)=f24;

  TMatrixD C(5,5);
  C(0,0)=fC00; 
  C(1,0)=fC10; C(1,1)=fC11; 
  C(2,0)=fC20; C(2,1)=fC21; C(2,2)=fC22;
  C(3,0)=fC30; C(3,1)=fC31; C(3,2)=fC32; C(3,3)=fC33;
  C(4,0)=fC40; C(4,1)=fC41; C(4,2)=fC42; C(4,3)=fC43; C(4,4)=fC44;

  C(0,1)=C(1,0);
  C(0,2)=C(2,0); C(1,2)=C(2,1);
  C(0,3)=C(3,0); C(1,3)=C(3,1); C(2,3)=C(3,2);
  C(0,4)=C(4,0); C(1,4)=C(4,1); C(2,4)=C(4,2); C(3,4)=C(4,3);

  TMatrixD tmp(C,TMatrixD::kMult,TMatrixD(TMatrixD::kTransposed, F));
  C.Mult(F,tmp);

  fC00=C(0,0); 
  fC10=C(1,0); fC11=C(1,1); 
  fC20=C(2,0); fC21=C(2,1); fC22=C(2,2);
  fC30=C(3,0); fC31=C(3,1); fC32=C(3,2); fC33=C(3,3);
  fC40=C(4,0); fC41=C(4,1); fC42=C(4,2); fC43=C(4,3); fC44=C(4,4);

  pp2=(1.+ GetTgl()*GetTgl())/(Get1Pt()*Get1Pt());
  Double_t beta2=pp2/(pp2 + pm*pm);

  //Multiple scattering******************
  //x0=0.;
  if (x0!=0.) {
     x0*=TMath::Sqrt((1.+ fP3*fP3)/(1.- fP2*fP2));
     Double_t theta2=14.1*14.1/(beta2*pp2*1e6)*x0;
     fC22 += theta2*(1.- fP2*fP2)*(1. + fP3*fP3);
     fC33 += theta2*(1. + fP3*fP3)*(1. + fP3*fP3);
     fC43 += theta2*fP3*fP4*(1. + fP3*fP3);
     fC44 += theta2*fP3*fP4*fP3*fP4;
  }

  //Energy losses************************
  if (rho!=0.) {  
     rho*=TMath::Sqrt((1.+ fP3*fP3)/(1.- fP2*fP2));
     Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2)) - beta2)*rho;
     if (x1 < x2) dE=-dE;
     fP4*=(1.- sqrt(pp2+pm*pm)/pp2*dE);
  }

  if (!Invariant()) {
     fP0=p0; fP1=p1; fP2=p2; fP3=p3; fP4=p4;
     fC00=c00;
     fC10=c10; fC11=c11;
     fC20=c20; fC21=c21; fC22=c22;
     fC30=c30; fC31=c31; fC32=c32; fC33=c33;
     fC40=c40; fC41=c41; fC42=c42; fC43=c43; fC44=c44;
     return 0;
  }

  fX=x2;

  return 1;
}

Double_t AliITStrackV2::GetD() const {
  //------------------------------------------------------------------
  //This function calculates the transverse impact parameter
  //------------------------------------------------------------------
  Double_t sn=fP4*fX - fP2, cs=fP4*fP0 + TMath::Sqrt(1.- fP2*fP2);
  Double_t a=2*(fX*fP2 - fP0*TMath::Sqrt(1.- fP2*fP2))-fP4*(fX*fX + fP0*fP0);
  if (fP4<0) a=-a;
  return a/(1 + TMath::Sqrt(sn*sn + cs*cs));
} 


Int_t AliITStrackV2::Improve(Double_t x0,Double_t yv,Double_t zv) {
  //------------------------------------------------------------------
  //This function improves angular track parameters  
  //------------------------------------------------------------------
  Double_t dy=fP0-yv, dz=fP1-zv;
  Double_t r2=fX*fX+dy*dy;
  Double_t p2=(1.+ GetTgl()*GetTgl())/(Get1Pt()*Get1Pt());
  Double_t beta2=p2/(p2 + 0.14*0.14);
  x0*=TMath::Sqrt((1.+ GetTgl()*GetTgl())/(1.- GetSnp()*GetSnp()));
  Double_t theta2=14.1*14.1/(beta2*p2*1e6)*x0;

  Double_t par=0.5*(fP4*fX + dy*TMath::Sqrt(4/r2-fP4*fP4));
  Double_t sigma2 = theta2*(1.- GetSnp()*GetSnp())*(1. + GetTgl()*GetTgl());
  sigma2 += fC00/r2*(1.- dy*dy/r2)*(1.- dy*dy/r2);
  sigma2 += kSigmaYV*kSigmaYV/r2;
  sigma2 += 0.25*fC44*fX*fX;
  Double_t eps2=sigma2/(fC22+sigma2), eps=TMath::Sqrt(eps2);
  if (10*r2*fC44<fC22) {
     fP2 = eps2*fP2 + (1-eps2)*par;
     fC22*=eps2; fC21*=eps; fC20*=eps; fC32*=eps; fC42*=eps;
  }

  par=0.5*fP4*dz/TMath::ASin(0.5*fP4*TMath::Sqrt(r2));
  sigma2=theta2;
  sigma2 += fC11/r2+fC00*dy*dy*dz*dz/(r2*r2*r2);
  sigma2 += kSigmaZV*kSigmaZV/r2;
  eps2=sigma2/(fC33+sigma2); eps=TMath::Sqrt(eps2);
  Double_t tgl=fP3;
  fP3 = eps2*fP3 + (1-eps2)*par;
  fC33*=eps2; fC32*=eps; fC31*=eps; fC30*=eps; fC43*=eps;

  eps=TMath::Sqrt((1+fP3*fP3)/(1+tgl*tgl));
  fP4*=eps;
  fC44*=eps*eps; fC43*=eps;fC42*=eps; fC41*=eps; fC40*=eps;

  if (!Invariant()) return 0;
  return 1;
} 

/*
Int_t AliITStrackV2::Improve(Double_t x0,Double_t xv,Double_t yv) {
  //------------------------------------------------------------------
  //This function improves angular track parameters  
  //------------------------------------------------------------------
  TMatrixD I(5,5);
  TVectorD v(5); v(0)=fP0; v(1)=fP1; v(2)=fP2; v(3)=fP3; v(4)=fP4;

  Double_t r2=fX*fX+fP0*fP0;
  Double_t p2=(1.+ GetTgl()*GetTgl())/(Get1Pt()*Get1Pt());
  Double_t beta2=p2/(p2 + 0.14*0.14);
  x0*=TMath::Sqrt((1.+ GetTgl()*GetTgl())/(1.- GetSnp()*GetSnp()));
  Double_t theta2=14.1*14.1/(beta2*p2*1e6)*x0;

  v(2)=0.5*(fP4*fX + fP0*TMath::Sqrt(4/r2-fP4*fP4));
  Double_t sigma2 = theta2*(1.- GetSnp()*GetSnp())*(1. + GetTgl()*GetTgl());
  sigma2 += fC00/r2*(1.- fP0*fP0/r2)*(1.- fP0*fP0/r2);
  sigma2 += kSigmaYV*kSigmaYV/r2;
  I(2,2)=1/sigma2;

  v(3)=0.5*fP4*fP1/TMath::ASin(0.5*fP4*TMath::Sqrt(r2));
  sigma2=theta2;
  sigma2 += fC11/r2+fC00*fP0*fP0*fP1*fP1/(r2*r2*r2);
  sigma2 += kSigmaZV*kSigmaZV/r2;
  I(3,3)=1/sigma2;

  Double_t tgl=fP3;

  TVectorD x(5); x(0)=fP0; x(1)=fP1; x(2)=fP2; x(3)=fP3; x(4)=fP4;
  TMatrixD C(5,5);
  C(0,0)=fC00; 
  C(1,0)=fC10; C(1,1)=fC11; 
  C(2,0)=fC20; C(2,1)=fC21; C(2,2)=fC22;
  C(3,0)=fC30; C(3,1)=fC31; C(3,2)=fC32; C(3,3)=fC33;
  C(4,0)=fC40; C(4,1)=fC41; C(4,2)=fC42; C(4,3)=fC43; C(4,4)=fC44;

  C(0,1)=C(1,0);
  C(0,2)=C(2,0); C(1,2)=C(2,1);
  C(0,3)=C(3,0); C(1,3)=C(3,1); C(2,3)=C(3,2);
  C(0,4)=C(4,0); C(1,4)=C(4,1); C(2,4)=C(4,2); C(3,4)=C(4,3);

  TMatrixD tmp(I,TMatrixD::kMult,C),U(5,5); U.UnitMatrix();
  U+=tmp;
  U.Invert();
  TMatrixD W1(U);
  TMatrixD W2(tmp,TMatrixD::kMult,W1);

  v*=W2; x*=W1; x+=v;

  C*=W1;


  fP0=x(0); fP1=x(1); fP2=x(2); fP3=x(3); fP4=x(4);
  fC00=C(0,0); 
  fC10=C(1,0); fC11=C(1,1); 
  fC20=C(2,0); fC21=C(2,1); fC22=C(2,2);
  fC30=C(3,0); fC31=C(3,1); fC32=C(3,2); fC33=C(3,3);
  fC40=C(4,0); fC41=C(4,1); fC42=C(4,2); fC43=C(4,3); fC44=C(4,4);

  eps=TMath::Sqrt((1+fP3*fP3)/(1+tgl*tgl));
  fP4*=eps;
  fC44*=eps*eps; fC43*=eps;fC42*=eps; fC41*=eps; fC40*=eps;

  if (!Invariant()) return 0;
  return 1;
} 
*/

