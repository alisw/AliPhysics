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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Implementation of the external track parameterisation class.              //
//                                                                           //
// This parameterisation is used to exchange tracks between the detectors.   //
// A set of functions returning the position and the momentum of tracks      //
// in the global coordinate system as well as the track impact parameters    //
// are implemented.
// Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                            //
///////////////////////////////////////////////////////////////////////////////
#include "AliExternalTrackParam.h"
#include "AliESDVertex.h"
#include "AliLog.h"

ClassImp(AliExternalTrackParam)

//_____________________________________________________________________________
AliExternalTrackParam::AliExternalTrackParam() :
  TObject(),
  fX(0),
  fAlpha(0)
{
  //
  // default constructor
  //
  for (Int_t i = 0; i < 5; i++) fP[i] = 0;
  for (Int_t i = 0; i < 15; i++) fC[i] = 0;
}

//_____________________________________________________________________________
AliExternalTrackParam::AliExternalTrackParam(const AliExternalTrackParam &track):
  TObject(track),
  fX(track.fX),
  fAlpha(track.fAlpha)
{
  //
  // copy constructor
  //
  for (Int_t i = 0; i < 5; i++) fP[i] = track.fP[i];
  for (Int_t i = 0; i < 15; i++) fC[i] = track.fC[i];
}

//_____________________________________________________________________________
AliExternalTrackParam::AliExternalTrackParam(Double_t x, Double_t alpha, 
					     const Double_t param[5], 
					     const Double_t covar[15]) :
  TObject(),
  fX(x),
  fAlpha(alpha)
{
  //
  // create external track parameters from given arguments
  //
  for (Int_t i = 0; i < 5; i++)  fP[i] = param[i];
  for (Int_t i = 0; i < 15; i++) fC[i] = covar[i];
}

//_____________________________________________________________________________
void AliExternalTrackParam::Set(Double_t x, Double_t alpha,
				const Double_t p[5], const Double_t cov[15]) {
  //
  //  Sets the parameters
  //
  fX=x;
  fAlpha=alpha;
  for (Int_t i = 0; i < 5; i++)  fP[i] = p[i];
  for (Int_t i = 0; i < 15; i++) fC[i] = cov[i];
}

//_____________________________________________________________________________
void AliExternalTrackParam::Reset() {
  //
  // Resets all the parameters to 0 
  //
  fX=fAlpha=0.;
  for (Int_t i = 0; i < 5; i++) fP[i] = 0;
  for (Int_t i = 0; i < 15; i++) fC[i] = 0;
}

Double_t AliExternalTrackParam::GetP() const {
  //---------------------------------------------------------------------
  // This function returns the track momentum
  // Results for (nearly) straight tracks are meaningless !
  //---------------------------------------------------------------------
  if (TMath::Abs(fP[4])<=kAlmost0) return kVeryBig;
  return TMath::Sqrt(1.+ fP[3]*fP[3])/TMath::Abs(fP[4]);
}

Double_t AliExternalTrackParam::Get1P() const {
  //---------------------------------------------------------------------
  // This function returns the 1/(track momentum)
  //---------------------------------------------------------------------
  return TMath::Abs(fP[4])/TMath::Sqrt(1.+ fP[3]*fP[3]);
}

//_______________________________________________________________________
Double_t AliExternalTrackParam::GetD(Double_t x,Double_t y,Double_t b) const {
  //------------------------------------------------------------------
  // This function calculates the transverse impact parameter
  // with respect to a point with global coordinates (x,y)
  // in the magnetic field "b" (kG)
  //------------------------------------------------------------------
  if (TMath::Abs(b) < kAlmost0Field) return GetLinearD(x,y);
  Double_t rp4=GetC(b);

  Double_t xt=fX, yt=fP[0];

  Double_t sn=TMath::Sin(fAlpha), cs=TMath::Cos(fAlpha);
  Double_t a = x*cs + y*sn;
  y = -x*sn + y*cs; x=a;
  xt-=x; yt-=y;

  sn=rp4*xt - fP[2]; cs=rp4*yt + TMath::Sqrt(1.- fP[2]*fP[2]);
  a=2*(xt*fP[2] - yt*TMath::Sqrt(1.- fP[2]*fP[2]))-rp4*(xt*xt + yt*yt);
  return  -a/(1 + TMath::Sqrt(sn*sn + cs*cs));
}

//_______________________________________________________________________
void AliExternalTrackParam::
GetDZ(Double_t x, Double_t y, Double_t z, Double_t b, Float_t dz[2]) const {
  //------------------------------------------------------------------
  // This function calculates the transverse and longitudinal impact parameters
  // with respect to a point with global coordinates (x,y)
  // in the magnetic field "b" (kG)
  //------------------------------------------------------------------
  Double_t f1 = fP[2], r1 = TMath::Sqrt(1. - f1*f1);
  Double_t xt=fX, yt=fP[0];
  Double_t sn=TMath::Sin(fAlpha), cs=TMath::Cos(fAlpha);
  Double_t a = x*cs + y*sn;
  y = -x*sn + y*cs; x=a;
  xt-=x; yt-=y;

  Double_t rp4=GetC(b);
  if ((TMath::Abs(b) < kAlmost0Field) || (TMath::Abs(rp4) < kAlmost0)) {
     dz[0] = -(xt*f1 - yt*r1);
     dz[1] = fP[1] + (dz[0]*f1 - xt)/r1*fP[3] - z;
     return;
  }

  sn=rp4*xt - f1; cs=rp4*yt + r1;
  a=2*(xt*f1 - yt*r1)-rp4*(xt*xt + yt*yt);
  Double_t rr=TMath::Sqrt(sn*sn + cs*cs);
  dz[0] = -a/(1 + rr);
  Double_t f2 = -sn/rr, r2 = TMath::Sqrt(1. - f2*f2);
  dz[1] = fP[1] + fP[3]/rp4*TMath::ASin(f2*r1 - f1*r2) - z;
}

//_______________________________________________________________________
Double_t AliExternalTrackParam::GetLinearD(Double_t xv,Double_t yv) const {
  //------------------------------------------------------------------
  // This function calculates the transverse impact parameter
  // with respect to a point with global coordinates (xv,yv)
  // neglecting the track curvature.
  //------------------------------------------------------------------
  Double_t sn=TMath::Sin(fAlpha), cs=TMath::Cos(fAlpha);
  Double_t x= xv*cs + yv*sn;
  Double_t y=-xv*sn + yv*cs;

  Double_t d = (fX-x)*fP[2] - (fP[0]-y)*TMath::Sqrt(1.- fP[2]*fP[2]);

  return -d;
}

Bool_t AliExternalTrackParam::
CorrectForMaterial(Double_t d,  Double_t x0, Double_t mass) {
  //------------------------------------------------------------------
  // This function corrects the track parameters for the crossed material
  // "d"    - the thickness (fraction of the radiation length)
  // "x0"   - the radiation length (g/cm^2) 
  // "mass" - the mass of this particle (GeV/c^2)
  //------------------------------------------------------------------
  Double_t &fP2=fP[2];
  Double_t &fP3=fP[3];
  Double_t &fP4=fP[4];

  Double_t &fC22=fC[5];
  Double_t &fC33=fC[9];
  Double_t &fC43=fC[13];
  Double_t &fC44=fC[14];

  Double_t p=GetP();
  Double_t p2=p*p;
  Double_t beta2=p2/(p2 + mass*mass);
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
  if (x0!=0. && beta2<1) {
     d*=x0;
     Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2)) - beta2)*d;
     if (beta2/(1-beta2)>3.5*3.5)
       dE=0.153e-3/beta2*(log(3.5*5940)+0.5*log(beta2/(1-beta2)) - beta2)*d;

     fP4*=(1.- TMath::Sqrt(p2 + mass*mass)/p2*dE);
  }

  return kTRUE;
}

Bool_t AliExternalTrackParam::Rotate(Double_t alpha) {
  //------------------------------------------------------------------
  // Transform this track to the local coord. system rotated
  // by angle "alpha" (rad) with respect to the global coord. system. 
  //------------------------------------------------------------------
  if (TMath::Abs(fP[2]) >= kAlmost1) {
     AliError(Form("Precondition is not satisfied: |sin(phi)|>1 ! %f",fP[2])); 
     return kFALSE;
  }
 
  if      (alpha < -TMath::Pi()) alpha += 2*TMath::Pi();
  else if (alpha >= TMath::Pi()) alpha -= 2*TMath::Pi();

  Double_t &fP0=fP[0];
  Double_t &fP2=fP[2];
  Double_t &fC00=fC[0];
  Double_t &fC10=fC[1];
  Double_t &fC20=fC[3];
  Double_t &fC21=fC[4];
  Double_t &fC22=fC[5];
  Double_t &fC30=fC[6];
  Double_t &fC32=fC[8];
  Double_t &fC40=fC[10];
  Double_t &fC42=fC[12];

  Double_t x=fX;
  Double_t ca=TMath::Cos(alpha-fAlpha), sa=TMath::Sin(alpha-fAlpha);
  Double_t sf=fP2, cf=TMath::Sqrt(1.- fP2*fP2);

  Double_t tmp=sf*ca - cf*sa;
  if (TMath::Abs(tmp) >= kAlmost1) return kFALSE;

  fAlpha = alpha;
  fX =  x*ca + fP0*sa;
  fP0= -x*sa + fP0*ca;
  fP2=  tmp;

  if (TMath::Abs(cf)<kAlmost0) {
    AliError(Form("Too small cosine value %f",cf)); 
    cf = kAlmost0;
  } 

  Double_t rr=(ca+sf/cf*sa);  

  fC00 *= (ca*ca);
  fC10 *= ca;
  fC20 *= ca*rr;
  fC21 *= rr;
  fC22 *= rr*rr;
  fC30 *= ca;
  fC32 *= rr;
  fC40 *= ca;
  fC42 *= rr;

  return kTRUE;
}

Bool_t AliExternalTrackParam::PropagateTo(Double_t xk, Double_t b) {
  //----------------------------------------------------------------
  // Propagate this track to the plane X=xk (cm) in the field "b" (kG)
  //----------------------------------------------------------------
  Double_t dx=xk-fX;
  if (TMath::Abs(dx)<=kAlmost0)  return kTRUE;

  Double_t crv=GetC(b);
  if (TMath::Abs(b) < kAlmost0Field) crv=0.;

  Double_t f1=fP[2], f2=f1 + crv*dx;
  if (TMath::Abs(f1) >= kAlmost1) return kFALSE;
  if (TMath::Abs(f2) >= kAlmost1) return kFALSE;

  Double_t &fP0=fP[0], &fP1=fP[1], &fP2=fP[2], &fP3=fP[3], &fP4=fP[4];
  Double_t 
  &fC00=fC[0],
  &fC10=fC[1],   &fC11=fC[2],  
  &fC20=fC[3],   &fC21=fC[4],   &fC22=fC[5],
  &fC30=fC[6],   &fC31=fC[7],   &fC32=fC[8],   &fC33=fC[9],  
  &fC40=fC[10],  &fC41=fC[11],  &fC42=fC[12],  &fC43=fC[13], &fC44=fC[14];

  Double_t r1=TMath::Sqrt(1.- f1*f1), r2=TMath::Sqrt(1.- f2*f2);

  fX=xk;
  fP0 += dx*(f1+f2)/(r1+r2);
  fP1 += dx*(r2 + f2*(f1+f2)/(r1+r2))*fP3;  // Many thanks to P.Hristov !
  fP2 += dx*crv;

  //f = F - 1
   
  Double_t f02=    dx/(r1*r1*r1);            Double_t cc=crv/fP4;
  Double_t f04=0.5*dx*dx/(r1*r1*r1);         f04*=cc;
  Double_t f12=    dx*fP3*f1/(r1*r1*r1);
  Double_t f14=0.5*dx*dx*fP3*f1/(r1*r1*r1);  f14*=cc;
  Double_t f13=    dx/r1;
  Double_t f24=    dx;                       f24*=cc;
  
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

  return kTRUE;
}

Double_t 
AliExternalTrackParam::GetPredictedChi2(Double_t p[2],Double_t cov[3]) const {
  //----------------------------------------------------------------
  // Estimate the chi2 of the space point "p" with the cov. matrix "cov"
  //----------------------------------------------------------------
  Double_t sdd = fC[0] + cov[0]; 
  Double_t sdz = fC[1] + cov[1];
  Double_t szz = fC[2] + cov[2];
  Double_t det = sdd*szz - sdz*sdz;

  if (TMath::Abs(det) < kAlmost0) return kVeryBig;

  Double_t d = fP[0] - p[0];
  Double_t z = fP[1] - p[1];

  return (d*szz*d - 2*d*sdz*z + z*sdd*z)/det;
}

Bool_t AliExternalTrackParam::Update(Double_t p[2], Double_t cov[3]) {
  //------------------------------------------------------------------
  // Update the track parameters with the space point "p" having
  // the covariance matrix "cov"
  //------------------------------------------------------------------
  Double_t &fP0=fP[0], &fP1=fP[1], &fP2=fP[2], &fP3=fP[3], &fP4=fP[4];
  Double_t 
  &fC00=fC[0],
  &fC10=fC[1],   &fC11=fC[2],  
  &fC20=fC[3],   &fC21=fC[4],   &fC22=fC[5],
  &fC30=fC[6],   &fC31=fC[7],   &fC32=fC[8],   &fC33=fC[9],  
  &fC40=fC[10],  &fC41=fC[11],  &fC42=fC[12],  &fC43=fC[13], &fC44=fC[14];

  Double_t r00=cov[0], r01=cov[1], r11=cov[2];
  r00+=fC00; r01+=fC10; r11+=fC11;
  Double_t det=r00*r11 - r01*r01;

  if (TMath::Abs(det) < kAlmost0) return kFALSE;


  Double_t tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;
 
  Double_t k00=fC00*r00+fC10*r01, k01=fC00*r01+fC10*r11;
  Double_t k10=fC10*r00+fC11*r01, k11=fC10*r01+fC11*r11;
  Double_t k20=fC20*r00+fC21*r01, k21=fC20*r01+fC21*r11;
  Double_t k30=fC30*r00+fC31*r01, k31=fC30*r01+fC31*r11;
  Double_t k40=fC40*r00+fC41*r01, k41=fC40*r01+fC41*r11;

  Double_t dy=p[0] - fP0, dz=p[1] - fP1;
  Double_t sf=fP2 + k20*dy + k21*dz;
  if (TMath::Abs(sf) > kAlmost1) return kFALSE;  
  
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

  return kTRUE;
}

void 
AliExternalTrackParam::GetHelixParameters(Double_t hlx[6], Double_t b) const {
  //--------------------------------------------------------------------
  // External track parameters -> helix parameters 
  // "b" - magnetic field (kG)
  //--------------------------------------------------------------------
  Double_t cs=TMath::Cos(fAlpha), sn=TMath::Sin(fAlpha);
  
  hlx[0]=fP[0]; hlx[1]=fP[1]; hlx[2]=fP[2]; hlx[3]=fP[3];

  hlx[5]=fX*cs - hlx[0]*sn;               // x0
  hlx[0]=fX*sn + hlx[0]*cs;               // y0
//hlx[1]=                                 // z0
  hlx[2]=TMath::ASin(hlx[2]) + fAlpha;    // phi0
//hlx[3]=                                 // tgl
  hlx[4]=GetC(b);                         // C
}


static void Evaluate(const Double_t *h, Double_t t,
                     Double_t r[3],  //radius vector
                     Double_t g[3],  //first defivatives
                     Double_t gg[3]) //second derivatives
{
  //--------------------------------------------------------------------
  // Calculate position of a point on a track and some derivatives
  //--------------------------------------------------------------------
  Double_t phase=h[4]*t+h[2];
  Double_t sn=TMath::Sin(phase), cs=TMath::Cos(phase);

  r[0] = h[5] + (sn - h[6])/h[4];
  r[1] = h[0] - (cs - h[7])/h[4];  
  r[2] = h[1] + h[3]*t;

  g[0] = cs; g[1]=sn; g[2]=h[3];
  
  gg[0]=-h[4]*sn; gg[1]=h[4]*cs; gg[2]=0.;
}

Double_t AliExternalTrackParam::GetDCA(const AliExternalTrackParam *p, 
Double_t b, Double_t &xthis, Double_t &xp) const {
  //------------------------------------------------------------
  // Returns the (weighed !) distance of closest approach between 
  // this track and the track "p".
  // Other returned values:
  //   xthis, xt - coordinates of tracks' reference planes at the DCA 
  //-----------------------------------------------------------
  Double_t dy2=GetSigmaY2() + p->GetSigmaY2();
  Double_t dz2=GetSigmaZ2() + p->GetSigmaZ2();
  Double_t dx2=dy2; 

  //dx2=dy2=dz2=1.;

  Double_t p1[8]; GetHelixParameters(p1,b);
  p1[6]=TMath::Sin(p1[2]); p1[7]=TMath::Cos(p1[2]);
  Double_t p2[8]; p->GetHelixParameters(p2,b);
  p2[6]=TMath::Sin(p2[2]); p2[7]=TMath::Cos(p2[2]);


  Double_t r1[3],g1[3],gg1[3]; Double_t t1=0.;
  Evaluate(p1,t1,r1,g1,gg1);
  Double_t r2[3],g2[3],gg2[3]; Double_t t2=0.;
  Evaluate(p2,t2,r2,g2,gg2);

  Double_t dx=r2[0]-r1[0], dy=r2[1]-r1[1], dz=r2[2]-r1[2];
  Double_t dm=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;

  Int_t max=27;
  while (max--) {
     Double_t gt1=-(dx*g1[0]/dx2 + dy*g1[1]/dy2 + dz*g1[2]/dz2);
     Double_t gt2=+(dx*g2[0]/dx2 + dy*g2[1]/dy2 + dz*g2[2]/dz2);
     Double_t h11=(g1[0]*g1[0] - dx*gg1[0])/dx2 + 
                  (g1[1]*g1[1] - dy*gg1[1])/dy2 +
                  (g1[2]*g1[2] - dz*gg1[2])/dz2;
     Double_t h22=(g2[0]*g2[0] + dx*gg2[0])/dx2 + 
                  (g2[1]*g2[1] + dy*gg2[1])/dy2 +
                  (g2[2]*g2[2] + dz*gg2[2])/dz2;
     Double_t h12=-(g1[0]*g2[0]/dx2 + g1[1]*g2[1]/dy2 + g1[2]*g2[2]/dz2);

     Double_t det=h11*h22-h12*h12;

     Double_t dt1,dt2;
     if (TMath::Abs(det)<1.e-33) {
        //(quasi)singular Hessian
        dt1=-gt1; dt2=-gt2;
     } else {
        dt1=-(gt1*h22 - gt2*h12)/det; 
        dt2=-(h11*gt2 - h12*gt1)/det;
     }

     if ((dt1*gt1+dt2*gt2)>0) {dt1=-dt1; dt2=-dt2;}

     //check delta(phase1) ?
     //check delta(phase2) ?

     if (TMath::Abs(dt1)/(TMath::Abs(t1)+1.e-3) < 1.e-4)
     if (TMath::Abs(dt2)/(TMath::Abs(t2)+1.e-3) < 1.e-4) {
        if ((gt1*gt1+gt2*gt2) > 1.e-4/dy2/dy2) 
	  AliWarning(" stopped at not a stationary point !");
        Double_t lmb=h11+h22; lmb=lmb-TMath::Sqrt(lmb*lmb-4*det);
        if (lmb < 0.) 
	  AliWarning(" stopped at not a minimum !");
        break;
     }

     Double_t dd=dm;
     for (Int_t div=1 ; ; div*=2) {
        Evaluate(p1,t1+dt1,r1,g1,gg1);
        Evaluate(p2,t2+dt2,r2,g2,gg2);
        dx=r2[0]-r1[0]; dy=r2[1]-r1[1]; dz=r2[2]-r1[2];
        dd=dx*dx/dx2 + dy*dy/dy2 + dz*dz/dz2;
	if (dd<dm) break;
        dt1*=0.5; dt2*=0.5;
        if (div>512) {
           AliWarning(" overshoot !"); break;
        }   
     }
     dm=dd;

     t1+=dt1;
     t2+=dt2;

  }

  if (max<=0) AliWarning(" too many iterations !");

  Double_t cs=TMath::Cos(GetAlpha());
  Double_t sn=TMath::Sin(GetAlpha());
  xthis=r1[0]*cs + r1[1]*sn;

  cs=TMath::Cos(p->GetAlpha());
  sn=TMath::Sin(p->GetAlpha());
  xp=r2[0]*cs + r2[1]*sn;

  return TMath::Sqrt(dm*TMath::Sqrt(dy2*dz2));
}
 
Double_t AliExternalTrackParam::
PropagateToDCA(AliExternalTrackParam *p, Double_t b) {
  //--------------------------------------------------------------
  // Propagates this track and the argument track to the position of the
  // distance of closest approach.
  // Returns the (weighed !) distance of closest approach.
  //--------------------------------------------------------------
  Double_t xthis,xp;
  Double_t dca=GetDCA(p,b,xthis,xp);

  if (!PropagateTo(xthis,b)) {
    //AliWarning(" propagation failed !");
    return 1e+33;
  }

  if (!p->PropagateTo(xp,b)) {
    //AliWarning(" propagation failed !";
    return 1e+33;
  }

  return dca;
}




Bool_t AliExternalTrackParam::PropagateToDCA(const AliESDVertex *vtx, Double_t b, Double_t maxd){
  //
  // Try to relate this track to the vertex "vtx", 
  // if the (rough) transverse impact parameter is not bigger then "maxd". 
  //            Magnetic field is "b" (kG).
  //
  // a) The track gets extapolated to the DCA to the vertex.
  // b) The impact parameters and their covariance matrix are calculated.
  //
  //    In the case of success, the returned value is kTRUE
  //    (otherwise, it's kFALSE)
  //  
  Double_t alpha=GetAlpha();
  Double_t sn=TMath::Sin(alpha), cs=TMath::Cos(alpha);
  Double_t x=GetX(), y=GetParameter()[0], snp=GetParameter()[2];
  Double_t xv= vtx->GetXv()*cs + vtx->GetYv()*sn;
  Double_t yv=-vtx->GetXv()*sn + vtx->GetYv()*cs;
  x-=xv; y-=yv;

  //Estimate the impact parameter neglecting the track curvature
  Double_t d=TMath::Abs(x*snp - y*TMath::Sqrt(1.- snp*snp));
  if (d > maxd) return kFALSE; 

  //Propagate to the DCA
  Double_t crv=0.299792458e-3*b*GetParameter()[4];
  Double_t tgfv=-(crv*x - snp)/(crv*y + TMath::Sqrt(1.-snp*snp));
  sn=tgfv/TMath::Sqrt(1.+ tgfv*tgfv); cs=TMath::Sqrt(1.- sn*sn);

  x = xv*cs + yv*sn;
  yv=-xv*sn + yv*cs; xv=x;

  if (!Propagate(alpha+TMath::ASin(sn),xv,b)) return kFALSE;
  return kTRUE;
}




Bool_t Local2GlobalMomentum(Double_t p[3],Double_t alpha) {
  //----------------------------------------------------------------
  // This function performs local->global transformation of the
  // track momentum.
  // When called, the arguments are:
  //    p[0] = 1/pt of the track;
  //    p[1] = sine of local azim. angle of the track momentum;
  //    p[2] = tangent of the track momentum dip angle;
  //   alpha - rotation angle. 
  // The result is returned as:
  //    p[0] = px
  //    p[1] = py
  //    p[2] = pz
  // Results for (nearly) straight tracks are meaningless !
  //----------------------------------------------------------------
  if (TMath::Abs(p[0])<=kAlmost0) return kFALSE;
  if (TMath::Abs(p[1])> kAlmost1) return kFALSE;

  Double_t pt=1./TMath::Abs(p[0]);
  Double_t cs=TMath::Cos(alpha), sn=TMath::Sin(alpha);
  Double_t r=TMath::Sqrt(1 - p[1]*p[1]);
  p[0]=pt*(r*cs - p[1]*sn); p[1]=pt*(p[1]*cs + r*sn); p[2]=pt*p[2];

  return kTRUE;
}

Bool_t Local2GlobalPosition(Double_t r[3],Double_t alpha) {
  //----------------------------------------------------------------
  // This function performs local->global transformation of the
  // track position.
  // When called, the arguments are:
  //    r[0] = local x
  //    r[1] = local y
  //    r[2] = local z
  //   alpha - rotation angle. 
  // The result is returned as:
  //    r[0] = global x
  //    r[1] = global y
  //    r[2] = global z
  //----------------------------------------------------------------
  Double_t cs=TMath::Cos(alpha), sn=TMath::Sin(alpha), x=r[0];
  r[0]=x*cs - r[1]*sn; r[1]=x*sn + r[1]*cs;

  return kTRUE;
}

Bool_t AliExternalTrackParam::GetPxPyPz(Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the global track momentum components
  // Results for (nearly) straight tracks are meaningless !
  //---------------------------------------------------------------------
  p[0]=fP[4]; p[1]=fP[2]; p[2]=fP[3];
  return Local2GlobalMomentum(p,fAlpha);
}

Bool_t AliExternalTrackParam::GetXYZ(Double_t *r) const {
  //---------------------------------------------------------------------
  // This function returns the global track position
  //---------------------------------------------------------------------
  r[0]=fX; r[1]=fP[0]; r[2]=fP[1];
  return Local2GlobalPosition(r,fAlpha);
}

Bool_t AliExternalTrackParam::GetCovarianceXYZPxPyPz(Double_t cv[21]) const {
  //---------------------------------------------------------------------
  // This function returns the global covariance matrix of the track params
  // 
  // Cov(x,x) ... :   cv[0]
  // Cov(y,x) ... :   cv[1]  cv[2]
  // Cov(z,x) ... :   cv[3]  cv[4]  cv[5]
  // Cov(px,x)... :   cv[6]  cv[7]  cv[8]  cv[9]
  // Cov(py,x)... :   cv[10] cv[11] cv[12] cv[13] cv[14]
  // Cov(pz,x)... :   cv[15] cv[16] cv[17] cv[18] cv[19] cv[20]
  //
  // Results for (nearly) straight tracks are meaningless !
  //---------------------------------------------------------------------
  if (TMath::Abs(fP[4])<=kAlmost0) {
     for (Int_t i=0; i<21; i++) cv[i]=0.;
     return kFALSE;
  }
  if (TMath::Abs(fP[2]) > kAlmost1) {
     for (Int_t i=0; i<21; i++) cv[i]=0.;
     return kFALSE;
  }
  Double_t pt=1./TMath::Abs(fP[4]);
  Double_t cs=TMath::Cos(fAlpha), sn=TMath::Sin(fAlpha);
  Double_t r=TMath::Sqrt(1-fP[2]*fP[2]);

  Double_t m00=-sn, m10=cs;
  Double_t m23=-pt*(sn + fP[2]*cs/r), m43=-pt*pt*(r*cs - fP[2]*sn);
  Double_t m24= pt*(cs - fP[2]*sn/r), m44=-pt*pt*(r*sn + fP[2]*cs);
  Double_t m35=pt, m45=-pt*pt*fP[3];

  cv[0 ] = fC[0]*m00*m00;
  cv[1 ] = fC[0]*m00*m10; 
  cv[2 ] = fC[0]*m10*m10;
  cv[3 ] = fC[1]*m00; 
  cv[4 ] = fC[1]*m10; 
  cv[5 ] = fC[2];
  cv[6 ] = m00*(fC[3]*m23 + fC[10]*m43); 
  cv[7 ] = m10*(fC[3]*m23 + fC[10]*m43); 
  cv[8 ] = fC[4]*m23 + fC[11]*m43; 
  cv[9 ] = m23*(fC[5]*m23 + fC[12]*m43)  +  m43*(fC[12]*m23 + fC[14]*m43);
  cv[10] = m00*(fC[3]*m24 + fC[10]*m44); 
  cv[11] = m10*(fC[3]*m24 + fC[10]*m44); 
  cv[12] = fC[4]*m24 + fC[11]*m44; 
  cv[13] = m23*(fC[5]*m24 + fC[12]*m44)  +  m43*(fC[12]*m24 + fC[14]*m44);
  cv[14] = m24*(fC[5]*m24 + fC[12]*m44)  +  m44*(fC[12]*m24 + fC[14]*m44);
  cv[15] = m00*(fC[6]*m35 + fC[10]*m45); 
  cv[16] = m10*(fC[6]*m35 + fC[10]*m45); 
  cv[17] = fC[7]*m35 + fC[11]*m45; 
  cv[18] = m23*(fC[8]*m35 + fC[12]*m45)  +  m43*(fC[13]*m35 + fC[14]*m45);
  cv[19] = m24*(fC[8]*m35 + fC[12]*m45)  +  m44*(fC[13]*m35 + fC[14]*m45); 
  cv[20] = m35*(fC[9]*m35 + fC[13]*m45)  +  m45*(fC[13]*m35 + fC[14]*m45);

  return kTRUE;
}


Bool_t 
AliExternalTrackParam::GetPxPyPzAt(Double_t x, Double_t b, Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the global track momentum extrapolated to
  // the radial position "x" (cm) in the magnetic field "b" (kG)
  //---------------------------------------------------------------------
  p[0]=fP[4]; 
  p[1]=fP[2]+(x-fX)*GetC(b); 
  p[2]=fP[3];
  return Local2GlobalMomentum(p,fAlpha);
}

Bool_t 
AliExternalTrackParam::GetYAt(Double_t x, Double_t b, Double_t &y) const {
  //---------------------------------------------------------------------
  // This function returns the local Y-coordinate of the intersection 
  // point between this track and the reference plane "x" (cm). 
  // Magnetic field "b" (kG)
  //---------------------------------------------------------------------
  Double_t dx=x-fX;
  if(TMath::Abs(dx)<=kAlmost0) {y=fP[0]; return kTRUE;}

  Double_t f1=fP[2], f2=f1 + dx*GetC(b);

  if (TMath::Abs(f1) >= kAlmost1) return kFALSE;
  if (TMath::Abs(f2) >= kAlmost1) return kFALSE;
  
  Double_t r1=TMath::Sqrt(1.- f1*f1), r2=TMath::Sqrt(1.- f2*f2);
  y = fP[0] + dx*(f1+f2)/(r1+r2);
  return kTRUE;
}

Bool_t 
AliExternalTrackParam::GetZAt(Double_t x, Double_t b, Double_t &z) const {
  //---------------------------------------------------------------------
  // This function returns the local Z-coordinate of the intersection 
  // point between this track and the reference plane "x" (cm). 
  // Magnetic field "b" (kG)
  //---------------------------------------------------------------------
  Double_t dx=x-fX;
  if(TMath::Abs(dx)<=kAlmost0) {z=fP[1]; return kTRUE;}

  Double_t f1=fP[2], f2=f1 + dx*fP[4]*b*kB2C;

  if (TMath::Abs(f1) >= kAlmost1) return kFALSE;
  if (TMath::Abs(f2) >= kAlmost1) return kFALSE;
  
  Double_t r1=sqrt(1.- f1*f1), r2=sqrt(1.- f2*f2);
  z = fP[1] + dx*(r2 + f2*(f1+f2)/(r1+r2))*fP[3]; // Many thanks to P.Hristov !
  return kTRUE;
}

Bool_t 
AliExternalTrackParam::GetXYZAt(Double_t x, Double_t b, Double_t *r) const {
  //---------------------------------------------------------------------
  // This function returns the global track position extrapolated to
  // the radial position "x" (cm) in the magnetic field "b" (kG)
  //---------------------------------------------------------------------
  Double_t dx=x-fX;
  if(TMath::Abs(dx)<=kAlmost0) return GetXYZ(r);

  Double_t f1=fP[2], f2=f1 + dx*GetC(b);

  if (TMath::Abs(f1) >= kAlmost1) return kFALSE;
  if (TMath::Abs(f2) >= kAlmost1) return kFALSE;
  
  Double_t r1=TMath::Sqrt(1.- f1*f1), r2=TMath::Sqrt(1.- f2*f2);
  r[0] = x;
  r[1] = fP[0] + dx*(f1+f2)/(r1+r2);
  r[2] = fP[1] + dx*(f1+f2)/(f1*r2 + f2*r1)*fP[3];
  return Local2GlobalPosition(r,fAlpha);
}

//_____________________________________________________________________________
void AliExternalTrackParam::Print(Option_t* /*option*/) const
{
// print the parameters and the covariance matrix

  printf("AliExternalTrackParam: x = %-12g  alpha = %-12g\n", fX, fAlpha);
  printf("  parameters: %12g %12g %12g %12g %12g\n",
	 fP[0], fP[1], fP[2], fP[3], fP[4]);
  printf("  covariance: %12g\n", fC[0]);
  printf("              %12g %12g\n", fC[1], fC[2]);
  printf("              %12g %12g %12g\n", fC[3], fC[4], fC[5]);
  printf("              %12g %12g %12g %12g\n", 
	 fC[6], fC[7], fC[8], fC[9]);
  printf("              %12g %12g %12g %12g %12g\n", 
	 fC[10], fC[11], fC[12], fC[13], fC[14]);
}

Double_t AliExternalTrackParam::GetSnpAt(Double_t x,Double_t b) const {
  //
  // Get sinus at given x
  //
  Double_t crv=GetC(b);
  if (TMath::Abs(b) < kAlmost0Field) crv=0.;
  Double_t dx = x-fX;
  Double_t res = fP[2]+dx*crv;
  return res;
}
