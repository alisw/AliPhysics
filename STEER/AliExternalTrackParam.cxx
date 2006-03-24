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
#include "AliKalmanTrack.h"
#include "AliTracker.h"


ClassImp(AliExternalTrackParam)

//_____________________________________________________________________________
AliExternalTrackParam::AliExternalTrackParam() :
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
AliExternalTrackParam::AliExternalTrackParam(Double_t x, Double_t alpha, 
					     const Double_t param[5], 
					     const Double_t covar[15]) :
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
AliExternalTrackParam::AliExternalTrackParam(const AliKalmanTrack& track) :
  fAlpha(track.GetAlpha())
{
  //
  //
  track.GetExternalParameters(fX,fP);
  track.GetExternalCovariance(fC);
}

//_____________________________________________________________________________
void AliExternalTrackParam::Set(const AliKalmanTrack& track) {
  //
  //
  fAlpha=track.GetAlpha();
  track.GetExternalParameters(fX,fP);
  track.GetExternalCovariance(fC);
}

//_____________________________________________________________________________
void AliExternalTrackParam::Reset() {
  fX=fAlpha=0.;
  for (Int_t i = 0; i < 5; i++) fP[i] = 0;
  for (Int_t i = 0; i < 15; i++) fC[i] = 0;
}

Double_t AliExternalTrackParam::GetP() const {
  //---------------------------------------------------------------------
  // This function returns the track momentum
  // Results for (nearly) straight tracks are meaningless !
  //---------------------------------------------------------------------
  if (TMath::Abs(fP[4])<=0) return 0;
  return TMath::Sqrt(1.+ fP[3]*fP[3])/TMath::Abs(fP[4]);
}

//_______________________________________________________________________
Double_t AliExternalTrackParam::GetD(Double_t b,Double_t x,Double_t y) const {
  //------------------------------------------------------------------
  // This function calculates the transverse impact parameter
  // with respect to a point with global coordinates (x,y)
  // in the magnetic field "b" (kG)
  //------------------------------------------------------------------
  Double_t rp4=kB2C*b*fP[4];

  Double_t xt=fX, yt=fP[0];

  Double_t sn=TMath::Sin(fAlpha), cs=TMath::Cos(fAlpha);
  Double_t a = x*cs + y*sn;
  y = -x*sn + y*cs; x=a;
  xt-=x; yt-=y;

  sn=rp4*xt - fP[2]; cs=rp4*yt + TMath::Sqrt(1.- fP[2]*fP[2]);
  a=2*(xt*fP[2] - yt*TMath::Sqrt(1.- fP[2]*fP[2]))-rp4*(xt*xt + yt*yt);
  if (rp4<0) a=-a;
  return a/(1 + TMath::Sqrt(sn*sn + cs*cs));
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

  return d;
}

Bool_t AliExternalTrackParam::Rotate(Double_t alpha) {
  //------------------------------------------------------------------
  // Transform this track to the local coord. system rotated
  // by angle "alpha" (rad) with respect to the global coord. system. 
  //------------------------------------------------------------------
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

  fAlpha = alpha;
  fX =  x*ca + fP0*sa;
  fP0= -x*sa + fP0*ca;
  fP2=  sf*ca - cf*sa;

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
  Double_t crv=kB2C*b*fP[4];
  Double_t dx=xk-fX;
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
  fP1 += dx*(f1+f2)/(f1*r2 + f2*r1)*fP3;
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
  if (TMath::Abs(p[0])<=0)        return kFALSE;
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
  if (TMath::Abs(fP[4])<=0) {
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
  p[1]=fP[2]+(x-fX)*fP[4]*b*kB2C; 
  p[2]=fP[3];
  return Local2GlobalMomentum(p,fAlpha);
}

Bool_t 
AliExternalTrackParam::GetXYZAt(Double_t x, Double_t b, Double_t *r) const {
  //---------------------------------------------------------------------
  // This function returns the global track position extrapolated to
  // the radial position "x" (cm) in the magnetic field "b" (kG)
  //---------------------------------------------------------------------
  Double_t dx=x-fX;
  Double_t f1=fP[2], f2=f1 + dx*fP[4]*b*kB2C;

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


Bool_t AliExternalTrackParam::PropagateTo(Double_t xToGo, Double_t mass, Double_t maxStep, Bool_t rotateTo){
  //----------------------------------------------------------------
  // Propagate this track to the plane X=xk (cm) 
  // correction for unhomogenity of the magnetic field and the
  // the correction for the material is included
  //
  //  Require acces to magnetic field and geomanager
  //
  // mass     - mass used in propagation - used for energy loss correction
  // maxStep  - maximal step for propagation
  //----------------------------------------------------------------
  const Double_t kEpsilon = 0.00001;
  Double_t xpos     = GetX();
  Double_t dir      = (xpos<xToGo) ? 1.:-1.;
  //
  while ( (xToGo-xpos)*dir > kEpsilon){
    Double_t step = dir*TMath::Min(TMath::Abs(xToGo-xpos), maxStep);
    Double_t x    = xpos+step;
    Double_t xyz0[3],xyz1[3],param[7];
    GetXYZ(xyz0);   //starting global position
    Float_t  pos0[3] = {xyz0[0],xyz0[1],xyz0[2]};
    Double_t magZ = AliTracker::GetBz(pos0);
    if (!GetXYZAt(x,magZ,xyz1)) return kFALSE;   // no prolongation
    AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);	
    if (!PropagateTo(x,magZ))  return kFALSE;
    Double_t distance = param[4];
    if (!CorrectForMaterial(distance,param[1],param[0],mass)) return kFALSE;
    if (rotateTo){
      GetXYZ(xyz0);   // global position
      Double_t alphan = TMath::ATan2(xyz0[1], xyz0[0]);
      if (!Rotate(alphan)) return kFALSE;
    }
    xpos = GetX();
  }
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliExternalTrackParam::CorrectForMaterial(Double_t d, Double_t x0, Double_t rho, Double_t mass)
{
  //
  // Take into account material effects assuming:
  // x0  - mean rad length
  // rho - mean density

  //
  // multiple scattering
  //
  if (mass<=0) {
    AliError("Non-positive mass");
    return kFALSE;
  }
  Double_t p2=(1.+ fP[3]*fP[3])/(fP[4]*fP[4]);
  Double_t beta2=p2/(p2 + mass*mass);
  Double_t theta2=14.1*14.1/(beta2*p2*1e6)*d/x0*rho;
  //
  fC[5] += theta2*(1.- fP[2]*fP[2])*(1. + fP[3]*fP[3]);
  fC[9] += theta2*(1. + fP[3]*fP[3])*(1. + fP[3]*fP[3]);
  fC[13] += theta2*fP[3]*fP[4]*(1. + fP[3]*fP[3]);
  fC[14] += theta2*fP[3]*fP[4]*fP[3]*fP[4];
  //
  Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2+1e-10)) - beta2)*d*rho;  
  fP[4] *=(1.- TMath::Sqrt(p2+mass*mass)/p2*dE);
  //
  Double_t sigmade = 0.02*TMath::Sqrt(TMath::Abs(dE));   // energy loss fluctuation 
  Double_t sigmac2 = sigmade*sigmade*fP[4]*fP[4]*(p2+mass*mass)/(p2*p2);
  fC[14] += sigmac2;
  return kTRUE;
}


