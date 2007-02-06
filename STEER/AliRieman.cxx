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
// Class for global helix fit of a track
// Author: M.Ivanov
// The method uses the following idea:
//------------------------------------------------------
// XY direction
//
//  (x-x0)^2+(y-y0)^2-R^2=0 ===>
//
//  (x^2+y^2 -2*x*x0 - 2*y*y0+ x0^2 -y0^2 -R^2 =0;  ==>
//
//   substitution t = 1/(x^2+y^2),   u = 2*x*t, v = 2*y*t,  D0 = R^2 - x0^2- y0^2
//
//  1 - u*x0 - v*y0 - t *D0 =0 ;  - linear equation
//     
//  next substition   a = 1/y0    b = -x0/y0   c = -D0/y0
//
//  final linear equation :   a + u*b +t*c - v =0;
//
// Minimization :
//
// sum( (a + ui*b +ti*c - vi)^2)/(sigmai)^2 = min;
//
// where sigmai is the error of  maesurement  (a + ui*b +ti*c - vi)
//
// neglecting error of xi, and supposing  xi>>yi    sigmai ~ sigmaVi ~ 2*sigmay*t  


#include <TMatrixDSym.h>
#include <TMath.h>
#include <TMatrixD.h>

#include "AliRieman.h"

ClassImp(AliRieman)



AliRieman::AliRieman() :
  TObject(),
  fCapacity(0),
  fN(0),
  fX(0),
  fY(0),
  fZ(0),
  fSy(0),
  fSz(0),
  fCovar(0),
  fCovarPolY(0),
  fCovarPolZ(0),
  fSumZZ(0),
  fChi2(0),
  fChi2Y(0),
  fChi2Z(0),
  fConv(kFALSE)
{
  //
  // default constructor
  //
  for (Int_t i=0;i<6;i++) fParams[i] = 0;
  for (Int_t i=0;i<9;i++) fSumXY[i] = 0;
  for (Int_t i=0;i<9;i++) fSumXZ[i] = 0;
  for (Int_t i=0;i<6;i++) {
    fSumPolY[i]=0;
    fSumPolZ[i]=0;
  }
}


AliRieman::AliRieman(Int_t capacity) :
  TObject(),
  fCapacity(capacity),
  fN(0),
  fX(new Double_t[fCapacity]),
  fY(new Double_t[fCapacity]),
  fZ(new Double_t[fCapacity]),
  fSy(new Double_t[fCapacity]),
  fSz(new Double_t[fCapacity]),
  fCovar(new TMatrixDSym(6)),
  fCovarPolY(new TMatrixDSym(3)),
  fCovarPolZ(new TMatrixDSym(2)),
  fSumZZ(0),
  fChi2(0),
  fChi2Y(0),
  fChi2Z(0),
  fConv(kFALSE)
{
  //
  // default constructor
  //
  for (Int_t i=0;i<6;i++) fParams[i] = 0;
  for (Int_t i=0;i<9;i++) fSumXY[i] = 0;
  for (Int_t i=0;i<9;i++) fSumXZ[i] = 0;
  for (Int_t i=0;i<6;i++) {
    fSumPolY[i]=0;
    fSumPolZ[i]=0;
  }
}

AliRieman::AliRieman(const AliRieman &rieman):
  TObject(rieman),
  fCapacity(rieman.fN),
  fN(rieman.fN),
  fX(new Double_t[fN]),
  fY(new Double_t[fN]),
  fZ(new Double_t[fN]),
  fSy(new Double_t[fN]),
  fSz(new Double_t[fN]),
  fCovar(new TMatrixDSym(*(rieman.fCovar))), 
  fCovarPolY(new TMatrixDSym(*(rieman.fCovarPolY))), 
  fCovarPolZ(new TMatrixDSym(*(rieman.fCovarPolZ))), 
  fSumZZ(rieman.fSumZZ),
  fChi2(rieman.fChi2),
  fChi2Y(rieman.fChi2Y),
  fChi2Z(rieman.fChi2Z),
  fConv(rieman.fConv)

{
  //
  // copy constructor
  // 
  for (Int_t i=0;i<6;i++) fParams[i] = rieman.fParams[i];
  for (Int_t i=0;i<9;i++) fSumXY[i]  = rieman.fSumXY[i];
  for (Int_t i=0;i<9;i++) fSumXZ[i]  = rieman.fSumXZ[i];
  for (Int_t i=0;i<6;i++) {
    fSumPolY[i]=rieman.fSumPolY[i];
    fSumPolZ[i]=rieman.fSumPolZ[i];
  }
  for (Int_t i=0;i<rieman.fN;i++){
    fX[i]   = rieman.fX[i];
    fY[i]   = rieman.fY[i];
    fZ[i]   = rieman.fZ[i];
    fSy[i]  = rieman.fSy[i];
    fSz[i]  = rieman.fSz[i];
  }
}

AliRieman::~AliRieman()
{
  //
  // Destructor
  //
  delete[]fX;
  delete[]fY;
  delete[]fZ;
  delete[]fSy;
  delete[]fSz;
  delete fCovar;
  delete fCovarPolY;
  delete fCovarPolZ;
}

void AliRieman::Reset()
{
  //
  // Reset all the data members
  //
  fN=0;
  for (Int_t i=0;i<6;i++) fParams[i] = 0;
  for (Int_t i=0;i<9;i++) fSumXY[i] = 0;
  for (Int_t i=0;i<9;i++) fSumXZ[i] = 0; 
  for (Int_t i=0;i<6;i++) {
    fSumPolY[i]=0;
    fSumPolZ[i]=0;
  }
  fSumZZ =0;
  fConv =kFALSE;
}


void AliRieman::AddPoint(Double_t x, Double_t y, Double_t z, Double_t sy, Double_t sz)
{
  //
  //  Rieman update
  //
  //------------------------------------------------------
  // XY direction
  //
  //  (x-x0)^2+(y-y0)^2-R^2=0 ===>
  //
  //  (x^2+y^2 -2*x*x0 - 2*y*y0+ x0^2 -y0^2 -R^2 =0;  ==>
  //
  //   substitution t = 1/(x^2+y^2),   u = 2*x*t, v = 2*y*t,  D0 = R^2 - x0^2- y0^2
  //
  //  1 - u*x0 - v*y0 - t *D0 =0 ;  - linear equation
  //     
  //  next substition   a = 1/y0    b = -x0/y0   c = -D0/y0
  //
  //  final linear equation :   a + u*b +t*c - v =0;
  //
  // Minimization :
  //
  // sum( (a + ui*b +ti*c - vi)^2)/(sigmai)^2 = min;
  //
  // where sigmai is the error of  maesurement  (a + ui*b +ti*c - vi)
  //
  // neglecting error of xi, and supposing  xi>>yi    sigmai ~ sigmaVi ~ 2*sigmay*t  
  //
  if (fN==fCapacity-1) return;  // out of capacity
  fX[fN] = x; fY[fN]=y; fZ[fN]=z; fSy[fN]=sy; fSz[fN]=sz;
  //
  // XY part
  //
  Double_t  t  =  x*x+y*y;
  if (t<2) return;
  t            = 1./t;
  Double_t  u  =  2.*x*t;
  Double_t  v  =  2.*y*t;
  Double_t  error = 2.*sy*t;
  error *=error;
  Double_t weight = 1./error;
  fSumXY[0] +=weight;
  fSumXY[1] +=u*weight;      fSumXY[2]+=v*weight;  fSumXY[3]+=t*weight;
  fSumXY[4] +=u*u*weight;    fSumXY[5]+=t*t*weight;
  fSumXY[6] +=u*v*weight;    fSumXY[7]+=u*t*weight; fSumXY[8]+=v*t*weight;
  //
  // XZ part
  //
  weight = 1./(sz*sz);
  fSumXZ[0] +=weight;
  Double_t r = TMath::Sqrt(x*x+y*y);
  fSumXZ[1] +=weight*r;   fSumXZ[2] +=weight*r*r; fSumXZ[3] +=weight*z; fSumXZ[4] += weight*r*z;
  // Now the auxulary sums
  fSumXZ[5] +=weight*r*r*r/24; fSumXZ[6] +=weight*r*r*r*r/12; fSumXZ[7] +=weight*r*r*r*z/24;
  fSumXZ[8] +=weight*r*r*r*r*r*r/(24*24);
  fSumZZ += z*z*weight;
  //
  // sum accumulation for rough error estimates of the track extrapolation error
  //
  Double_t maty = 1./(sy*sy);
  Double_t matz = 1./(sz*sz);
  for (Int_t i=0;i<5; i++){
    fSumPolY[i] += maty;
    fSumPolZ[i] += matz;
    maty *= x;
    matz *= x;
  }
  fN++;
}


void AliRieman::UpdatePol(){
  //
  //  Rieman update
  //
  //
  if (fN==0) return;
  for (Int_t i=0;i<6;i++)fParams[i]=0;
  Int_t conv=0;
  //
  // XY part
  //
  TMatrixDSym     smatrix(3);
  TMatrixD        sums(1,3);
  //
  //   smatrix(0,0) = s0; smatrix(1,1)=su2; smatrix(2,2)=st2;
  //   smatrix(0,1) = su; smatrix(0,2)=st; smatrix(1,2)=sut;
  //   sums(0,0)    = sv; sums(0,1)=suv; sums(0,2)=svt;

  smatrix(0,0) = fSumXY[0]; smatrix(1,1)=fSumXY[4]; smatrix(2,2)=fSumXY[5];
  smatrix(0,1) = fSumXY[1]; smatrix(0,2)=fSumXY[3]; smatrix(1,2)=fSumXY[7];
  smatrix(1,0) = fSumXY[1]; smatrix(2,0)=fSumXY[3]; smatrix(2,1)=fSumXY[7];
  sums(0,0)    = fSumXY[2]; sums(0,1)   =fSumXY[6]; sums(0,2)   =fSumXY[8];
  smatrix.Invert();
  if (smatrix.IsValid()){
    for (Int_t i=0;i<3;i++)
      for (Int_t j=0;j<=i;j++){
	(*fCovar)(i,j)=smatrix(i,j);
      }
    TMatrixD  res = sums*smatrix;
    fParams[0] = res(0,0);
    fParams[1] = res(0,1);
    fParams[2] = res(0,2);
    conv++;
  }
  //
  // XZ part
  //
  Double_t rm1  = fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1); 
  
//   fSumXZ[1] += fSumXZ[5]*rm1*rm1;
//   fSumXZ[2] += fSumXZ[6]*rm1*rm1 + fSumXZ[8]*rm1*rm1*rm1*rm1;
//   fSumXZ[4] += fSumXZ[7]*rm1*rm1;
  Double_t sum1 = fSumXZ[1] + fSumXZ[5]*rm1*rm1;
  Double_t sum2 = fSumXZ[2] + fSumXZ[6]*rm1*rm1 + fSumXZ[8]*rm1*rm1*rm1*rm1;
  Double_t sum4 = fSumXZ[4] + fSumXZ[7]*rm1*rm1;
  
  TMatrixDSym     smatrixz(2);
  //  smatrixz(0,0) = fSumXZ[0]; smatrixz(0,1) = fSumXZ[1]; smatrixz(1,1) = fSumXZ[2];
  //smatrixz(1,0) = fSumXZ[1];
  smatrixz(0,0) = fSumXZ[0]; smatrixz(0,1) = sum1; smatrixz(1,1) = sum2;
  smatrixz(1,0) = sum1;
  smatrixz.Invert();
  TMatrixD        sumsxz(1,2);
  if (smatrixz.IsValid()){
    sumsxz(0,0)    = fSumXZ[3];
    //    sumsxz(0,1)    = fSumXZ[4];
    sumsxz(0,1)    = sum4;
    TMatrixD res = sumsxz*smatrixz;
    fParams[3] = res(0,0);
    fParams[4] = res(0,1);
    fParams[5] = 0;
    for (Int_t i=0;i<2;i++)
      for (Int_t j=0;j<=i;j++){
	(*fCovar)(i+3,j+3)=smatrixz(i,j);
      }
    conv++;
  }
  UpdateCovariancePol();
  //  (x-x0)^2+(y-y0)^2-R^2=0 ===>
  //
  //  (x^2+y^2 -2*x*x0 - 2*y*y0+ x0^2 -y0^2 -R^2 =0;  ==>
  //   substitution t = 1/(x^2+y^2),   u = 2*x*t, y = 2*y*t,  D0 = R^2 - x0^2- y0^2
  //
  //  1 - u*x0 - v*y0 - t *D0 =0 ;  - linear equation
  //     
  //  next substition   a = 1/y0    b = -x0/y0   c = -D0/y0
  //  final linear equation :   a + u*b +t*c - v =0;
  //
  //  fParam[0]  = 1/y0
  //  fParam[1]  = -x0/y0
  //  fParam[2]  = - (R^2 - x0^2 - y0^2)/y0
  if (conv>1) fConv =kTRUE;
  else
    fConv=kFALSE;
}

void AliRieman::Update(){
  //
  //  Rieman update
  //
  //
  if (fN==0) return;
  for (Int_t i=0;i<6;i++)fParams[i]=0;
  Int_t conv=0;
  //
  // XY part
  //
  TMatrixDSym     smatrix(3);
  TMatrixD        sums(1,3);
  //
  //   smatrix(0,0) = s0; smatrix(1,1)=su2; smatrix(2,2)=st2;
  //   smatrix(0,1) = su; smatrix(0,2)=st; smatrix(1,2)=sut;
  //   sums(0,0)    = sv; sums(0,1)=suv; sums(0,2)=svt;

  smatrix(0,0) = fSumXY[0]; smatrix(1,1)=fSumXY[4]; smatrix(2,2)=fSumXY[5];
  smatrix(0,1) = fSumXY[1]; smatrix(0,2)=fSumXY[3]; smatrix(1,2)=fSumXY[7];
  //
  smatrix(1,0) = fSumXY[1];
  smatrix(2,0) = fSumXY[3];
  smatrix(2,1) = fSumXY[7];

  sums(0,0)    = fSumXY[2]; sums(0,1)   =fSumXY[6]; sums(0,2)   =fSumXY[8];
  //TDecompChol decomp(smatrix);
  //decomp.SetTol(1.0e-14);
  //smatrix = 
  //decomp.Invert();
  smatrix.Invert();
  if (smatrix.IsValid()){
    for (Int_t i=0;i<3;i++)
      for (Int_t j=0;j<3;j++){
	(*fCovar)(i,j)=smatrix(i,j);
      }
    TMatrixD  res = sums*smatrix;
    fParams[0] = res(0,0);
    fParams[1] = res(0,1);
    fParams[2] = res(0,2);
    conv++;
  }
  if (conv==0){
    fConv = kFALSE;   //non converged
    return;
  }
  if (-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1.<0){
    fConv = kFALSE;   //non converged
    return;
  }
  //
  // XZ part
  //
  Double_t x0 = -fParams[1]/fParams[0];
  Double_t rm1  = fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1.); 
  Double_t sumXZ[9];

  for (Int_t i=0;i<9;i++) sumXZ[i]=0;
  for (Int_t i=0;i<fN;i++){
    Double_t phi  = TMath::ASin((fX[i]-x0)*rm1);
    Double_t phi0 = TMath::ASin((0.-x0)*rm1);
    Double_t weight = 1/fSz[i];
    weight *=weight;
    Double_t dphi = (phi-phi0)/rm1;
    sumXZ[0] +=weight;
    sumXZ[1] +=weight*dphi;
    sumXZ[2] +=weight*dphi*dphi;
    sumXZ[3] +=weight*fZ[i];
    sumXZ[4] +=weight*dphi*fZ[i];

  }

  TMatrixDSym     smatrixz(2);
  TMatrixD        sumsz(1,2);
  smatrixz(0,0) = sumXZ[0]; smatrixz(0,1) = sumXZ[1]; smatrixz(1,1) = sumXZ[2];
  smatrixz(1,0) = sumXZ[1];
  smatrixz.Invert();
  if (smatrixz.IsValid()){
    sumsz(0,0)    = sumXZ[3];
    sumsz(0,1)    = sumXZ[4];
    TMatrixD res = sumsz*smatrixz;
    fParams[3] = res(0,0);
    fParams[4] = res(0,1);
    for (Int_t i=0;i<2;i++)
      for (Int_t j=0;j<2;j++){
	(*fCovar)(i+3,j+3)=smatrixz(i,j);
    }
    conv++;
  }
  UpdateCovariancePol();
  //  (x-x0)^2+(y-y0)^2-R^2=0 ===>
  //
  //  (x^2+y^2 -2*x*x0 - 2*y*y0+ x0^2 -y0^2 -R^2 =0;  ==>
  //   substitution t = 1/(x^2+y^2),   u = 2*x*t, v = 2*y*t,  D0 = R^2 - x0^2- y0^2
  //
  //  1 - u*x0 - v*y0 - t *D0 =0 ;  - linear equation
  //     
  //  next substition   a = 1/y0    b = -x0/y0   c = -D0/y0
  //  final linear equation :   a + u*b +t*c - v =0;
  //
  //  fParam[0]  = 1/y0
  //  fParam[1]  = -x0/y0
  //  fParam[2]  = - (R^2 - x0^2 - y0^2)/y0
  if (conv>1) fConv =kTRUE;
  else
    fConv=kFALSE;
  fChi2Y = CalcChi2Y();
  fChi2Z = CalcChi2Z();
  fChi2  = fChi2Y +fChi2Z;
}

void AliRieman::UpdateCovariancePol(){
  //
  // covariance matrices for rough extrapolation error estimate
  // covariance matrix - get error estimates at point x = 0 
  // ! NOTE - error estimates is very rough and it is valid only if dl<<R
  // when dl is the distance between first and last point  
  //
  //
  (*fCovarPolY)(0,0) = fSumPolY[0]; (*fCovarPolY)(0,1) = fSumPolY[1]; (*fCovarPolY)(0,2) = fSumPolY[2];
  (*fCovarPolY)(1,0) = fSumPolY[1]; (*fCovarPolY)(1,1) = fSumPolY[2]; (*fCovarPolY)(1,2) = fSumPolY[3];
  (*fCovarPolY)(2,0) = fSumPolY[2]; (*fCovarPolY)(2,1) = fSumPolY[3]; (*fCovarPolY)(2,2) = fSumPolY[4];
  fCovarPolY->Invert();
  //

  (*fCovarPolZ)(0,0) = fSumPolZ[0]; (*fCovarPolZ)(0,1) = fSumPolZ[1];
  (*fCovarPolZ)(1,0) = fSumPolZ[1]; (*fCovarPolZ)(1,1) = fSumPolZ[2];
  fCovarPolZ->Invert();
  //
}

Double_t  AliRieman::GetErrY(Double_t x) const{
  //
  //    P0'  = P0 + P1 * x +  P2 * x^2 
  //    P1'  =      P1     +  P2 * x
  //    P2'  =             +  P2
  TMatrixD trans(3,3);
  trans(0,0) = 1;
  trans(0,1) = x;
  trans(0,2) = x*x;
  trans(1,1) = 1;
  trans(1,2) = x;
  trans(2,2) = 1;
  //
  TMatrixD covarX(trans,TMatrixD::kMult,*fCovarPolY);
  covarX*=trans.T();
  return covarX(0,0);
}

Double_t  AliRieman::GetErrZ(Double_t x) const{
  //
  //    assumption error of curvature determination neggligible
  //
  //    P0'  = P0 + P1 * x  
  //    P1'  =      P1    
  TMatrixD trans(2,2);
  trans(0,0) = 1;
  trans(0,1) = x;
  trans(1,1) = 1;
  //
  TMatrixD covarX(trans,TMatrixD::kMult,*fCovarPolZ);
  covarX*=trans.T();
  return covarX(0,0);
}

Bool_t AliRieman::GetExternalParameters(Double_t xref, Double_t *params, Double_t * covar){
  //
  // Get external parameters
  // + approximative covariance  
  //  0  1  2  3  4      // y   - local y
  //     5  6  7  8      // z   - local z
  //        9  10 11     // snp - local sin fi  
  //           12 13     // tgl - deep angle
  //              14     // C   - curvature
  Double_t sign = (GetC()>0) ? 1.:-1.;

  params[0] = GetYat(xref);
  params[1] = GetZat(xref);
  params[2] = TMath::Sin(TMath::ATan(GetDYat(xref)));
  params[3] = sign*fParams[4];
  params[4] = GetC();
  //
  // covariance matrix in y 
  //
  TMatrixD transY(3,3);
  transY(0,0) = 1;
  transY(0,1) = xref;
  transY(0,2) = xref*xref;
  transY(1,1) = 1;
  transY(1,2) = xref;
  transY(2,2) = 1;
  TMatrixD covarY(transY,TMatrixD::kMult,*fCovarPolY);
  covarY*=transY.T();
  //
  TMatrixD transZ(2,2);
  transZ(0,0) = 1;
  transZ(0,1) = xref;
  transZ(1,1) = 1;
  TMatrixD covarZ(transZ,TMatrixD::kMult,*fCovarPolZ);
  covarZ*=transZ.T();
  //
  // C ~ 2*P2 - in rotated frame
  //
  covar[0]  = covarY(0,0);  covar[1] = 0; covar[2] = covarY(0,1); covar[3] = 0; covar[4] = TMath::Sqrt(2.)*covarY(0,2);
  covar[5]  = covarZ(0,0);  covar[6] = 0; covar[7] = sign*covarZ(0,1); covar[8] = 0;
  covar[9]  = covarY(1,1);  covar[10]= 0; covar[11]= TMath::Sqrt(2.)*covarY(1,2);  
  covar[12] = covarZ(1,1);  covar[13]= 0; 
  covar[14] = 2.*covarY(2,2);
  //
  return fConv;
}




Double_t AliRieman::GetYat(Double_t x) const {
  //
  // get y at given x position
  //
  if (!fConv) return 0.;
  Double_t res2 = (x*fParams[0]+fParams[1]);
  res2*=res2;
  res2 = 1.-fParams[2]*fParams[0]+fParams[1]*fParams[1]-res2;
  if (res2<0) return 0;
  res2 = TMath::Sqrt(res2);
  res2 = (1-res2)/fParams[0];
  return res2;
}

Double_t AliRieman::GetDYat(Double_t x) const {
  //
  // get dy/dx at given x position
  //
  if (!fConv) return 0.;
  Double_t x0 = -fParams[1]/fParams[0];
  if (-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1<0) return 0;
  Double_t rm1  = fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1); 
  if ( 1./(rm1*rm1)-(x-x0)*(x-x0)<=0) return 0;
  Double_t res = (x-x0)/TMath::Sqrt(1./(rm1*rm1)-(x-x0)*(x-x0));
  if (fParams[0]<0) res*=-1.;
  return res;
}



Double_t AliRieman::GetZat(Double_t x) const {
  //
  // get z at given x position
  //
  if (!fConv) return 0.;
  Double_t x0 = -fParams[1]/fParams[0];
  if (-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1<=0) return 0;
  Double_t rm1  = fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1); 
  Double_t phi  = TMath::ASin((x-x0)*rm1);
  Double_t phi0 = TMath::ASin((0.-x0)*rm1);
  Double_t dphi = (phi-phi0);
  Double_t res = fParams[3]+fParams[4]*dphi/rm1;
  return res;
}

Double_t AliRieman::GetDZat(Double_t x) const {
  //
  // get dz/dx at given x postion
  //
  if (!fConv) return 0.;
  Double_t x0 = -fParams[1]/fParams[0]; 
  if  (-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1<=0) return 0;
  Double_t rm1  = fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1); 
  Double_t res = fParams[4]/TMath::Cos(TMath::ASin((x-x0)*rm1));
  return res;
}


//_____________________________________________________________________________
Bool_t AliRieman::GetXYZat(Double_t r, Double_t alpha, Float_t *xyz) const
{
  //
  // Returns position given radius
  //
  if (!fConv) return kFALSE;
  Double_t res2 = (r*fParams[0]+fParams[1]);
  res2*=res2;
  res2 = 1.-fParams[2]*fParams[0]+fParams[1]*fParams[1]-res2;
  if (res2<0) return kFALSE;
  res2 = TMath::Sqrt(res2);
  res2 = (1-res2)/fParams[0];

  if (!fConv) return kFALSE;
  Double_t x0 = -fParams[1]/fParams[0];
  if (-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1<=0) return 0;
  Double_t rm1  = fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1); 
  Double_t phi  = TMath::ASin((r-x0)*rm1);
  Double_t phi0 = TMath::ASin((0.-x0)*rm1);
  Double_t dphi = (phi-phi0);
  Double_t res = fParams[3]+fParams[4]*dphi/rm1;

  Double_t sin = TMath::Sin(alpha);
  Double_t cos = TMath::Cos(alpha);
  xyz[0] = r*cos - res2*sin;
  xyz[1] = res2*cos + r*sin;
  xyz[2] = res;

  return kTRUE;
}


Double_t AliRieman::GetC() const{
  //
  // get curvature
  //
  return fParams[0]/TMath::Sqrt(-fParams[2]*fParams[0]+fParams[1]*fParams[1]+1);
}

Double_t AliRieman::CalcChi2Y() const{ 
  //
  // calculate chi2 for Y
  //
  Double_t sumchi2 = 0;
  for (Int_t i=0;i<fN;i++){
    Double_t chi2 = (fY[i] - GetYat(fX[i]))/fSy[i];
    sumchi2+=chi2*chi2;    
  }  
  return sumchi2;
}


Double_t AliRieman::CalcChi2Z() const{
  //
  // calculate chi2 for Z
  //
  Double_t sumchi2 = 0;
  for (Int_t i=0;i<fN;i++){
    Double_t chi2 = (fZ[i] - GetZat(fX[i]))/fSy[i];
    sumchi2+=chi2*chi2;    
  }  
  return sumchi2;
}

Double_t AliRieman::CalcChi2() const{
  //
  // sum chi2 in both coord - supposing Y and Z independent
  //
  return CalcChi2Y()+CalcChi2Z();
}

AliRieman * AliRieman::MakeResiduals() const{
  //
  // create residual structure - ONLY for debug purposes
  //
  AliRieman *rieman = new AliRieman(fN);  
  for (Int_t i=0;i<fN;i++){
    rieman->AddPoint(fX[i],fY[i]-GetYat(fX[i]),fZ[i]-GetZat(fX[i]), fSy[i],fSz[i]);
  }
  return rieman;
}


