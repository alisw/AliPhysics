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


#include "TMatrixDSym.h"
//#include "TDecompChol.h"
#include "TMatrixD.h"
#include "AliRieman.h"

ClassImp(AliRieman)



AliRieman::AliRieman(){
  //
  // default constructor
  //
  for (Int_t i=0;i<6;i++) fParams[i] = 0;
  for (Int_t i=0;i<9;i++) fSumXY[i] = 0;
  for (Int_t i=0;i<9;i++) fSumXZ[i] = 0;

  fCapacity = 0;
  fN =0;
  fX  = 0;
  fY  = 0;
  fZ  = 0;
  fSy = 0;
  fSz = 0;
  fChi2  = 0;
  fChi2Y = 0;
  fChi2Z = 0;
  fCovar = 0;
  fConv = kFALSE;
}


AliRieman::AliRieman(Int_t capacity){
  //
  // default constructor
  //
  for (Int_t i=0;i<6;i++) fParams[i] = 0;
  for (Int_t i=0;i<9;i++) fSumXY[i] = 0;
  for (Int_t i=0;i<9;i++) fSumXZ[i] = 0;

  fCapacity = capacity;
  fN =0;
  fX  = new Double_t[fCapacity];
  fY  = new Double_t[fCapacity];
  fZ  = new Double_t[fCapacity];
  fSy = new Double_t[fCapacity];
  fSz = new Double_t[fCapacity];
  fCovar = new TMatrixDSym(6);
  fChi2  = 0;
  fChi2Y = 0;
  fChi2Z = 0;
  fConv = kFALSE;
}

AliRieman::AliRieman(const AliRieman &rieman):TObject(rieman){
  //
  // copy constructor
  // 
  for (Int_t i=0;i<6;i++) fParams[i] = rieman.fParams[i];
  for (Int_t i=0;i<9;i++) fSumXY[i]  = rieman.fSumXY[i];
  for (Int_t i=0;i<9;i++) fSumXZ[i]  = rieman.fSumXZ[i];
  fCapacity = rieman.fN;
  fN =rieman.fN;
  fCovar = new TMatrixDSym(*(rieman.fCovar)); 
  fConv = rieman.fConv;
  fX  = new Double_t[fN];
  fY  = new Double_t[fN];
  fZ  = new Double_t[fN];
  fSy = new Double_t[fN];
  fSz = new Double_t[fN];
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
  weight = 1./sz;
  fSumXZ[0] +=weight;
  fSumXZ[1] +=weight*x;   fSumXZ[2] +=weight*x*x; fSumXZ[3] +=weight*x*x*x; fSumXZ[4] += weight*x*x*x*x;
  fSumXZ[5] +=weight*z;   fSumXZ[6] +=weight*x*z; fSumXZ[7] +=weight*x*x*z;
  //
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
  TMatrixDSym     smatrixz(3);
  smatrixz(0,0) = fSumXZ[0]; smatrixz(0,1) = fSumXZ[1]; smatrixz(0,2) = fSumXZ[2];
  smatrixz(1,1) = fSumXZ[2]; smatrixz(1,2) = fSumXZ[3];
  smatrixz(2,2) = fSumXZ[4];
  smatrixz.Invert();
  if (smatrixz.IsValid()){
    sums(0,0)    = fSumXZ[5];
    sums(0,1)    = fSumXZ[6];
    sums(0,2)    = fSumXZ[7];
    TMatrixD res = sums*smatrixz;
    fParams[3] = res(0,0);
    fParams[4] = res(0,1);
    fParams[5] = res(0,2);
    for (Int_t i=0;i<3;i++)
      for (Int_t j=0;j<=i;j++){
	(*fCovar)(i+2,j+2)=smatrixz(i,j);
    }
    conv++;
  }

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


