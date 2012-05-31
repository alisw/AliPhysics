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

///////////////////////////////////////////////////////////////////////////////
//
//                      Kalman-Filter-like fit 
//   to a straight-line crossing a set of arbitrarily oriented planes.
//   The resulting line is given by the equation:
//                  (x-x0)/vx = (y-y0)/1 = (z-z0)/vz
//   Parameters of the fit are:
//        x0,y0,z0 (a point on the line) and
//        vx,1,vz  (a vector collinear with the line)
//
//   LIMITATION:  The line must not be perpendicular to the Y axis. 
//
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//
//////////////////////////////////////////////////////////////////////////////

#include <TMatrixD.h>
#include <TMatrixDSym.h>

#include "AliLog.h"
#include "AliTrackFitterKalman.h"

ClassImp(AliTrackFitterKalman)

const Double_t AliTrackFitterKalman::fgkMaxChi2=33.;

AliTrackFitterKalman::
AliTrackFitterKalman(AliTrackPointArray *array, Bool_t owner):
  AliTrackFitter(array,owner),
  fMaxChi2(fgkMaxChi2)
{
  //
  // Constructor
  //
}

Bool_t AliTrackFitterKalman::Begin(Int_t first, Int_t last) {
  //
  // Make a seed out of the track points with the indices "first" and "last". 
  // This is the "default" seed.
  //
  fPoints->Sort();
  AliTrackPoint fp,lp;
  fPoints->GetPoint(fp,first); fPoints->GetPoint(lp,last);
  return MakeSeed(&fp,&lp);
}


Bool_t AliTrackFitterKalman::
MakeSeed(const AliTrackPoint *p0, const AliTrackPoint *p1) {
  //
  // Make a seed out of the two track points 
  //
  Double_t x0=p0->GetX(), y0=p0->GetY(), z0=p0->GetZ();
  Double_t dx=p1->GetX()-x0, dy=p1->GetY()-y0, dz=p1->GetZ()-z0;
  if (dy==0.) { 
     AliError("Seeds perpendicular to Y axis are not allowed !"); 
     return kFALSE; 
  }
  Double_t vx=dx/dy;
  Double_t vz=dz/dy;
  Double_t par[5]={x0,y0,z0,vx,vz};

  const Float_t *cv0=p0->GetCov();
  const Float_t *cv1=p1->GetCov();
  Double_t rdy2=(cv0[3]+cv1[3])/dy/dy;
  Double_t svx2=(cv0[0]+cv1[0])/dy/dy + vx*vx*rdy2;     
  Double_t svz2=(cv0[5]+cv1[5])/dy/dy + vz*vz*rdy2;     
  Double_t cov[15]={
     cv0[0],
     cv0[1],cv0[3],
     cv0[2],cv0[4],cv0[5],
     0.,0.,0.,svx2,
     0.,0.,0.,0.,svz2
  };

  SetSeed(par,cov);

  return kTRUE;
}

Bool_t AliTrackFitterKalman::AddPoint(const AliTrackPoint *p)
{
  //
  // Add a point to the fit
  //

  if (!Propagate(p))   return kFALSE;
  Double_t chi2=GetPredictedChi2(p);
  if (chi2>fMaxChi2)   return kFALSE;
  if (!Update(p,chi2)) return kFALSE;
  
  return kTRUE;
}


Double_t AliTrackFitterKalman::GetPredictedChi2(const AliTrackPoint *p) const {
  //
  // Calculate the predicted chi2 increment.
  //

  Float_t txyz[3]={GetParam()[0],GetParam()[1],GetParam()[2]};
  TMatrixDSym &cv=*fCov;
  Float_t tcov[6]={
    cv(0,0),cv(1,0),cv(2,0),
            cv(1,1),cv(2,1),
                    cv(2,2)
  };
  AliTrackPoint tp(txyz,tcov,p->GetVolumeID());

  Double_t chi2=p->GetResidual(tp,kTRUE);

  return chi2;
}


Bool_t AliTrackFitterKalman::Propagate(const AliTrackPoint *p) {
  //
  // Propagate the track towards the measured point "p"
  //

  TMatrixDSym &cv=*fCov;
  Double_t s=p->GetY() - fParams[1];
  Double_t sig2=s*s/12.; 

  Double_t vx=fParams[3], vz=fParams[4];
  fParams[0] += s*vx;
  fParams[1] += s;
  fParams[2] += s*vz;

  Double_t 
  c00 = cv(0,0) + 2*s*cv(3,0) + s*s*cv(3,3) + vx*vx*sig2,

  c10 = cv(1,0) + s*cv(1,3) + vx*sig2, c11=cv(1,1) + sig2,
  
  c20 = cv(2,0) + s*(cv(4,0) + cv(2,3)) + s*s*cv(4,3) + vx*vz*sig2,
  c21 = cv(2,1) + s*cv(4,1) + vz*sig2,
  c22 = cv(2,2) + 2*s*cv(4,2) + s*s*cv(4,4) + vz*vz*sig2,

  c30 = cv(3,0) + s*cv(3,3), c31 = cv(3,1),
  c32 = cv(3,2) + s*cv(3,4), c33 = cv(3,3),

  c40 = cv(4,0) + s*cv(4,3), c41 = cv(4,1),
  c42 = cv(4,2) + s*cv(4,4), c43 = cv(4,3), c44 = cv(4,4);


  cv(0,0)=c00; cv(0,1)=c10; cv(0,2)=c20; cv(0,3)=c30; cv(0,4)=c40;
  cv(1,0)=c10; cv(1,1)=c11; cv(1,2)=c21; cv(1,3)=c31; cv(1,4)=c41;
  cv(2,0)=c20; cv(2,1)=c21; cv(2,2)=c22; cv(2,3)=c32; cv(2,4)=c42;
  cv(3,0)=c30; cv(3,1)=c31; cv(3,2)=c32; cv(3,3)=c33; cv(3,4)=c43;
  cv(4,0)=c40; cv(4,1)=c41; cv(4,2)=c42; cv(4,3)=c43; cv(4,4)=c44;

  return kTRUE;
}

Bool_t AliTrackFitterKalman::Update(const AliTrackPoint *p, Double_t chi2) {
  //
  // Update the track params using the measured point "p"
  //

  TMatrixDSym &c=*fCov;
  const Float_t *cov=p->GetCov();

  TMatrixDSym v(3);
  v(0,0)=cov[0]+c(0,0); v(0,1)=cov[1]+c(0,1); v(0,2)=cov[2]+c(0,2);
  v(1,0)=cov[1]+c(1,0); v(1,1)=cov[3]+c(1,1); v(1,2)=cov[4]+c(1,2);
  v(2,0)=cov[2]+c(2,0); v(2,1)=cov[4]+c(2,1); v(2,2)=cov[5]+c(2,2);
  if (TMath::Abs(v.Determinant()) < 1.e-33) return kFALSE;
  v.Invert();

  TMatrixD ch(5,3);
  ch(0,0)=c(0,0); ch(0,1)=c(0,1); ch(0,2)=c(0,2);
  ch(1,0)=c(1,0); ch(1,1)=c(1,1); ch(1,2)=c(1,2);
  ch(2,0)=c(2,0); ch(2,1)=c(2,1); ch(2,2)=c(2,2);
  ch(3,0)=c(3,0); ch(3,1)=c(3,1); ch(3,2)=c(3,2);
  ch(4,0)=c(4,0); ch(4,1)=c(4,1); ch(4,2)=c(4,2);

  TMatrixD k(ch,TMatrixD::kMult,v);

  TMatrixD d(3,1);
  d(0,0) = p->GetX() - fParams[0];
  d(1,0) = p->GetY() - fParams[1];
  d(2,0) = p->GetZ() - fParams[2];

  TMatrixD x(k,TMatrixD::kMult,d);

  fParams[0]+=x(0,0);
  fParams[1]+=x(1,0);
  fParams[2]+=x(2,0);
  fParams[3]+=x(3,0);
  fParams[4]+=x(4,0);

  TMatrixD hc(3,5);
  hc(0,0)=c(0,0);hc(0,1)=c(0,1);hc(0,2)=c(0,2);hc(0,3)=c(0,3);hc(0,4)=c(0,4);
  hc(1,0)=c(1,0);hc(1,1)=c(1,1);hc(1,2)=c(1,2);hc(1,3)=c(1,3);hc(1,4)=c(1,4);
  hc(2,0)=c(2,0);hc(2,1)=c(2,1);hc(2,2)=c(2,2);hc(2,3)=c(2,3);hc(2,4)=c(2,4);
  
  TMatrixD s(k,TMatrixD::kMult,hc);

  c(0,0)-=s(0,0);c(0,1)-=s(0,1);c(0,2)-=s(0,2);c(0,3)-=s(0,3);c(0,4)-=s(0,4);
  c(1,0)-=s(1,0);c(1,1)-=s(1,1);c(1,2)-=s(1,2);c(1,3)-=s(1,3);c(1,4)-=s(1,4);
  c(2,0)-=s(2,0);c(2,1)-=s(2,1);c(2,2)-=s(2,2);c(2,3)-=s(2,3);c(2,4)-=s(2,4);
  c(3,0)-=s(3,0);c(3,1)-=s(3,1);c(3,2)-=s(3,2);c(3,3)-=s(3,3);c(3,4)-=s(3,4);
  c(4,0)-=s(4,0);c(4,1)-=s(4,1);c(4,2)-=s(4,2);c(4,3)-=s(4,3);c(4,4)-=s(4,4);

  fChi2 += chi2;
  fNdf  += 2;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliTrackFitterKalman::GetPCA(const AliTrackPoint &p, AliTrackPoint &i) const
{
  //
  // Get the intersection point "i" between the track and the plane
  // the point "p" belongs to.
  //

  TMatrixD t(3,1);
  Double_t s=p.GetY() - fParams[1];
  Double_t vx=fParams[3], vz=fParams[4];
  t(0,0) = fParams[0] + s*vx;
  t(1,0) = fParams[1] + s;
  t(2,0) = fParams[2] + s*vz;
  
  TMatrixDSym tC(3);
  {
  Double_t sig2=s*s/12.;

  tC(0,0) = vx*vx*sig2;
  tC(1,0) = vx*sig2; 
  tC(1,1) = sig2;
  tC(2,0) = vx*vz*sig2;
  tC(2,1) = vz*sig2;
  tC(2,2) = vz*vz*sig2;

  tC(0,1) = tC(1,0); tC(0,2) = tC(2,0);
                     tC(1,2) = tC(2,1);
  }

  TMatrixD m(3,1);
  m(0,0)=p.GetX();
  m(1,0)=p.GetY();
  m(2,0)=p.GetZ();
 
  TMatrixDSym mC(3);
  {
  const Float_t *cv=p.GetCov();
  mC(0,0)=cv[0]; mC(0,1)=cv[1]; mC(0,2)=cv[2];
  mC(1,0)=cv[1]; mC(1,1)=cv[3]; mC(1,2)=cv[4];
  mC(2,0)=cv[2]; mC(2,1)=cv[4]; mC(2,2)=cv[5];
  }

  TMatrixDSym tmW(tC);
  tmW+=mC;
  if (TMath::Abs(tmW.Determinant()) < 1.e-33) return kFALSE;
  tmW.Invert();

  TMatrixD mW(tC,TMatrixD::kMult,tmW);
  TMatrixD tW(mC,TMatrixD::kMult,tmW);

  TMatrixD mi(mW,TMatrixD::kMult,m);
  TMatrixD ti(tW,TMatrixD::kMult,t);
  ti+=mi;

  TMatrixD iC(tC,TMatrixD::kMult,tmW);
  iC*=mC;
  Double_t sig2=1.;  // Releasing the matrix by 1 cm along the track direction 
  iC(0,0) += vx*vx*sig2;
  iC(0,1) += vx*sig2; 
  iC(1,1) += sig2;
  iC(0,2) += vx*vz*sig2;
  iC(1,2) += vz*sig2;
  iC(2,2) += vz*vz*sig2;


  Float_t cov[6]={
    iC(0,0), iC(0,1), iC(0,2),
             iC(1,1), iC(1,2),
                      iC(2,2)
  };
  i.SetXYZ(ti(0,0),ti(1,0),ti(2,0),cov);
  UShort_t id=p.GetVolumeID();
  i.SetVolumeID(id);

  return kTRUE;
}


//_____________________________________________________________________________
void 
AliTrackFitterKalman::SetSeed(const Double_t par[5], const Double_t cov[15]) {
  //
  // Set the initial approximation for the track parameters
  //
  for (Int_t i=0; i<5; i++) fParams[i]=par[i];
  fParams[5]=0.;

  delete fCov;
  fCov=new TMatrixDSym(5);
  TMatrixDSym &cv=*fCov;

  cv(0,0)=cov[0 ];
  cv(1,0)=cov[1 ]; cv(1,1)=cov[2 ];
  cv(2,0)=cov[3 ]; cv(2,1)=cov[4 ]; cv(2,2)=cov[5 ];
  cv(3,0)=cov[6 ]; cv(3,1)=cov[7 ]; cv(3,2)=cov[8 ]; cv(3,3)=cov[9 ];
  cv(4,0)=cov[10]; cv(4,1)=cov[11]; cv(4,2)=cov[12]; cv(4,3)=cov[13]; cv(4,4)=cov[14];
 
  cv(0,1)=cv(1,0);
  cv(0,2)=cv(2,0); cv(1,2)=cv(2,1);
  cv(0,3)=cv(3,0); cv(1,3)=cv(3,1); cv(2,3)=cv(3,2);
  cv(0,4)=cv(4,0); cv(1,4)=cv(4,1); cv(2,4)=cv(4,2); cv(3,4)=cv(4,3);

}
