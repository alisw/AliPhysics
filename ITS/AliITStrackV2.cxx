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

///////////////////////////////////////////////////////////////////////////
//                Implementation of the ITS track class
//
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//     dEdx analysis by: Boris Batyunya, JINR, Boris.Batiounia@cern.ch
///////////////////////////////////////////////////////////////////////////
#include <TMath.h>

#include "AliCluster.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliITSReconstructor.h"
#include "AliITStrackV2.h"
#include "AliTracker.h"

const Int_t AliITStrackV2::fgkWARN = 5;

ClassImp(AliITStrackV2)


//____________________________________________________________________________
AliITStrackV2::AliITStrackV2() : AliKalmanTrack(),
  fdEdx(0),
  fESDtrack(0)
{
  for(Int_t i=0; i<2*AliITSgeomTGeo::kNLayers; i++) {fIndex[i]=-1; fModule[i]=-1;}
  for(Int_t i=0; i<4; i++) fdEdxSample[i]=0;
}


//____________________________________________________________________________
AliITStrackV2::AliITStrackV2(AliESDtrack& t,Bool_t c) throw (const Char_t *) :
  AliKalmanTrack(),
  fdEdx(t.GetITSsignal()),
  fESDtrack(&t)
{
  //------------------------------------------------------------------
  // Conversion ESD track -> ITS track.
  // If c==kTRUE, create the ITS track out of the constrained params.
  //------------------------------------------------------------------
  const AliExternalTrackParam *par=&t;
  if (c) {
    par=t.GetConstrainedParam();
    if (!par) throw "AliITStrackV2: conversion failed !\n";
  }
  Set(par->GetX(),par->GetAlpha(),par->GetParameter(),par->GetCovariance());

  //if (!Invariant()) throw "AliITStrackV2: conversion failed !\n";

  SetLabel(t.GetLabel());
  SetMass(t.GetMass());
  SetNumberOfClusters(t.GetITSclusters(fIndex));

  if (t.GetStatus()&AliESDtrack::kTIME) {
    StartTimeIntegral();
    Double_t times[10]; t.GetIntegratedTimes(times); SetIntegratedTimes(times);
    SetIntegratedLength(t.GetIntegratedLength());
  }

  for(Int_t i=0; i<4; i++) fdEdxSample[i]=0;
}

//____________________________________________________________________________
void AliITStrackV2::ResetClusters() {
  //------------------------------------------------------------------
  // Reset the array of attached clusters.
  //------------------------------------------------------------------
  for (Int_t i=0; i<2*AliITSgeomTGeo::kNLayers; i++) fIndex[i]=-1;
  SetChi2(0.); 
  SetNumberOfClusters(0);
} 

void AliITStrackV2::UpdateESDtrack(ULong_t flags) const {
  fESDtrack->UpdateTrackParams(this,flags);
  // copy the module indices
  for(Int_t i=0;i<12;i++) {
    //   printf("     %d\n",GetModuleIndex(i));
    fESDtrack->SetITSModuleIndex(i,GetModuleIndex(i));
  }
}

//____________________________________________________________________________
AliITStrackV2::AliITStrackV2(const AliITStrackV2& t) : 
  AliKalmanTrack(t),
  fdEdx(t.fdEdx),
  fESDtrack(t.fESDtrack) 
{
  //------------------------------------------------------------------
  //Copy constructor
  //------------------------------------------------------------------
  Int_t i;
  for (i=0; i<4; i++) fdEdxSample[i]=t.fdEdxSample[i];
  for (i=0; i<2*AliITSgeomTGeo::GetNLayers(); i++) {
    fIndex[i]=t.fIndex[i];
    fModule[i]=t.fModule[i];
  }
}

//_____________________________________________________________________________
Int_t AliITStrackV2::Compare(const TObject *o) const {
  //-----------------------------------------------------------------
  // This function compares tracks according to the their curvature
  //-----------------------------------------------------------------
  AliITStrackV2 *t=(AliITStrackV2*)o;
  //Double_t co=OneOverPt();
  //Double_t c =OneOverPt();
  Double_t co=t->GetSigmaY2()*t->GetSigmaZ2();
  Double_t c =GetSigmaY2()*GetSigmaZ2();
  if (c>co) return 1;
  else if (c<co) return -1;
  return 0;
}

//____________________________________________________________________________
Bool_t 
AliITStrackV2::PropagateToVertex(const AliESDVertex *v,Double_t d,Double_t x0) 
{
  //------------------------------------------------------------------
  //This function propagates a track to the minimal distance from the origin
  //------------------------------------------------------------------  
  Double_t bz=GetBz();
  if (PropagateToDCA(v,bz,kVeryBig)) {
    Double_t xOverX0,xTimesRho; 
    xOverX0 = d; xTimesRho = d*x0;
    if (CorrectForMeanMaterial(xOverX0,xTimesRho,kTRUE)) return kTRUE;
  }
  return kFALSE;
}

//____________________________________________________________________________
Bool_t AliITStrackV2::
GetGlobalXYZat(Double_t xloc, Double_t &x, Double_t &y, Double_t &z) const {
  //------------------------------------------------------------------
  //This function returns a track position in the global system
  //------------------------------------------------------------------
  Double_t r[3];
  Bool_t rc=GetXYZAt(xloc, GetBz(), r);
  x=r[0]; y=r[1]; z=r[2]; 
  return rc;
}

//_____________________________________________________________________________
Double_t AliITStrackV2::GetPredictedChi2(const AliCluster *c) const {
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  Double_t p[2]={c->GetY(), c->GetZ()};
  Double_t cov[3]={c->GetSigmaY2(), 0., c->GetSigmaZ2()};
  return AliExternalTrackParam::GetPredictedChi2(p,cov);
}

//____________________________________________________________________________
Bool_t AliITStrackV2::PropagateTo(Double_t xk, Double_t d, Double_t x0) {
  //------------------------------------------------------------------
  //This function propagates a track
  //------------------------------------------------------------------

  Double_t oldX=GetX(), oldY=GetY(), oldZ=GetZ();
  
  Double_t bz=GetBz();
  if (!AliExternalTrackParam::PropagateTo(xk,bz)) return kFALSE;
  Double_t xOverX0,xTimesRho; 
  xOverX0 = d; xTimesRho = d*x0;
  if (!CorrectForMeanMaterial(xOverX0,xTimesRho,kTRUE)) return kFALSE;

  Double_t x=GetX(), y=GetY(), z=GetZ();
  if (IsStartedTimeIntegral() && x>oldX) {
    Double_t l2 = (x-oldX)*(x-oldX) + (y-oldY)*(y-oldY) + (z-oldZ)*(z-oldZ);
    AddTimeStep(TMath::Sqrt(l2));
  }

  return kTRUE;
}

//____________________________________________________________________________
Bool_t AliITStrackV2::PropagateToTGeo(Double_t xToGo, Int_t nstep, Double_t &xOverX0, Double_t &xTimesRho, Bool_t addTime) {
  //-------------------------------------------------------------------
  //  Propagates the track to a reference plane x=xToGo in n steps.
  //  These n steps are only used to take into account the curvature.
  //  The material is calculated with TGeo. (L.Gaudichet)
  //-------------------------------------------------------------------
  
  Double_t startx = GetX(), starty = GetY(), startz = GetZ();
  Double_t sign = (startx<xToGo) ? -1.:1.;
  Double_t step = (xToGo-startx)/TMath::Abs(nstep);

  Double_t start[3], end[3], mparam[7], bz = GetBz();
  Double_t x = startx;
  
  for (Int_t i=0; i<nstep; i++) {
    
    GetXYZ(start);   //starting global position
    x += step;
    if (!GetXYZAt(x, bz, end)) return kFALSE;
    if (!AliExternalTrackParam::PropagateTo(x, bz)) return kFALSE;
    AliTracker::MeanMaterialBudget(start, end, mparam);
    if (mparam[1]<900000) {
      xTimesRho = sign*mparam[4]*mparam[0];
      xOverX0   = mparam[1];
      if (!AliExternalTrackParam::CorrectForMeanMaterial(xOverX0,
				      xTimesRho,GetMass())) return kFALSE;
    }
  }

  if (addTime && IsStartedTimeIntegral() && GetX()>startx) {
    Double_t l2 = ( (GetX()-startx)*(GetX()-startx) +
		    (GetY()-starty)*(GetY()-starty) +
		    (GetZ()-startz)*(GetZ()-startz) );
    AddTimeStep(TMath::Sqrt(l2));
  }

  return kTRUE;
}

//____________________________________________________________________________
Bool_t AliITStrackV2::Update(const AliCluster* c, Double_t chi2, Int_t index) 
{
  //------------------------------------------------------------------
  //This function updates track parameters
  //------------------------------------------------------------------
  Double_t p[2]={c->GetY(), c->GetZ()};
  Double_t cov[3]={c->GetSigmaY2(), 0., c->GetSigmaZ2()};

  if (!AliExternalTrackParam::Update(p,cov)) return kFALSE;

  Int_t n=GetNumberOfClusters();
  if (!Invariant()) {
     if (n>fgkWARN) AliWarning("Wrong invariant !");
     return kFALSE;
  }

  if (chi2<0) return kTRUE;

  // fill residuals for ITS+TPC tracks 
  if (fESDtrack->GetStatus()&AliESDtrack::kTPCin) {
    AliTracker::FillResiduals(this,p,cov,c->GetVolumeId());
  }

  fIndex[n]=index;
  SetNumberOfClusters(n+1);
  SetChi2(GetChi2()+chi2);

  return kTRUE;
}

Bool_t AliITStrackV2::Invariant() const {
  //------------------------------------------------------------------
  // This function is for debugging purpose only
  //------------------------------------------------------------------
  Int_t n=GetNumberOfClusters();

  // take into account the misalignment error
  Float_t maxMisalErrY2=0,maxMisalErrZ2=0;
  for (Int_t lay=0; lay<AliITSgeomTGeo::kNLayers; lay++) {
    maxMisalErrY2 = TMath::Max(maxMisalErrY2,AliITSReconstructor::GetRecoParam()->GetClusterMisalErrorY(lay));
    maxMisalErrZ2 = TMath::Max(maxMisalErrZ2,AliITSReconstructor::GetRecoParam()->GetClusterMisalErrorZ(lay));
  }
  maxMisalErrY2 *= maxMisalErrY2;
  maxMisalErrZ2 *= maxMisalErrZ2;
  // this is because when we reset before refitting, we multiply the
  // matrix by 10
  maxMisalErrY2 *= 10.; 
  maxMisalErrZ2 *= 10.;

  Double_t sP2=GetParameter()[2];
  if (TMath::Abs(sP2) >= kAlmost1){
     if (n>fgkWARN) Warning("Invariant","fP2=%f\n",sP2);
     return kFALSE;
  }
  Double_t sC00=GetCovariance()[0];
  if (sC00<=0 || sC00>(9.+maxMisalErrY2)) {
     if (n>fgkWARN) Warning("Invariant","fC00=%f\n",sC00); 
     return kFALSE;
  }
  Double_t sC11=GetCovariance()[2];
  if (sC11<=0 || sC11>(9.+maxMisalErrZ2)) {
     if (n>fgkWARN) Warning("Invariant","fC11=%f\n",sC11); 
     return kFALSE;
  }
  Double_t sC22=GetCovariance()[5];
  if (sC22<=0 || sC22>1.) {
     if (n>fgkWARN) Warning("Invariant","fC22=%f\n",sC22); 
     return kFALSE;
  }
  Double_t sC33=GetCovariance()[9];
  if (sC33<=0 || sC33>1.) {
     if (n>fgkWARN) Warning("Invariant","fC33=%f\n",sC33); 
     return kFALSE;
  }
  Double_t sC44=GetCovariance()[14];
  if (sC44<=0 /*|| sC44>6e-5*/) {
     if (n>fgkWARN) Warning("Invariant","fC44=%f\n",sC44);
     return kFALSE;
  }

  return kTRUE;
}

//____________________________________________________________________________
Bool_t AliITStrackV2::Propagate(Double_t alp,Double_t xk) {
  //------------------------------------------------------------------
  //This function propagates a track
  //------------------------------------------------------------------
  Double_t bz=GetBz();
  if (!AliExternalTrackParam::Propagate(alp,xk,bz)) return kFALSE;

  if (!Invariant()) {
    Int_t n=GetNumberOfClusters();
    if (n>fgkWARN) AliWarning("Wrong invariant !");
    return kFALSE;
  }

  return kTRUE;
}

Bool_t AliITStrackV2::MeanBudgetToPrimVertex(Double_t xyz[3], Double_t step, Double_t &d) const {

  //-------------------------------------------------------------------
  //  Get the mean material budget between the actual point and the
  //  primary vertex. (L.Gaudichet)
  //-------------------------------------------------------------------

  Double_t cs=TMath::Cos(GetAlpha()), sn=TMath::Sin(GetAlpha());
  Double_t vertexX = xyz[0]*cs + xyz[1]*sn;

  Int_t nstep = Int_t((GetX()-vertexX)/step);
  if (nstep<1) nstep = 1;
  step = (GetX()-vertexX)/nstep;

  //  Double_t mparam[7], densMean=0, radLength=0, length=0;
  Double_t mparam[7];
  Double_t p1[3], p2[3], x = GetX(), bz = GetBz();
  GetXYZ(p1);

  d=0.;

  for (Int_t i=0; i<nstep; i++) {
    x  += step;
    if (!GetXYZAt(x, bz, p2)) return kFALSE;
    AliTracker::MeanMaterialBudget(p1, p2, mparam);
    if (mparam[1]>900000) return kFALSE;
    d  += mparam[1];

    p1[0] = p2[0];
    p1[1] = p2[1];
    p1[2] = p2[2];
  }

  return kTRUE;
}

Bool_t AliITStrackV2::Improve(Double_t x0,Double_t xyz[3],Double_t ers[3]) {
  //------------------------------------------------------------------
  //This function improves angular track parameters
  //------------------------------------------------------------------
  Double_t cs=TMath::Cos(GetAlpha()), sn=TMath::Sin(GetAlpha());
//Double_t xv = xyz[0]*cs + xyz[1]*sn; // vertex
  Double_t yv =-xyz[0]*sn + xyz[1]*cs; // in the
  Double_t zv = xyz[2];                // local frame

  Double_t dy = Par(0) - yv, dz = Par(1) - zv;
  Double_t r2=GetX()*GetX() + dy*dy;
  Double_t p2=(1.+ GetTgl()*GetTgl())/(GetSigned1Pt()*GetSigned1Pt());
  Double_t beta2=p2/(p2 + GetMass()*GetMass());
  x0*=TMath::Sqrt((1.+ GetTgl()*GetTgl())/(1.- GetSnp()*GetSnp()));
  Double_t theta2=14.1*14.1/(beta2*p2*1e6)*x0;
  //Double_t theta2=1.0259e-6*14*14/28/(beta2*p2)*x0*9.36*2.33;

  Double_t cnv=GetBz()*kB2C;
  {
    Double_t dummy = 4/r2 - GetC()*GetC();
    if (dummy < 0) return kFALSE;
    Double_t parp = 0.5*(GetC()*GetX() + dy*TMath::Sqrt(dummy));
    Double_t sigma2p = theta2*(1.- GetSnp()*GetSnp())*(1. + GetTgl()*GetTgl());
    sigma2p += Cov(0)/r2*(1.- dy*dy/r2)*(1.- dy*dy/r2);
    sigma2p += ers[1]*ers[1]/r2;
    sigma2p += 0.25*Cov(14)*cnv*cnv*GetX()*GetX();
    Double_t eps2p=sigma2p/(Cov(5) + sigma2p);
    Par(0) += Cov(3)/(Cov(5) + sigma2p)*(parp - GetSnp());
    Par(2)  = eps2p*GetSnp() + (1 - eps2p)*parp;
    Cov(5) *= eps2p;
    Cov(3) *= eps2p;
  }
  {
    Double_t parl=0.5*GetC()*dz/TMath::ASin(0.5*GetC()*TMath::Sqrt(r2));
    Double_t sigma2l=theta2;
    sigma2l += Cov(2)/r2 + Cov(0)*dy*dy*dz*dz/(r2*r2*r2);
    sigma2l += ers[2]*ers[2]/r2;
    Double_t eps2l = sigma2l/(Cov(9) + sigma2l);
    Par(1) += Cov(7 )/(Cov(9) + sigma2l)*(parl - Par(3));
    Par(4) += Cov(13)/(Cov(9) + sigma2l)*(parl - Par(3));
    Par(3)  = eps2l*Par(3) + (1-eps2l)*parl;
    Cov(9) *= eps2l; 
    Cov(13)*= eps2l; 
    Cov(7) *= eps2l; 
  }
  if (!Invariant()) return kFALSE;

  return kTRUE;
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

//____________________________________________________________________________
Bool_t AliITStrackV2::
GetPhiZat(Double_t r, Double_t &phi, Double_t &z) const {
  //------------------------------------------------------------------
  // This function returns the global cylindrical (phi,z) of the track 
  // position estimated at the radius r. 
  // The track curvature is neglected.
  //------------------------------------------------------------------
  Double_t d=GetD(0.,0.);
  if (TMath::Abs(d) > r) return kFALSE; 

  Double_t rcurr=TMath::Sqrt(GetX()*GetX() + GetY()*GetY());
  if (TMath::Abs(d) > rcurr) return kFALSE; 
  Double_t phicurr=GetAlpha()+TMath::ASin(GetSnp());

  phi=phicurr+TMath::ASin(d/r)-TMath::ASin(d/rcurr);
  z=GetZ()+GetTgl()*(TMath::Sqrt((r-d)*(r+d))-TMath::Sqrt((rcurr-d)*(rcurr+d)));
  return kTRUE;
}
//____________________________________________________________________________
Bool_t AliITStrackV2::
GetLocalXat(Double_t r,Double_t &xloc) const {
  //------------------------------------------------------------------
  // This function returns the local x of the track 
  // position estimated at the radius r. 
  // The track curvature is neglected.
  //------------------------------------------------------------------
  Double_t d=GetD(0.,0.);
  if (TMath::Abs(d) > r) return kFALSE; 

  Double_t rcurr=TMath::Sqrt(GetX()*GetX() + GetY()*GetY());
  Double_t phicurr=GetAlpha()+TMath::ASin(GetSnp());
  Double_t phi=phicurr+TMath::ASin(d/r)-TMath::ASin(d/rcurr);

  xloc=r*(TMath::Cos(phi)*TMath::Cos(GetAlpha())
         +TMath::Sin(phi)*TMath::Sin(GetAlpha())); 

  return kTRUE;
}
