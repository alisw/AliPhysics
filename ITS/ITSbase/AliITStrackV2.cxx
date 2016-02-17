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
#include "AliESDVertex.h"
#include "AliITSReconstructor.h"
#include "AliITStrackV2.h"
#include "AliTracker.h"
#include "AliLog.h"
#include "AliPID.h"

const Int_t AliITStrackV2::fgkWARN = 5;

ClassImp(AliITStrackV2)


//____________________________________________________________________________
AliITStrackV2::AliITStrackV2() : AliKalmanTrack(),
  fCheckInvariant(kTRUE),
  fdEdx(0),
  fESDtrack(0)
{
  for(Int_t i=0; i<2*AliITSgeomTGeo::kNLayers; i++) {fIndex[i]=-1; fModule[i]=-1;}
  for(Int_t i=0; i<AliITSgeomTGeo::kNLayers; i++) {fSharedWeight[i]=0;}
  for(Int_t i=0; i<4; i++) fdEdxSample[i]=0;
}


//____________________________________________________________________________
AliITStrackV2::AliITStrackV2(AliESDtrack& t,Bool_t c):
  AliKalmanTrack(),
  fCheckInvariant(kTRUE),
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
    if (!par) AliError("AliITStrackV2: conversion failed !\n");
  }
  Set(par->GetX(),par->GetAlpha(),par->GetParameter(),par->GetCovariance());

  SetLabel(t.GetITSLabel());
  SetMass(t.GetMassForTracking());
  SetNumberOfClusters(t.GetITSclusters(fIndex));

  if (t.GetStatus()&AliESDtrack::kTIME) {
    StartTimeIntegral();
    Double_t times[AliPID::kSPECIESC]; 
    t.GetIntegratedTimes(times,AliPID::kSPECIESC); 
    SetIntegratedTimes(times);
    SetIntegratedLength(t.GetIntegratedLength());
  }

  for(Int_t i=0; i<AliITSgeomTGeo::kNLayers; i++) {fSharedWeight[i]=0;}
  for(Int_t i=0; i<4; i++) fdEdxSample[i]=0;
}

//____________________________________________________________________________
void AliITStrackV2::ResetClusters() {
  //------------------------------------------------------------------
  // Reset the array of attached clusters.
  //------------------------------------------------------------------
  for (Int_t i=0; i<2*AliITSgeomTGeo::kNLayers; i++) fIndex[i]=-1;
  for (Int_t i=0; i<AliITSgeomTGeo::kNLayers; i++) {fSharedWeight[i]=0;}
  SetChi2(0.); 
  SetNumberOfClusters(0);
} 

//____________________________________________________________________________
void AliITStrackV2::UpdateESDtrack(ULong_t flags) const {
  // Update track params
  fESDtrack->UpdateTrackParams(this,flags);
  //
  // set correctly the global label
  if (fESDtrack->IsOn(AliESDtrack::kTPCin)) { 
    // for global track the GetLabel should be negative if
    // 1) GetTPCLabel<0
    // 2) this->GetLabel()<0
    // 3) GetTPCLabel() != this->GetLabel()
    int label = fESDtrack->GetTPCLabel();
    int itsLabel = GetLabel();
    if (label<0 || itsLabel<0 || itsLabel!=label) label = -TMath::Abs(label);
    fESDtrack->SetLabel(label);
  }
  //
  // copy the module indices
  Int_t i;
  for(i=0;i<2*AliITSgeomTGeo::kNLayers;i++) {
    //   printf("     %d\n",GetModuleIndex(i));
    fESDtrack->SetITSModuleIndex(i,GetModuleIndex(i));
  }
  // copy the map of shared clusters
  if(flags==AliESDtrack::kITSin) {
    UChar_t itsSharedMap=0;
    for(i=0;i<AliITSgeomTGeo::kNLayers;i++) {
      if(fSharedWeight[i]>0) SETBIT(itsSharedMap,i);
      
    }
    fESDtrack->SetITSSharedMap(itsSharedMap);
  }

  // copy the 4 dedx samples
  Double_t sdedx[4]={0.,0.,0.,0.};
  for(i=0; i<4; i++) sdedx[i]=fdEdxSample[i];
  fESDtrack->SetITSdEdxSamples(sdedx);
}

//____________________________________________________________________________
AliITStrackV2::AliITStrackV2(const AliITStrackV2& t) : 
  AliKalmanTrack(t),
  fCheckInvariant(t.fCheckInvariant),
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
  for (i=0; i<AliITSgeomTGeo::kNLayers; i++) {fSharedWeight[i]=t.fSharedWeight[i];}
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
  //  Double_t bz=GetBz(); // RS: in the ITS constant field can be used
  Double_t bz=AliTrackerBase::GetBz();
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
  //  Bool_t rc=GetXYZAt(xloc, GetBz(), r);
  Bool_t rc=GetXYZAt(xloc, AliTrackerBase::GetBz(), r);
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
  const double kXThresh = 50.;
  Double_t oldX=GetX(), oldY=GetY(), oldZ=GetZ();

  if (oldX>kXThresh || xk>kXThresh) {
    Double_t b[3]; GetBxByBz(b);
    if (!AliExternalTrackParam::PropagateToBxByBz(xk,b)) return kFALSE;
  }
  else {
    Double_t bz=AliTrackerBase::GetBz();
    if (!AliExternalTrackParam::PropagateTo(xk,bz)) return kFALSE; //RS revert fast propagation
  }
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

  Double_t start[3], end[3], mparam[7];
  //Double_t bz = GetBz();
  Double_t b[3]; GetBxByBz(b);
  Double_t bz = b[2];

  Double_t x = startx;
  
  for (Int_t i=0; i<nstep; i++) {
    
    GetXYZ(start);   //starting global position
    x += step;
    //    if (!GetXYZAt(x, bz, end)) return kFALSE;
    //if (!AliExternalTrackParam::PropagateTo(x, bz)) return kFALSE;
    if (!AliExternalTrackParam::PropagateToBxByBz(x, b)) return kFALSE;
    GetXYZ(end);
    AliTracker::MeanMaterialBudget(start, end, mparam);
    xTimesRho = sign*mparam[4]*mparam[0];
    xOverX0   = mparam[1];
    if (mparam[1]<900000) {
      if (!AliExternalTrackParam::CorrectForMeanMaterial(xOverX0,
			   xTimesRho,GetMass())) return kFALSE;
    } else { // this happens when MeanMaterialBudget cannot cross a boundary
      return kFALSE;
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
  Double_t cov[3]={c->GetSigmaY2(), c->GetSigmaYZ(), c->GetSigmaZ2()};

  if (!AliExternalTrackParam::Update(p,cov)) return kFALSE;

  Int_t n=GetNumberOfClusters();
  if (!Invariant()) {
    if (n>fgkWARN) AliDebug(1,"Wrong invariant !");
     return kFALSE;
  }

  if (chi2<0) return kTRUE;

  // fill residuals for ITS+TPC tracks 
  if (fESDtrack) {
    if (fESDtrack->GetStatus()&AliESDtrack::kTPCin) {
      AliTracker::FillResiduals(this,p,cov,c->GetVolumeId());
    }
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
  if(!fCheckInvariant) return kTRUE;

  Int_t n=GetNumberOfClusters();
  static Float_t bz = GetBz();
  // take into account the misalignment error
  Float_t maxMisalErrY2=0,maxMisalErrZ2=0;
  //RS
  const AliITSRecoParam* recopar = AliITSReconstructor::GetRecoParam();
  if (!recopar) recopar = AliITSRecoParam::GetHighFluxParam();

  for (Int_t lay=0; lay<AliITSgeomTGeo::kNLayers; lay++) {
    maxMisalErrY2 = TMath::Max(maxMisalErrY2,recopar->GetClusterMisalErrorY(lay,bz));
    maxMisalErrZ2 = TMath::Max(maxMisalErrZ2,recopar->GetClusterMisalErrorZ(lay,bz));
  }
  maxMisalErrY2 *= maxMisalErrY2;
  maxMisalErrZ2 *= maxMisalErrZ2;
  // this is because when we reset before refitting, we multiply the
  // matrix by 10
  maxMisalErrY2 *= 10.; 
  maxMisalErrZ2 *= 10.;

  Double_t sP2=GetParameter()[2];
  if (TMath::Abs(sP2) >= kAlmost1){
    if (n>fgkWARN) AliDebug(1,Form("fP2=%f\n",sP2));
     return kFALSE;
  }
  Double_t sC00=GetCovariance()[0];
  if (sC00<=0 || sC00>(9.+maxMisalErrY2)) {
    if (n>fgkWARN) AliDebug(1,Form("fC00=%f\n",sC00)); 
     return kFALSE;
  }
  Double_t sC11=GetCovariance()[2];
  if (sC11<=0 || sC11>(9.+maxMisalErrZ2)) {
    if (n>fgkWARN) AliDebug(1,Form("fC11=%f\n",sC11)); 
     return kFALSE;
  }
  Double_t sC22=GetCovariance()[5];
  if (sC22<=0 || sC22>1.) {
    if (n>fgkWARN) AliDebug(1,Form("fC22=%f\n",sC22)); 
     return kFALSE;
  }
  Double_t sC33=GetCovariance()[9];
  if (sC33<=0 || sC33>1.) {
    if (n>fgkWARN) AliDebug(1,Form("fC33=%f\n",sC33)); 
     return kFALSE;
  }
  Double_t sC44=GetCovariance()[14];
  if (sC44<=0 /*|| sC44>6e-5*/) {
    if (n>fgkWARN) AliDebug(1,Form("fC44=%f\n",sC44));
     return kFALSE;
  }

  return kTRUE;
}

//____________________________________________________________________________
Bool_t AliITStrackV2::Propagate(Double_t alp,Double_t xk) {
  //------------------------------------------------------------------
  //This function propagates a track
  //------------------------------------------------------------------
  //Double_t bz=GetBz();
  //if (!AliExternalTrackParam::Propagate(alp,xk,bz)) return kFALSE;
  Double_t b[3]; GetBxByBz(b);
  if (!AliExternalTrackParam::PropagateBxByBz(alp,xk,b)) return kFALSE;

  if (!Invariant()) {
    Int_t n=GetNumberOfClusters();
    if (n>fgkWARN) AliDebug(1,"Wrong invariant !");
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
  //Store the initail track parameters

  Double_t x = GetX();
  Double_t alpha = GetAlpha();
  Double_t par[] = {GetY(),GetZ(),GetSnp(),GetTgl(),GetSigned1Pt()};
  Double_t cov[] = {
    GetSigmaY2(),
    GetSigmaZY(),
    GetSigmaZ2(),
    GetSigmaSnpY(),
    GetSigmaSnpZ(),
    GetSigmaSnp2(),
    GetSigmaTglY(),
    GetSigmaTglZ(),
    GetSigmaTglSnp(),
    GetSigmaTgl2(),
    GetSigma1PtY(),
    GetSigma1PtZ(),
    GetSigma1PtSnp(),
    GetSigma1PtTgl(),
    GetSigma1Pt2()
  }; 


  Double_t cs=TMath::Cos(GetAlpha()), sn=TMath::Sin(GetAlpha());
  Double_t xv = xyz[0]*cs + xyz[1]*sn; // vertex
  Double_t yv =-xyz[0]*sn + xyz[1]*cs; // in the
  Double_t zv = xyz[2];                // local frame

  Double_t dx = x - xv, dy = par[0] - yv, dz = par[1] - zv;
  Double_t r2=dx*dx + dy*dy;
  Double_t p2=(1.+ GetTgl()*GetTgl())/(GetSigned1Pt()*GetSigned1Pt());
  if (GetMass()<0) p2 *= 4; // q=2
  Double_t beta2=p2/(p2 + GetMass()*GetMass());
  x0*=TMath::Sqrt((1.+ GetTgl()*GetTgl())/(1.- GetSnp()*GetSnp()));
  Double_t theta2=14.1*14.1/(beta2*p2*1e6)*x0;
  //Double_t theta2=1.0259e-6*14*14/28/(beta2*p2)*x0*9.36*2.33;

  Double_t bz=GetBz();
  Double_t cnv=bz*kB2C;
  Double_t curv=GetC(bz);
  {
    Double_t dummy = 4/r2 - curv*curv;
    if (dummy < 0) return kFALSE;
    Double_t parp = 0.5*(curv*dx + dy*TMath::Sqrt(dummy));
    Double_t sigma2p = theta2*(1.-GetSnp())*(1.+GetSnp())*(1. + GetTgl()*GetTgl());
    Double_t ovSqr2 = 1./TMath::Sqrt(r2);
    Double_t tfact = ovSqr2*(1.-dy*ovSqr2)*(1.+dy*ovSqr2);
    sigma2p += cov[0]*tfact*tfact;
    sigma2p += ers[1]*ers[1]/r2;
    sigma2p += 0.25*cov[14]*cnv*cnv*dx*dx;
    Double_t eps2p=sigma2p/(cov[5] + sigma2p);
    par[0] += cov[3]/(cov[5] + sigma2p)*(parp - GetSnp());
    par[2]  = eps2p*GetSnp() + (1 - eps2p)*parp;
    cov[5] *= eps2p;
    cov[3] *= eps2p;
  }
  {
    Double_t parl=0.5*curv*dz/TMath::ASin(0.5*curv*TMath::Sqrt(r2));
    Double_t sigma2l=theta2;
    sigma2l += cov[2]/r2 + cov[0]*dy*dy*dz*dz/(r2*r2*r2);
    sigma2l += ers[2]*ers[2]/r2;
    Double_t eps2l = sigma2l/(cov[9] + sigma2l);
    par[1] += cov[7 ]/(cov[9] + sigma2l)*(parl - par[3]);
    par[4] += cov[13]/(cov[9] + sigma2l)*(parl - par[3]);
    par[3]  = eps2l*par[3] + (1-eps2l)*parl;
    cov[9] *= eps2l; 
    cov[13]*= eps2l; 
    cov[7] *= eps2l; 
  }

  Set(x,alpha,par,cov);

  if (!Invariant()) return kFALSE;

  return kTRUE;
}

void AliITStrackV2::CookdEdx(Double_t /*low*/, Double_t /*up*/) {
  //-----------------------------------------------------------------
  // This function calculates dE/dX within the "low" and "up" cuts.
  // Origin: Boris Batyunya, JINR, Boris.Batiounia@cern.ch 
  // Updated: F. Prino 8-June-2009
  //-----------------------------------------------------------------
  // The cluster order is: SDD-1, SDD-2, SSD-1, SSD-2

  Int_t nc=0;
  Float_t dedx[4];
  for (Int_t il=0; il<4; il++) { // count good (>0) dE/dx values
    if(fdEdxSample[il]>0.){
      dedx[nc]= fdEdxSample[il];
      nc++;
    }
  }
  if(nc<1){
    SetdEdx(0.);
    return;
  }

  Int_t swap; // sort in ascending order
  do {
    swap=0;
    for (Int_t i=0; i<nc-1; i++) {
      if (dedx[i]<=dedx[i+1]) continue;
      Float_t tmp=dedx[i];
      dedx[i]=dedx[i+1]; 
      dedx[i+1]=tmp;
      swap++;
    }
  } while (swap);


  Double_t sumamp=0,sumweight=0;
  Double_t weight[4]={1.,1.,0.,0.};
  if(nc==3) weight[1]=0.5;
  else if(nc<3) weight[1]=0.;
  for (Int_t i=0; i<nc; i++) {
    sumamp+= dedx[i]*weight[i];
    sumweight+=weight[i];
  }
  SetdEdx(sumamp/sumweight);
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
  if (TMath::Abs(d) > r) {
    if (r>1e-1) return kFALSE;
    r = TMath::Abs(d);
  }

  Double_t rcurr=TMath::Sqrt(GetX()*GetX() + GetY()*GetY());
  if (TMath::Abs(d) > rcurr) return kFALSE;
  Double_t globXYZcurr[3]; GetXYZ(globXYZcurr); 
  Double_t phicurr=TMath::ATan2(globXYZcurr[1],globXYZcurr[0]);

  if (GetX()>=0.) {
    phi=phicurr+TMath::ASin(d/r)-TMath::ASin(d/rcurr);
  } else {
    phi=phicurr+TMath::ASin(d/r)+TMath::ASin(d/rcurr)-TMath::Pi();
  }

  // return a phi in [0,2pi[ 
  if (phi<0.) phi+=2.*TMath::Pi();
  else if (phi>=2.*TMath::Pi()) phi-=2.*TMath::Pi();
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
  if (TMath::Abs(d) > r) { 
    if (r>1e-1) return kFALSE; 
    r = TMath::Abs(d); 
  } 

  Double_t rcurr=TMath::Sqrt(GetX()*GetX() + GetY()*GetY());
  Double_t globXYZcurr[3]; GetXYZ(globXYZcurr); 
  Double_t phicurr=TMath::ATan2(globXYZcurr[1],globXYZcurr[0]);
  Double_t phi;
  if (GetX()>=0.) {
    phi=phicurr+TMath::ASin(d/r)-TMath::ASin(d/rcurr);
  } else {
    phi=phicurr+TMath::ASin(d/r)+TMath::ASin(d/rcurr)-TMath::Pi();
  }

  xloc=r*(TMath::Cos(phi)*TMath::Cos(GetAlpha())
         +TMath::Sin(phi)*TMath::Sin(GetAlpha())); 

  return kTRUE;
}

//____________________________________________________________________________
Bool_t AliITStrackV2::ImproveKalman(Double_t xyz[3],Double_t ers[3], const Double_t* xlMS, const Double_t* x2X0MS, Int_t nMS)
{
  // Substitute the state of the track (p_{k|k},C_{k|k}) at the k-th measumerent by its
  // smoothed value from the k-th measurement + measurement at the vertex.
  // Account for the MS on nMS layers at x-postions xlMS with x/x0 = x2X0MS
  // p_{k|kv} = p_{k|k} + C_{k|k}*D^Tr_{k+1} B^{-1}_{k+1} ( vtx - D_{k+1}*p_{k|k})
  // C_{k|kv} = C_{k|k}*( I - D^Tr_{k+1} B^{-1}_{k+1} D_{k+1} C^Tr_{k|k})
  // 
  // where D_{k} = H_{k} F_{k} with H being the matrix converting the tracks parameters
  // to measurements m_{k} = H_{k} p_{k} and F_{k} the matrix propagating the track between the
  // the point k-1 and k:  p_{k|k-1} = F_{k} p_{k-1|k-1}
  //
  // B_{k+1} = V_{k+1} + H_{k+1} C_{k+1|k} H^Tr_{k+1} with V_{k+1} being the error of the measurment
  // at point k+1 (i.e. vertex), and C_{k+1|k} - error matrix extrapolated from k-th measurement to
  // k+1 (vtx) and accounting for the MS inbetween
  //
  // H = {{1,0,0,0,0},{0,1,0,0,0}}
  //
  double covc[15], *cori = (double*) GetCovariance(),par[5] = {GetY(),GetZ(),GetSnp(),GetTgl(),GetSigned1Pt()},
    &c00=cori[0],
    &c01=cori[1],&c11=cori[2],
    &c02=cori[3],&c12=cori[4],&c22=cori[5],
    &c03=cori[6],&c13=cori[7],&c23=cori[8],&c33=cori[9],
    &c04=cori[10],&c14=cori[11],&c24=cori[12],&c34=cori[13],&c44=cori[14],
    // for smoothed cov matrix 
    &cov00=covc[0],
    &cov01=covc[1],&cov11=covc[2],
    &cov02=covc[3],&cov12=covc[4],&cov22=covc[5],
    &cov03=covc[6],&cov13=covc[7],&cov23=covc[8],&cov33=covc[9],
    &cov04=covc[10],&cov14=covc[11],&cov24=covc[12],&cov34=covc[13],&cov44=covc[14];
  //
  double x = GetX(), alpha = GetAlpha();
  // vertex in the track frame
  double cs=TMath::Cos(alpha), sn=TMath::Sin(alpha);
  double xv = xyz[0]*cs + xyz[1]*sn, yv =-xyz[0]*sn + xyz[1]*cs, zv = xyz[2];
  double dx = xv - GetX();
  if (TMath::Abs(dx)<=kAlmost0)  return kTRUE;
  //
  double cnv=GetBz()*kB2C, x2r=cnv*par[4]*dx, f1=par[2], f2=f1+x2r;
  if (TMath::Abs(f1) >= kAlmost1 || TMath::Abs(f2) >= kAlmost1) {
    AliInfo(Form("Fail: %+e %+e",f1,f2));
    return kFALSE;
  }
  double r1=TMath::Sqrt((1.-f1)*(1.+f1)), r2=TMath::Sqrt((1.-f2)*(1.+f2)), dx2r=dx/(r1+r2);
  // elements of matrix F_{k+1} (1s on diagonal)
  double f02 = 2*dx2r, f04 = cnv*dx*dx2r, f13/*, f24 = cnv*dx*/;
  if (TMath::Abs(x2r)<0.05) f13 = dx*r2+f2*(f1+f2)*dx2r; // see AliExternalTrackParam::PropagateTo
  else {
    double dy2dx = (f1+f2)/(r1+r2);
    f13 = 2*TMath::ASin(0.5*TMath::Sqrt(1+dy2dx*dy2dx)*x2r)/(cnv*par[4]);
  }  
  // elements of matrix D_{k+1} = H_{k+1} * F_{k+1}
  // double d00 = 1., d11 = 1.;
  double &d02 = f02,  &d04 = f04, &d13 = f13;
  //
  // elements of matrix DC = D_{k+1}*C_{kk}^T
  double dc00 = c00+c02*d02+c04*d04,   dc10 = c01+c03*d13;
  double dc01 = c01+c12*d02+c14*d04,   dc11 = c11+c13*d13;
  double dc02 = c02+c22*d02+c24*d04,   dc12 = c12+c23*d13;
  double dc03 = c03+c23*d02+c34*d04,   dc13 = c13+c33*d13;
  double dc04 = c04+c24*d02+c44*d04,   dc14 = c14+c34*d13;
  //
  // difference between the vertex and the the track extrapolated to vertex
  yv -= par[0] + par[2]*d02 + par[4]*d04;
  zv -= par[1] + par[3]*d13;
  //
  // y,z part of the cov.matrix extrapolated to vtx (w/o MS contribution)
  // C_{k+1,k} = H F_{k+1} C_{k,k} F^Tr_{k+1} H^Tr = D C D^Tr
  double cv00 = dc00+dc02*d02+dc04*d04, cv01 = dc01+dc03*d13, cv11 = dc11+dc13*d13;
  //
  // add MS contribution layer by layer
  double xCurr = x;
  double p2Curr = par[2];
  //
  // precalculated factors of MS contribution matrix:
  double ms22t = (1. + par[3]*par[3]);
  double ms33t = ms22t*ms22t;
  double p34 = par[3]*par[4];
  double ms34t = p34*ms22t;
  double ms44t = p34*p34;
  //
  double p2=(1.+ par[3]*par[3])/(par[4]*par[4]);
  if (GetMass()<0) p2 *= 4; // q=2
  double beta2 = p2/(p2+GetMass()*GetMass());
  double theta2t = 14.1*14.1/(beta2*p2*1e6) * (1. + par[3]*par[3]);
  //
  // account for the MS in the layers between the last measurement and the vertex
  for (int il=0;il<nMS;il++) {
    double dfx = xlMS[il] - xCurr;
    xCurr = xlMS[il];
    p2Curr += dfx*cnv*par[4];   // p2 at the scattering layer
    double dxL=xv - xCurr;    // distance from scatering layer to vtx
    double x2rL=cnv*par[4]*dxL, f1L=p2Curr, f2L=f1L+x2rL;
    if (TMath::Abs(f1L) >= kAlmost1 || TMath::Abs(f2L) >= kAlmost1) {
      AliInfo(Form("FailMS at step %d of %d: dfx:%e dxL:%e %e %e",il,nMS,dfx,dxL,f1L,f2L));
      return kFALSE;
    }
    double r1L=TMath::Sqrt((1.-f1L)*(1.+f1L)), r2L=TMath::Sqrt((1.-f2L)*(1.+f2L)), dx2rL=dxL/(r1L+r2L);
    // elements of matrix for propagation from scatering layer to vertex
    double f02L = 2*dx2rL, f04L = cnv*dxL*dx2rL, f13L/*, f24L = cnv*dxL*/;
    if (TMath::Abs(x2rL)<0.05) f13L = dxL*r2L+f2L*(f1L+f2L)*dx2rL; // see AliExternalTrackParam::PropagateTo
    else {
      double dy2dxL = (f1L+f2L)/(r1L+r2L);
      f13L = 2*TMath::ASin(0.5*TMath::Sqrt(1+dy2dxL*dy2dxL)*x2rL)/(cnv*par[4]);
    }
    // MS contribution matrix:
    double theta2 = theta2t*TMath::Abs(x2X0MS[il]);
    double ms22 = theta2*(1.-p2Curr)*(1.+p2Curr)*ms22t;
    double ms33 = theta2*ms33t;
    double ms34 = theta2*ms34t;
    double ms44 = theta2*ms44t;
    //
    // add  H F MS F^Tr H^Tr to cv
    cv00 += f02L*f02L*ms22 + f04L*f04L*ms44;
    cv01 += f04L*f13L*ms34;
    cv11 += f13L*f13L*ms33;
  }
  //
  // inverse of matrix B
  double b11 = ers[1]*ers[1] + cv00;
  double b00 = ers[2]*ers[2] + cv11;
  double det = b11*b00 - cv01*cv01;
  if (TMath::Abs(det)<kAlmost0) {
    AliInfo(Form("Fail on det %e: %e %e %e",det,cv00,cv11,cv01));
    return kFALSE;
  }
  det = 1./det;
  b00 *= det; b11 *= det; 
  double b01 = -cv01*det;
  //
  // elements of matrix DC^Tr * B^-1
  double dcb00 = b00*dc00+b01*dc10, dcb01 = b01*dc00+b11*dc10;
  double dcb10 = b00*dc01+b01*dc11, dcb11 = b01*dc01+b11*dc11;
  double dcb20 = b00*dc02+b01*dc12, dcb21 = b01*dc02+b11*dc12;
  double dcb30 = b00*dc03+b01*dc13, dcb31 = b01*dc03+b11*dc13;
  double dcb40 = b00*dc04+b01*dc14, dcb41 = b01*dc04+b11*dc14;
  //
  // p_{k|k+1} = p_{k|k} +  C_{k|k}*D^Tr_{k+1} B^{-1}_{k+1} ( vtx - D_{k+1}*p_{k|k})
  par[0] += dcb00*yv + dcb01*zv;
  par[1] += dcb10*yv + dcb11*zv;
  par[2] += dcb20*yv + dcb21*zv;
  par[3] += dcb30*yv + dcb31*zv;
  par[4] += dcb40*yv + dcb41*zv;
  //
  // C_{k|kv} = C_{k|k} - C_{k|k} D^Tr_{k+1} B^{-1}_{k+1} D_{k+1} C^Tr_{k|k})
  //
  cov00 = c00 - (dc00*dcb00 + dc10*dcb01);
  cov01 = c01 - (dc01*dcb00 + dc11*dcb01);
  cov02 = c02 - (dc02*dcb00 + dc12*dcb01);
  cov03 = c03 - (dc03*dcb00 + dc13*dcb01);
  cov04 = c04 - (dc04*dcb00 + dc14*dcb01);
  //
  cov11 = c11 - (dc01*dcb10 + dc11*dcb11);
  cov12 = c12 - (dc02*dcb10 + dc12*dcb11);
  cov13 = c13 - (dc03*dcb10 + dc13*dcb11);
  cov14 = c14 - (dc04*dcb10 + dc14*dcb11);
  //
  cov22 = c22 - (dc02*dcb20 + dc12*dcb21);
  cov23 = c23 - (dc03*dcb20 + dc13*dcb21);
  cov24 = c24 - (dc04*dcb20 + dc14*dcb21);
  //
  cov33 = c33 - (dc03*dcb30 + dc13*dcb31);
  cov34 = c34 - (dc04*dcb30 + dc14*dcb31);
  //
  cov44 = c44 - (dc04*dcb40 + dc14*dcb41);
  //
  Set(x,alpha,par,covc);
  if (!Invariant()) {
    AliInfo(Form("Fail on Invariant, X=%e",GetX()));
    return kFALSE;
  }
  return kTRUE;
  //
}
