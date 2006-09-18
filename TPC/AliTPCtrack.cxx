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

//-----------------------------------------------------------------
//           Implementation of the TPC track class
//        This class is used by the AliTPCtracker class
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include <Riostream.h>

#include "AliTPCtrack.h"
#include "AliCluster.h"
#include "AliTracker.h"
#include "AliESDtrack.h"

ClassImp(AliTPCtrack)

//_________________________________________________________________________
AliTPCtrack::AliTPCtrack(): 
  AliKalmanTrack(),
  fdEdx(0),
  fSdEdx(1e10),
  fNFoundable(0),
  fBConstrain(kFALSE),
  fLastPoint(-1),
  fFirstPoint(-1),
  fRemoval(0),
  fTrackType(0),
  fLab2(-1),
  fNShared(0),
  fReference()
{
  //-------------------------------------------------
  // default constructor
  //-------------------------------------------------
  for (Int_t i=0; i<kMaxRow;i++) fIndex[i]=-2;
  for (Int_t i=0; i<4;i++) fPoints[i]=0.;
  for (Int_t i=0; i<12;i++) fKinkPoint[i]=0.;
  for (Int_t i=0; i<3;i++) fKinkIndexes[i]=0;
  for (Int_t i=0; i<3;i++) fV0Indexes[i]=0;
}

//_________________________________________________________________________



AliTPCtrack::AliTPCtrack(Double_t x, Double_t alpha, const Double_t p[5],
		         const Double_t cov[15], Int_t index) :
  AliKalmanTrack(),
  fdEdx(0),
  fSdEdx(1e10),
  fNFoundable(0),
  fBConstrain(kFALSE),
  fLastPoint(0),
  fFirstPoint(0),
  fRemoval(0),
  fTrackType(0),
  fLab2(0),
  fNShared(0),
  fReference()
{
  //-----------------------------------------------------------------
  // This is the main track constructor.
  //-----------------------------------------------------------------
  Double_t cnv=1./(AliTracker::GetBz()*kB2C);

  Double_t pp[5]={
    p[0],
    p[1],
    x*p[4] - p[2],
    p[3],
    p[4]*cnv
  };

  Double_t c22 = x*x*cov[14] - 2*x*cov[12] + cov[5];
  Double_t c32 = x*cov[13] - cov[8];
  Double_t c20 = x*cov[10] - cov[3], 
           c21 = x*cov[11] - cov[4], c42 = x*cov[14] - cov[12];

  Double_t cc[15]={
    cov[0 ],
    cov[1 ],     cov[2 ],
    c20,         c21,         c22,
    cov[6 ],     cov[7 ],     c32,     cov[9 ],
    cov[10]*cnv, cov[11]*cnv, c42*cnv, cov[13]*cnv, cov[14]*cnv*cnv
  };

  Double_t p0=TMath::Sign(1/kMostProbablePt,pp[4]);
  Double_t w0=cc[14]/(cc[14] + p0*p0), w1=p0*p0/(cc[14] + p0*p0);
  pp[4] = w0*p0 + w1*pp[4]; 
  cc[10]*=w1; cc[11]*=w1; cc[12]*=w1; cc[13]*=w1; cc[14]*=w1;

  Set(x,alpha,pp,cc);

  SetNumberOfClusters(1);
  
  fIndex[0]=index;
  for (Int_t i=1; i<kMaxRow;i++) fIndex[i]=-2;
  for (Int_t i=0; i<4;i++) fPoints[i]=0.;
  for (Int_t i=0; i<12;i++) fKinkPoint[i]=0.;
  for (Int_t i=0; i<3;i++) fKinkIndexes[i]=0;
  for (Int_t i=0; i<3;i++) fV0Indexes[i]=0;
}

//_____________________________________________________________________________
AliTPCtrack::AliTPCtrack(const AliESDtrack& t) :
  AliKalmanTrack(),
  fdEdx(t.GetTPCsignal()),
  fSdEdx(1e10),
  fNFoundable(0),
  fBConstrain(kFALSE),
  fLastPoint(0),
  fFirstPoint(0),
  fRemoval(0),
  fTrackType(0),
  fLab2(0),
  fNShared(0),
  fReference()
{
  //-----------------------------------------------------------------
  // Conversion AliESDtrack -> AliTPCtrack.
  //-----------------------------------------------------------------
  SetNumberOfClusters(t.GetTPCclusters(fIndex));
  SetLabel(t.GetLabel());
  SetMass(t.GetMass());
  for (Int_t i=0; i<4;i++) fPoints[i]=0.;
  for (Int_t i=0; i<12;i++) fKinkPoint[i]=0.;
  for (Int_t i=0; i<3;i++) fKinkIndexes[i]=0;
  for (Int_t i=0; i<3;i++) fV0Indexes[i]=0;

  Set(t.GetX(),t.GetAlpha(),t.GetParameter(),t.GetCovariance());

  if ((t.GetStatus()&AliESDtrack::kTIME) == 0) return;
  StartTimeIntegral();
  Double_t times[10]; t.GetIntegratedTimes(times); SetIntegratedTimes(times);
  SetIntegratedLength(t.GetIntegratedLength());
}

//_____________________________________________________________________________
AliTPCtrack::AliTPCtrack(const AliTPCtrack& t) :
  AliKalmanTrack(t),
  fdEdx(t.fdEdx),
  fSdEdx(t.fSdEdx),
  fNFoundable(t.fNFoundable),
  fBConstrain(t.fBConstrain),
  fLastPoint(t.fLastPoint),
  fFirstPoint(t.fFirstPoint),
  fRemoval(t.fRemoval),
  fTrackType(t.fTrackType),
  fLab2(t.fLab2),
  fNShared(t.fNShared),
  fReference(t.fReference)

{
  //-----------------------------------------------------------------
  // This is a track copy constructor.
  //-----------------------------------------------------------------
  Set(t.GetX(),t.GetAlpha(),t.GetParameter(),t.GetCovariance());

  for (Int_t i=0; i<kMaxRow; i++) fIndex[i]=t.fIndex[i];
  for (Int_t i=0; i<4;i++) fPoints[i]=t.fPoints[i];
  for (Int_t i=0; i<12;i++) fKinkPoint[i]=t.fKinkPoint[i];
  for (Int_t i=0; i<3;i++) fKinkIndexes[i]=t.fKinkIndexes[i];
  for (Int_t i=0; i<3;i++) fV0Indexes[i]=t.fV0Indexes[i];
}

//_____________________________________________________________________________
Int_t AliTPCtrack::Compare(const TObject *o) const {
  //-----------------------------------------------------------------
  // This function compares tracks according to the their curvature
  //-----------------------------------------------------------------
  AliTPCtrack *t=(AliTPCtrack*)o;
  //Double_t co=TMath::Abs(t->Get1Pt());
  //Double_t c =TMath::Abs(Get1Pt());
  Double_t co=t->GetSigmaY2()*t->GetSigmaZ2();
  Double_t c =GetSigmaY2()*GetSigmaZ2();
  if (c>co) return 1;
  else if (c<co) return -1;
  return 0;
}

Double_t AliTPCtrack::GetPredictedChi2(const AliCluster *c) const {
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  Double_t p[2]={c->GetY(), c->GetZ()};
  Double_t cov[3]={c->GetSigmaY2(), 0., c->GetSigmaZ2()};
  return AliExternalTrackParam::GetPredictedChi2(p,cov);
}

//_____________________________________________________________________________
Bool_t AliTPCtrack::PropagateTo(Double_t xk,Double_t rho,Double_t x0) {
  //-----------------------------------------------------------------
  //  This function propagates a track to a reference plane x=xk.
  //  rho - density of the crossed matrial (g/cm^3)
  //  x0  - radiation length of the crossed material (g/cm^2) 
  //-----------------------------------------------------------------
  Double_t oldX=GetX(), oldY=GetY(), oldZ=GetZ();

  Double_t bz=GetBz();
  if (!AliExternalTrackParam::PropagateTo(xk,bz)) return kFALSE;

  Double_t d = TMath::Sqrt((GetX()-oldX)*(GetX()-oldX) + 
                           (GetY()-oldY)*(GetY()-oldY) + 
                           (GetZ()-oldZ)*(GetZ()-oldZ));
  if (IsStartedTimeIntegral() && GetX()>oldX) AddTimeStep(d);

  if (oldX < xk) d = -d;
  if (!AliExternalTrackParam::CorrectForMaterial(d*rho/x0,x0,GetMass())) 
     return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t 
AliTPCtrack::PropagateToVertex(const AliESDVertex *v,Double_t d,Double_t x0) 
{
  //-----------------------------------------------------------------
  // This function propagates tracks to the vertex.
  //-----------------------------------------------------------------
  Double_t bz=GetBz();
  if (PropagateToDCA(v,bz,kVeryBig))
   if (AliExternalTrackParam::CorrectForMaterial(d,x0,GetMass())) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliTPCtrack::Update(const AliCluster *c, Double_t chisq, Int_t index) {
  //-----------------------------------------------------------------
  // This function associates a cluster with this track.
  //-----------------------------------------------------------------
  Double_t p[2]={c->GetY(), c->GetZ()};
  Double_t cov[3]={c->GetSigmaY2(), 0., c->GetSigmaZ2()};

  if (!AliExternalTrackParam::Update(p,cov)) return kFALSE;

  Int_t n=GetNumberOfClusters();
  fIndex[n]=index;
  SetNumberOfClusters(n+1);
  SetChi2(GetChi2()+chisq);

  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
// MI ADDITION

Float_t AliTPCtrack::Density(Int_t row0, Int_t row1)
{
  //
  // calculate cluster density
  Int_t good  = 0;
  Int_t found = 0;
  //if (row0<fFirstPoint) row0 = fFirstPoint;
  if (row1>fLastPoint) row1 = fLastPoint;

  
  for (Int_t i=row0;i<=row1;i++){ 
    //    Int_t index = fClusterIndex[i];
    Int_t index = fIndex[i];
    if (index!=-1)  good++;
    if (index>0)    found++;
  }
  Float_t density=0;
  if (good>0) density = Float_t(found)/Float_t(good);
  return density;
}


Float_t AliTPCtrack::Density2(Int_t row0, Int_t row1)
{
  //
  // calculate cluster density
  Int_t good  = 0;
  Int_t found = 0;
  //  
  for (Int_t i=row0;i<=row1;i++){     
    Int_t index = fIndex[i];
    if (index!=-1)  good++;
    if (index>0)    found++;
  }
  Float_t density=0;
  if (good>0) density = Float_t(found)/Float_t(good);
  return density;
}

void  AliTPCtrack::UpdatePoints()
{
  //--------------------------------------------------
  //calculates first ,amx dens and last points
  //--------------------------------------------------
  Float_t density[160];
  for (Int_t i=0;i<160;i++) density[i]=-1.;
  fPoints[0]= 160;
  fPoints[1] = -1;
  //
  Int_t ngood=0;
  Int_t undeff=0;
  Int_t nall =0;
  Int_t range=20;
  for (Int_t i=0;i<160;i++){
    Int_t last = i-range;
    if (nall<range) nall++;
    if (last>=0){
      if (fIndex[last]>0&& (fIndex[last]&0x8000)==0) ngood--;
      if (fIndex[last]==-1) undeff--;
    }
    if (fIndex[i]>0&& (fIndex[i]&0x8000)==0)   ngood++;
    if (fIndex[i]==-1) undeff++;
    if (nall==range &&undeff<range/2) density[i-range/2] = Float_t(ngood)/Float_t(nall-undeff);
  }
  Float_t maxdens=0;
  Int_t indexmax =0;
  for (Int_t i=0;i<160;i++){
    if (density[i]<0) continue;
    if (density[i]>maxdens){
      maxdens=density[i];
      indexmax=i;
    }
  }
  //
  //max dens point
  fPoints[3] = maxdens;
  fPoints[1] = indexmax;
  //
  // last point
  for (Int_t i=indexmax;i<160;i++){
    if (density[i]<0) continue;
    if (density[i]<maxdens/2.) {
      break;
    }
    fPoints[2]=i;
  }
  //
  // first point
  for (Int_t i=indexmax;i>0;i--){
    if (density[i]<0) continue;
    if (density[i]<maxdens/2.) {
      break;
    }
    fPoints[0]=i;
  }
  //
}

Double_t AliTPCtrack::GetBz() const {
  //
  // returns Bz component of the magnetic field (kG)
  //
  if (AliTracker::UniformField()) return AliTracker::GetBz();
  Double_t r[3]; GetXYZ(r);
  return AliTracker::GetBz(r);
}

