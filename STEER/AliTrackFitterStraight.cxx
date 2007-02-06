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

#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMatrixD.h>

#include "AliTrackFitterStraight.h"

ClassImp(AliTrackFitterStraight)

AliTrackFitterStraight::AliTrackFitterStraight():
  AliTrackFitter(),
  fAlpha(0.),
  fSumYY(0),
  fSumZZ(0),
  fNUsed(0),
  fConv(kFALSE)
{
  //
  // default constructor
  //
  for (Int_t i=0;i<5;i++) fSumXY[i] = 0;
  for (Int_t i=0;i<5;i++) fSumXZ[i] = 0;
}


AliTrackFitterStraight::AliTrackFitterStraight(AliTrackPointArray *array, Bool_t owner):
  AliTrackFitter(array,owner),
  fAlpha(0.),
  fSumYY(0),
  fSumZZ(0),
  fNUsed(0),
  fConv(kFALSE)
{
  //
  // Constructor
  //
  for (Int_t i=0;i<5;i++) fSumXY[i] = 0;
  for (Int_t i=0;i<5;i++) fSumXZ[i] = 0;
}

AliTrackFitterStraight::AliTrackFitterStraight(const AliTrackFitterStraight &fitter):
  AliTrackFitter(fitter),
  fAlpha(fitter.fAlpha),
  fSumYY(fitter.fSumYY),
  fSumZZ(fitter.fSumZZ),
  fNUsed(fitter.fNUsed),
  fConv(fitter.fConv)
{
  //
  // copy constructor
  //
  for (Int_t i=0;i<5;i++) fSumXY[i]  = fitter.fSumXY[i];
  for (Int_t i=0;i<5;i++) fSumXZ[i]  = fitter.fSumXZ[i];
}

//_____________________________________________________________________________
AliTrackFitterStraight &AliTrackFitterStraight::operator =(const AliTrackFitterStraight& fitter)
{
  // assignment operator
  //
  if(this==&fitter) return *this;
  ((AliTrackFitter *)this)->operator=(fitter);

  fAlpha = fitter.fAlpha;
  for (Int_t i=0;i<5;i++) fSumXY[i]  = fitter.fSumXY[i];
  fSumYY = fitter.fSumYY;
  for (Int_t i=0;i<5;i++) fSumXZ[i]  = fitter.fSumXZ[i];
  fSumZZ = fitter.fSumZZ;
  fNUsed = fitter.fNUsed;
  fConv = fitter.fConv;

  return *this;
}

AliTrackFitterStraight::~AliTrackFitterStraight()
{
  // destructor
  //
}

void AliTrackFitterStraight::Reset()
{
  // Reset the track parameters and
  // sums
  AliTrackFitter::Reset();
  fAlpha = 0.;
  for (Int_t i=0;i<5;i++) fSumXY[i] = 0;
  fSumYY = 0;
  for (Int_t i=0;i<5;i++) fSumXZ[i] = 0;
  fSumZZ = 0;
  fNUsed = 0;
  fConv =kFALSE;
}

Bool_t AliTrackFitterStraight::Fit(const TArrayI *volIds,const TArrayI *volIdsFit,
				 AliAlignObj::ELayerID layerRangeMin,
				 AliAlignObj::ELayerID layerRangeMax)
{
  // Fit the track points. The method takes as an input
  // the set of id's (volids) of the volumes in which
  // one wants to calculate the residuals.
  // The following parameters are used to define the
  // range of volumes to be used in the fitting
  // As a result two AliTrackPointArray's obects are filled.
  // The first one contains the space points with
  // volume id's from volids list. The second array of points represents
  // the track extrapolations corresponding to the space points
  // in the first array. The two arrays can be used to find
  // the residuals in the volids and consequently construct a
  // chi2 function to be minimized during the alignment
  // procedures. For the moment the track extrapolation is taken
  // at the space-point reference plane. The reference plane is
  // found using the covariance matrix of the point
  // (assuming sigma(x)=0 at the reference coordinate system.

  Reset();

  Int_t npoints = fPoints->GetNPoints();
  if (npoints < 2) return kFALSE;

  Bool_t isAlphaCalc = kFALSE;
  AliTrackPoint p,plocal;

  Int_t npVolId = 0;
  fNUsed = 0;
  Int_t *pindex = new Int_t[npoints];
  for (Int_t ipoint = 0; ipoint < npoints; ipoint++)
    {
      fPoints->GetPoint(p,ipoint);
      UShort_t iVolId = p.GetVolumeID();
      if (FindVolId(volIds,iVolId)) {
	pindex[npVolId] = ipoint;
	npVolId++;
      }
      if (volIdsFit != 0x0) {
	if (!FindVolId(volIdsFit,iVolId)) continue;
      }
      else {
	if (iVolId < AliAlignObj::LayerToVolUID(layerRangeMin,0) ||
	    iVolId > AliAlignObj::LayerToVolUID(layerRangeMax,
						AliAlignObj::LayerSize(layerRangeMax))) continue;
      }
      if (!isAlphaCalc) {
	fAlpha = p.GetAngle();
	isAlphaCalc = kTRUE;
      }
      plocal = p.Rotate(fAlpha);
      AddPoint(plocal.GetX(),plocal.GetY(),plocal.GetZ(),
	       TMath::Sqrt(plocal.GetCov()[3]),TMath::Sqrt(plocal.GetCov()[5]));
      fNUsed++;
    }

  if (fNUsed < 2) {
    delete [] pindex;
    return kFALSE;
  }

  Update();

  if (!fConv) {
    delete [] pindex;
    return kFALSE;
  }

  if (fNUsed < fMinNPoints ) {
    delete [] pindex;
    return kFALSE;
  }

  fPVolId = new AliTrackPointArray(npVolId);
  fPTrack = new AliTrackPointArray(npVolId);
  AliTrackPoint p2;
  for (Int_t ipoint = 0; ipoint < npVolId; ipoint++)
    {
      Int_t index = pindex[ipoint];
      fPoints->GetPoint(p,index);
      if (GetPCA(p,p2)) {
	Float_t xyz[3],xyz2[3];
	p.GetXYZ(xyz); p2.GetXYZ(xyz2);
	//	printf("residuals %f %d %d %f %f %f %f %f %f\n",fChi2,fNUsed,fConv,xyz[0],xyz[1],xyz[2],xyz2[0]-xyz[0],xyz2[1]-xyz[1],xyz2[2]-xyz[2]);
	fPVolId->AddPoint(ipoint,&p);
	fPTrack->AddPoint(ipoint,&p2);
      }
    }  

  delete [] pindex;

  return kTRUE;
}

void AliTrackFitterStraight::AddPoint(Float_t x, Float_t y, Float_t z, Float_t sy, Float_t sz)
{
  // Straight track fitter
  // The method add a point to the sums
  // used to extract track parameters

  //
  // XY part
  //
  Double_t weight = 1./(sy*sy);
  fSumXY[0] +=weight;
  fSumXY[1] +=x*weight;      fSumXY[2] +=x*x*weight;
  fSumXY[3] +=y*weight;      fSumXY[4] +=x*y*weight;
  fSumYY += y*y*weight;
  //
  // XZ part
  //
  weight = 1./(sz*sz);
  fSumXZ[0] +=weight;
  fSumXZ[1] +=x*weight;      fSumXZ[2] +=x*x*weight;
  fSumXZ[3] +=z*weight;      fSumXZ[4] +=x*z*weight;
  fSumZZ += z*z*weight;
}

void AliTrackFitterStraight::Update(){
  //
  //  Track fitter update
  //
  //
  for (Int_t i=0;i<6;i++)fParams[i]=0;
  fChi2 = 0;
  fNdf = 0;
  Int_t conv=0;
  //
  // XY part
  //
  TMatrixDSym     smatrix(2);
  TMatrixD        sums(1,2);
  smatrix(0,0) = fSumXY[0]; smatrix(1,1)=fSumXY[2];
  smatrix(0,1) = fSumXY[1]; smatrix(1,0)=fSumXY[1];
  sums(0,0)    = fSumXY[3]; sums(0,1)   =fSumXY[4];
  smatrix.Invert();
  if (smatrix.IsValid()){
    for (Int_t i=0;i<2;i++)
      for (Int_t j=0;j<=i;j++){
	(*fCov)(i,j)=smatrix(i,j);
      }
    TMatrixD  res = sums*smatrix;
    fParams[0] = res(0,0);
    fParams[1] = res(0,1);
    TMatrixD  tmp = res*sums.T();
    fChi2 += fSumYY - tmp(0,0);
    fNdf  += fNUsed - 2;
    conv++;
  }
  //
  // XZ part
  //
  TMatrixDSym     smatrixz(2);
  TMatrixD        sumsxz(1,2);
  smatrixz(0,0) = fSumXZ[0]; smatrixz(1,1) = fSumXZ[2];
  smatrixz(0,1) = fSumXZ[1]; smatrixz(1,0) = fSumXZ[1];
  sumsxz(0,0)   = fSumXZ[3]; sumsxz(0,1)   = fSumXZ[4];
  smatrixz.Invert();
  if (smatrixz.IsValid()){
    TMatrixD res = sumsxz*smatrixz;
    fParams[2] = res(0,0);
    fParams[3] = res(0,1);
    fParams[4] = fParams[5] = 0;
    for (Int_t i=0;i<2;i++)
      for (Int_t j=0;j<=i;j++){
	(*fCov)(i+2,j+2)=smatrixz(i,j);
      }
    TMatrixD  tmp = res*sumsxz.T();
    fChi2 += fSumZZ - tmp(0,0);
    fNdf  += fNUsed - 2;
    conv++;
  }

  if (conv>1)
    fConv =kTRUE;
  else
    fConv=kFALSE;
}

Double_t AliTrackFitterStraight::GetYat(Double_t x) const {
  if (!fConv) return 0.;
  return (fParams[0]+x*fParams[1]);
}

Double_t AliTrackFitterStraight::GetDYat(Double_t x) const {
  if (!fConv) return 0.;
  return fParams[1]+0.*x;
}



Double_t AliTrackFitterStraight::GetZat(Double_t x) const {
  if (!fConv) return 0.;
  return (fParams[2]+x*fParams[3]);
}

Double_t AliTrackFitterStraight::GetDZat(Double_t x) const {
  if (!fConv) return 0.;
  return fParams[3]+0.*x;
}

Bool_t AliTrackFitterStraight::GetXYZat(Double_t r, Float_t *xyz) const {
  if (!fConv) return kFALSE;
  Double_t y = (fParams[0]+r*fParams[1]);
  Double_t z = (fParams[2]+r*fParams[3]);

  Double_t sin = TMath::Sin(fAlpha);
  Double_t cos = TMath::Cos(fAlpha);
  xyz[0] = r*cos - y*sin;
  xyz[1] = y*cos + r*sin;
  xyz[2] = z;

  return kTRUE;
}

Bool_t AliTrackFitterStraight::GetPCA(const AliTrackPoint &p, AliTrackPoint &p2) const
{
  // Get the closest to a given spacepoint track trajectory point
  // Look for details in the description of the Fit() method

  if (!fConv) return kFALSE;

  // First X and Y coordinates
  Double_t sin = TMath::Sin(fAlpha);
  Double_t cos = TMath::Cos(fAlpha);
  // Track parameters in the global coordinate system
  Double_t x0 = -fParams[0]*sin;
  Double_t y0 =  fParams[0]*cos;
  if ((cos - fParams[1]*sin) == 0) return kFALSE;
  Double_t dydx = (fParams[1]*cos + sin)/(cos - fParams[1]*sin);

  // Define space-point refence plane
  Double_t alphap = p.GetAngle();
  Double_t sinp = TMath::Sin(alphap);
  Double_t cosp = TMath::Cos(alphap);
  Double_t x  = p.GetX()*cosp + p.GetY()*sinp;
  //  Double_t y  = p.GetY()*cosp - p.GetX()*sinp;
  Double_t x0p= x0*cosp + y0*sinp;
  Double_t y0p= y0*cosp - x0*sinp;
  if ((cos + dydx*sin) == 0) return kFALSE;
  Double_t dydxp = (dydx*cos - sin)/(cos + dydx*sin);
  Double_t yprime = y0p + dydxp*(x-x0p);

  // Back to the global coordinate system
  Double_t xsecond = x*cosp - yprime*sinp;
  Double_t ysecond = yprime*cosp + x*sinp;

  // Now Z coordinate and track angles
  Double_t x2 = xsecond*cos + ysecond*sin;
  Double_t zsecond = GetZat(x2);
  Double_t dydx2 = fParams[1];
  Double_t dzdx = fParams[3];

  // Fill the cov matrix of the track extrapolation point
  Double_t cov[6] = {0,0,0,0,0,0};
  Double_t sigmax = 100*100.;
  cov[0] = sigmax;             cov[1] = sigmax*dydx2;      cov[2] = sigmax*dzdx;
  cov[3] = sigmax*dydx2*dydx2; cov[4] = sigmax*dydx2*dzdx;
  cov[5] = sigmax*dzdx*dzdx;

  Float_t  newcov[6];
  newcov[0] = cov[0]*cos*cos-
            2*cov[1]*sin*cos+
              cov[3]*sin*sin;
  newcov[1] = cov[1]*(cos*cos-sin*sin)-
             (cov[3]-cov[0])*sin*cos;
  newcov[2] = cov[2]*cos-
              cov[4]*sin;
  newcov[3] = cov[0]*sin*sin+
            2*cov[1]*sin*cos+
              cov[3]*cos*cos;
  newcov[4] = cov[4]*cos+
              cov[2]*sin;
  newcov[5] = cov[5];

  p2.SetXYZ(xsecond,ysecond,zsecond,newcov);

  return kTRUE;
}
