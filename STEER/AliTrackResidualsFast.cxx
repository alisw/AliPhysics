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
//   Implementation of the derived class for track residuals
//   based on linear chi2 minimization (in approximation of
//   small alignment angles and translations)
//
//-----------------------------------------------------------------

#include <TMath.h>
#include <TMinuit.h>
#include <TGeoMatrix.h>

#include "AliLog.h"
#include "AliAlignObj.h"
#include "AliTrackPointArray.h"
#include "AliTrackResidualsFast.h"

ClassImp(AliTrackResidualsFast)

//______________________________________________________________________________
AliTrackResidualsFast::AliTrackResidualsFast():
  AliTrackResiduals(),
  fSumR(0)
{
  // Default constructor
  for (Int_t i = 0; i < 27; i++) fSum[i] = 0;
}

//______________________________________________________________________________
AliTrackResidualsFast::AliTrackResidualsFast(Int_t ntracks):
  AliTrackResiduals(ntracks),
  fSumR(0)
{
  // Constructor
  for (Int_t i = 0; i < 27; i++) fSum[i] = 0;
}
 
//______________________________________________________________________________
AliTrackResidualsFast::AliTrackResidualsFast(const AliTrackResidualsFast &res):
  AliTrackResiduals(res),
  fSumR(res.fSumR)
{
  // Copy constructor
  for (Int_t i = 0; i < 27; i++) fSum[i] = res.fSum[i];
}

//______________________________________________________________________________
AliTrackResidualsFast &AliTrackResidualsFast::operator= (const AliTrackResidualsFast& res)
{
  // Assignment operator
 ((AliTrackResiduals *)this)->operator=(res);
 for (Int_t i = 0; i < 27; i++) fSum[i] = res.fSum[i];
 fSumR = res.fSumR;

 return *this;
}

//______________________________________________________________________________
Bool_t AliTrackResidualsFast::Minimize()
{
  // Implementation of fast linear Chi2
  // based minimization of track residuals sum
  for (Int_t i = 0; i < 27; i++) fSum[i] = 0;
  fSumR = 0;

  AliTrackPoint p1,p2;

  for (Int_t itrack = 0; itrack < fLast; itrack++) {
    if (!fVolArray[itrack] || !fTrackArray[itrack]) continue;
    for (Int_t ipoint = 0; ipoint < fVolArray[itrack]->GetNPoints(); ipoint++) {
      fVolArray[itrack]->GetPoint(p1,ipoint);
      fTrackArray[itrack]->GetPoint(p2,ipoint);
      AddPoints(p1,p2);
    }
  }

  return Update();

  // debug info
//   Float_t chi2 = 0;
//   for (Int_t itrack = 0; itrack < fLast; itrack++) {
//     if (!fVolArray[itrack] || !fTrackArray[itrack]) continue;
//     for (Int_t ipoint = 0; ipoint < fVolArray[itrack]->GetNPoints(); ipoint++) {
//       fVolArray[itrack]->GetPoint(p1,ipoint);
//       fAlignObj->Transform(p1);
//       fTrackArray[itrack]->GetPoint(p2,ipoint);
//       Float_t residual = p2.GetResidual(p1,kFALSE);
//       chi2 += residual;
//     }
//   }
//   printf("Final chi2 = %f\n",chi2);
}

//______________________________________________________________________________
void AliTrackResidualsFast::AddPoints(AliTrackPoint &p, AliTrackPoint &pprime)
{
  // Update the sums used for
  // the linear chi2 minimization
  Float_t xyz[3],xyzp[3];
  Float_t cov[6],covp[6];
  p.GetXYZ(xyz,cov); pprime.GetXYZ(xyzp,covp);
  TMatrixDSym mcov(3);
  mcov(0,0) = cov[0]; mcov(0,1) = cov[1]; mcov(0,2) = cov[2];
  mcov(1,0) = cov[1]; mcov(1,1) = cov[3]; mcov(1,2) = cov[4];
  mcov(2,0) = cov[2]; mcov(2,1) = cov[4]; mcov(2,2) = cov[5];
  TMatrixDSym mcovp(3);
  mcovp(0,0) = covp[0]; mcovp(0,1) = covp[1]; mcovp(0,2) = covp[2];
  mcovp(1,0) = covp[1]; mcovp(1,1) = covp[3]; mcovp(1,2) = covp[4];
  mcovp(2,0) = covp[2]; mcovp(2,1) = covp[4]; mcovp(2,2) = covp[5];
  TMatrixDSym msum = mcov + mcovp;
  msum.Invert();
  if (!msum.IsValid()) return;

  TMatrixD        sums(3,1);
  sums(0,0) = (xyzp[0]-xyz[0]); 
  sums(1,0) = (xyzp[1]-xyz[1]);
  sums(2,0) = (xyzp[2]-xyz[2]);   
  TMatrixD sumst = sums.T(); sums.T();

  TMatrixD        mf(3,6);
  mf(0,0) = 1;      mf(1,0) = 0;       mf(2,0) = 0;
  mf(0,1) = 0;      mf(1,1) = 1;       mf(2,1) = 0;
  mf(0,2) = 0;      mf(1,2) = 0;       mf(2,2) = 1;
  mf(0,3) = 0;      mf(1,3) =-xyz[2];  mf(2,3) = xyz[1];
  mf(0,4) = xyz[2]; mf(1,4) = 0;       mf(2,4) =-xyz[0];
  mf(0,5) =-xyz[1]; mf(1,5) = xyz[0];  mf(2,5) = 0;
  TMatrixD        mft = mf.T(); mf.T();
  TMatrixD sums2 = mft * msum * sums;

  TMatrixD smatrix = mft * msum * mf;

  fSum[0] += smatrix(0,0);
  fSum[1] += smatrix(0,1);
  fSum[2] += smatrix(0,2);
  fSum[3] += smatrix(0,3);
  fSum[4] += smatrix(0,4);
  fSum[5] += smatrix(0,5);
  fSum[6] += smatrix(1,1);
  fSum[7] += smatrix(1,2);
  fSum[8] += smatrix(1,3);
  fSum[9] += smatrix(1,4);
  fSum[10]+= smatrix(1,5);
  fSum[11]+= smatrix(2,2);
  fSum[12]+= smatrix(2,3);
  fSum[13]+= smatrix(2,4);
  fSum[14]+= smatrix(2,5);
  fSum[15]+= smatrix(3,3);
  fSum[16]+= smatrix(3,4);
  fSum[17]+= smatrix(3,5);
  fSum[18]+= smatrix(4,4);
  fSum[19]+= smatrix(4,5);
  fSum[20]+= smatrix(5,5);
  fSum[21] += sums2(0,0);
  fSum[22] += sums2(1,0);
  fSum[23] += sums2(2,0);
  fSum[24] += sums2(3,0);
  fSum[25] += sums2(4,0);
  fSum[26] += sums2(5,0);

  TMatrixD  tmp = sumst * msum * sums;
  fSumR += tmp(0,0);

  fNdf += 3;
}

//______________________________________________________________________________
Bool_t AliTrackResidualsFast::Update()
{
  // Find the alignment parameters
  // by using the already accumulated
  // sums
  TMatrixDSym     smatrix(6);
  TMatrixD        sums(1,6);

  smatrix(0,0) = fSum[0];
  smatrix(0,1) = smatrix(1,0) = fSum[1];
  smatrix(0,2) = smatrix(2,0) = fSum[2];
  smatrix(0,3) = smatrix(3,0) = fSum[3];
  smatrix(0,4) = smatrix(4,0) = fSum[4];
  smatrix(0,5) = smatrix(5,0) = fSum[5];
  smatrix(1,1) = fSum[6];
  smatrix(1,2) = smatrix(2,1) = fSum[7];
  smatrix(1,3) = smatrix(3,1) = fSum[8];
  smatrix(1,4) = smatrix(4,1) = fSum[9];
  smatrix(1,5) = smatrix(5,1) = fSum[10];
  smatrix(2,2) = fSum[11];
  smatrix(2,3) = smatrix(3,2) = fSum[12];
  smatrix(2,4) = smatrix(4,2) = fSum[13];
  smatrix(2,5) = smatrix(5,2) = fSum[14];
  smatrix(3,3) = fSum[15];
  smatrix(3,4) = smatrix(4,3) = fSum[16];
  smatrix(3,5) = smatrix(5,3) = fSum[17];
  smatrix(4,4) = fSum[18];
  smatrix(4,5) = smatrix(5,4) = fSum[19];
  smatrix(5,5) = fSum[20];

  sums(0,0)    = fSum[21]; sums(0,1) = fSum[22]; sums(0,2) = fSum[23];
  sums(0,3)    = fSum[24]; sums(0,4) = fSum[25]; sums(0,5) = fSum[26];

  smatrix.Invert();
  if (!smatrix.IsValid()) return kFALSE;

  TMatrixD res = sums*smatrix;
  fAlignObj->SetPars(res(0,0),res(0,1),res(0,2),
		     TMath::RadToDeg()*res(0,3),
		     TMath::RadToDeg()*res(0,4),
		     TMath::RadToDeg()*res(0,5));
  TMatrixD  tmp = res*sums.T();
  fChi2 = fSumR - tmp(0,0);
  fNdf -= 6;

  return kTRUE;
}
