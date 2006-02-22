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

//-----------------------------------------------------------------
//   Implementation of the derived class for track residuals
//   based on linear chi2 minimization (in approximation of
//   small alignment angles and translations)
//
//-----------------------------------------------------------------

#include <TMinuit.h>
#include <TGeoMatrix.h>

#include "AliAlignObj.h"
#include "AliTrackPointArray.h"
#include "AliTrackResidualsFast.h"

ClassImp(AliTrackResidualsFast)

//______________________________________________________________________________
AliTrackResidualsFast::AliTrackResidualsFast():AliTrackResiduals()
{
  // Default constructor
  for (Int_t i = 0; i < 16; i++) fSum[i] = 0;
  fSumR = 0;
}

//______________________________________________________________________________
AliTrackResidualsFast::AliTrackResidualsFast(Int_t ntracks):
  AliTrackResiduals(ntracks)
{
  // Constructor
  for (Int_t i = 0; i < 16; i++) fSum[i] = 0;
  fSumR = 0;
}
 
//______________________________________________________________________________
AliTrackResidualsFast::AliTrackResidualsFast(const AliTrackResidualsFast &res):
  AliTrackResiduals(res)
{
  // Copy constructor
  for (Int_t i = 0; i < 16; i++) fSum[i] = res.fSum[i];
  fSumR = res.fSumR;
}

//______________________________________________________________________________
AliTrackResidualsFast &AliTrackResidualsFast::operator= (const AliTrackResidualsFast& res)
{
  // Assignment operator
 ((AliTrackResiduals *)this)->operator=(res);
 for (Int_t i = 0; i < 16; i++) fSum[i] = res.fSum[i];
 fSumR = res.fSumR;

 return *this;
}

//______________________________________________________________________________
Bool_t AliTrackResidualsFast::Minimize()
{
  // Implementation of fast linear Chi2
  // based minimization of track residuals sum
  for (Int_t i = 0; i < 16; i++) fSum[i] = 0;
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
  p.GetXYZ(xyz); pprime.GetXYZ(xyzp);
  Double_t weight = 1;
  weight *= weight;
  fSum[0] += weight;
  fSum[1] += xyz[0]*weight;
  fSum[2] += xyz[1]*weight;
  fSum[3] += xyz[2]*weight;
  fSum[4] += (xyz[0]*xyz[0]+xyz[1]*xyz[1])*weight;
  fSum[5] += (xyz[0]*xyz[0]+xyz[2]*xyz[2])*weight;
  fSum[6] += (xyz[1]*xyz[1]+xyz[2]*xyz[2])*weight;
  fSum[7] -= xyz[0]*xyz[1]*weight;
  fSum[8] -= xyz[0]*xyz[2]*weight;
  fSum[9] -= xyz[1]*xyz[2]*weight;
  fSum[10] += (xyzp[0]-xyz[0])*weight;
  fSum[11] += (xyzp[1]-xyz[1])*weight;
  fSum[12] += (xyzp[2]-xyz[2])*weight;
  fSum[13] += (xyzp[2]*xyz[1]-xyzp[1]*xyz[2])*weight;
  fSum[14] += (xyzp[0]*xyz[2]-xyzp[2]*xyz[0])*weight;
  fSum[15] += (xyzp[1]*xyz[0]-xyzp[0]*xyz[1])*weight;

  fSumR += ((xyzp[0]-xyz[0])*(xyzp[0]-xyz[0])+
	    (xyzp[1]-xyz[1])*(xyzp[1]-xyz[1])+
	    (xyzp[2]-xyz[2])*(xyzp[2]-xyz[2]))*weight;

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

  smatrix(0,0) = fSum[0]; smatrix(1,1) = fSum[0]; smatrix(2,2)=fSum[0];
  smatrix(3,3) = fSum[6]; smatrix(4,4) = fSum[5]; smatrix(5,5)=fSum[4];

  smatrix(0,1) = 0; smatrix(0,2) = 0;
  smatrix(0,3) = 0; smatrix(0,4) = fSum[3]; smatrix(0,5) = -fSum[2];
  smatrix(1,2) = 0;
  smatrix(1,3) = -fSum[3]; smatrix(1,4) = 0; smatrix(1,5) = fSum[1];
  smatrix(2,3) = fSum[2]; smatrix(2,4) = -fSum[1]; smatrix(2,5) = 0;

  smatrix(3,4) = fSum[7]; smatrix(3,5) = fSum[8];
  smatrix(4,5) = fSum[9];

  sums(0,0)    = fSum[10]; sums(0,1) = fSum[11]; sums(0,2) = fSum[12];
  sums(0,3)    = fSum[13]; sums(0,4) = fSum[14]; sums(0,5) = fSum[15];

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
