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


////////////////////////////////////////////////////////////////////////////////
//
// AliBarrelTrack is a snapshot of the TPC or TRD track.
// The barrel track is stored by tracker every time a track is crossing 
// reference plane.
//  
// Barrel track have
// - state vactor in "external represenantion"
// - diagonal elements of covariance matric
// - auxiliary paramrters: 
//   chi2, number of clusters, number of sector srossed
// - methods to compare with track references
//
// Barrel track can be directly compared with AliTrackReferences,
// TPC and TRD tracks can be compared using the same macro.
//
// S. Radomski, [GSI] mail: S.Radomski@gsi.de
// 07.04.2003
// 
////////////////////////////////////////////////////////////////////////////////

#include "AliBarrelTrack.h"

ClassImp(AliBarrelTrack)

////////////////////////////////////////////////////////////////////////////////

AliBarrelTrack:: AliBarrelTrack() : TObject() {
  //
  // Standard constructor 
  // reset all data-members
  //

  fLabel = fRefPlane = fIsIn = 0;
  fX =  fAlpha = 0.0;
  fZ = fY = fTgLambda = fSnPhi = f1Pt = 0.0;

  for(Int_t i=0; i<5; i++)  fTimeHypothesis[i] = 0;
 
  fLength = 0.0;
  fNClusters = fNWrong = fNRotate = 0;
  fChi2 = 0.0; 
  fMass = fdEdX = 0;

  fCy = fCz = fCtg = fCpt = fCphi = 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void AliBarrelTrack::SetLabel(Int_t label) {
  //
  // Sets the label
  //

  fLabel = label;  
}

////////////////////////////////////////////////////////////////////////////////

void AliBarrelTrack::SetX(Double_t x, Double_t alpha) {
  //
  
  fX = x;
  fAlpha = alpha;
}

////////////////////////////////////////////////////////////////////////////////

void AliBarrelTrack::SetRefPlane(Int_t nRefPlane, Int_t isIn) 
{
  //
  // Define the reference plane
  //
  fRefPlane = nRefPlane;
  fIsIn = isIn;
}

////////////////////////////////////////////////////////////////////////////////

void AliBarrelTrack::SetNClusters(Int_t nClusters, Double_t chi2) 
{
  //
  // Set number of cluster
  //
  fNClusters = nClusters;
  fChi2 = chi2;
}

////////////////////////////////////////////////////////////////////////////////

void AliBarrelTrack::SetTime(Double_t time[5], Double_t length) 
{
  //
  // Set time for a track
  //
  for(Int_t i=0; i<5; i++)
    fTimeHypothesis[i] = time[i];

  fLength = length;
}

////////////////////////////////////////////////////////////////////////////////

void AliBarrelTrack::SetStateVector(Double_t vec[5]) {
  //
  // Set State Vector from external representation
  //


  fY = vec[0];
  fZ = vec[1];
  fTgLambda = vec[3];
  fSnPhi = vec[2]; 
  f1Pt = vec[4];
}

////////////////////////////////////////////////////////////////////////////////

void AliBarrelTrack::SetCovarianceMatrix(Double_t vec[15]) {
  //
  // Set Covariance Matrix from external represenatation
  // 

  fCy   = vec[0];
  fCz   = vec[2];
  fCtg  = vec[9];
  fCphi = vec[5];
  fCpt  = vec[14];
}

////////////////////////////////////////////////////////////////////////////////
