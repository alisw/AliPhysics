
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

void AliBarrelTrack::SetRefPlane(Int_t nRefPlane, Int_t isIn) {
  
  fRefPlane = nRefPlane;
  fIsIn = isIn;
}

////////////////////////////////////////////////////////////////////////////////

void AliBarrelTrack::SetNClusters(Int_t nClusters, Double_t chi2) {

  fNClusters = nClusters;
  fChi2 = chi2;
}

////////////////////////////////////////////////////////////////////////////////

void AliBarrelTrack::SetTime(Double_t time[5], Double_t length) {

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
