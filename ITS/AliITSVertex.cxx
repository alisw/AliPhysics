/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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
//           Implementation of the Primary Vertex class
//
// Origin: A.Dainese, Padova, andrea.dainese@pd.infn.it
//-----------------------------------------------------------------

//---- standard headers ----
#include "Riostream.h"
//---- Root headers --------
#include <TMath.h>
#include <TROOT.h>
#include <TNamed.h>
//---- AliRoot headers -----
#include "AliITSVertex.h"


ClassImp(AliITSVertex)

//--------------------------------------------------------------------------
AliITSVertex::AliITSVertex() {
//
// Default Constructor, set everything to 0
//
  SetToZero();
}
//--------------------------------------------------------------------------
AliITSVertex::AliITSVertex(Double_t positionZ,Double_t sigmaZ,
			   Int_t nContributors, const char *vtxName) {
  //
  // Constructor for vertex Z from pixels
  //

  SetToZero();

  fPosition[2]   = positionZ;
  fCovZZ         = sigmaZ*sigmaZ;
  fNContributors = nContributors;
  SetName(vtxName);

}
//------------------------------------------------------------------------- 
AliITSVertex::AliITSVertex(Double_t phi,
			   Double_t position[3],Double_t covmatrix[6],
			   Double_t chi2,Int_t nContributors, const char *vtxName) {
//
// Constructor for vertex in 3D from tracks
//

  SetToZero();
    fPosition[0]   = position[0];
    fPosition[1]   = position[1];
    fPosition[2]   = position[2];
    fCovXX         = covmatrix[0];
    fCovXY         = covmatrix[1];
    fCovYY         = covmatrix[2];
    fCovXZ         = covmatrix[3];
    fCovYZ         = covmatrix[4];
    fCovZZ         = covmatrix[5];


    fPhi           = phi;
    fChi2          = chi2;
    fNContributors = nContributors;

    SetName(vtxName);

}
//--------------------------------------------------------------------------
AliITSVertex::AliITSVertex(Double_t position[3],Double_t sigma[3],
			   const char *vtxName) {
//
// Constructor for smearing of true position
//

  SetToZero();
    fPosition[0]   = position[0];
    fPosition[1]   = position[1];
    fPosition[2]   = position[2];
    fCovXX         = sigma[0]*sigma[0];
    fCovXY         = 0;
    fCovYY         = sigma[1]*sigma[1];
    fCovXZ         = 0;
    fCovYZ         = 0;
    fCovZZ         = sigma[2]*sigma[2];


    SetName(vtxName);

}
//--------------------------------------------------------------------------
AliITSVertex::AliITSVertex(Double_t position[3],Double_t sigma[3],
			   Double_t snr[3], const char *vtxName) {
  //
  // Constructor for Pb-Pb
  //

  SetToZero();
  fPosition[0]   = position[0];
  fPosition[1]   = position[1];
  fPosition[2]   = position[2];
  fCovXX         = sigma[0]*sigma[0];
  fCovXY         = 0;
  fCovYY         = sigma[1]*sigma[1];
  fCovXZ         = 0;
  fCovYZ         = 0;
  fCovZZ         = sigma[2]*sigma[2];

  fSNR[0]        = snr[0];
  fSNR[1]        = snr[1];
  fSNR[2]        = snr[2];

  SetName(vtxName);

}
//--------------------------------------------------------------------------
void AliITSVertex::SetToZero() {
  //
  // Set some data members to 0. Used by constructors
  //
  for(Int_t i=0; i<3; i++){
    fPosition[i] = 0.;
    fTruePos[i] = 0;
    fSNR[i] = 0.;
  }
  fCovXX         = 0;
  fCovXY         = 0;
  fCovYY         = 0;
  fCovXZ         = 0;
  fCovYZ         = 0;
  fCovZZ         = 0;

  fPhi           = 0;
  fChi2          = 0;
  fNContributors = 0;

  SetDebug();
}
//--------------------------------------------------------------------------
AliITSVertex::~AliITSVertex() {
//  
// Default Destructor
//

}
//--------------------------------------------------------------------------
void AliITSVertex::GetXYZ(Double_t position[3]) const {
//
// Return position of the vertex in global frame
//
  position[0] = fPosition[0]*TMath::Cos(fPhi)-fPosition[1]*TMath::Sin(fPhi);
  position[1] = fPosition[0]*TMath::Sin(fPhi)+fPosition[1]*TMath::Cos(fPhi);
  position[2] = fPosition[2];

  return;
}
//--------------------------------------------------------------------------
void AliITSVertex::GetXYZ_ThrustFrame(Double_t position[3]) const {
//
// Return position of the vertex in thrust frame
//
  position[0] = fPosition[0];
  position[1] = fPosition[1];
  position[2] = fPosition[2];

  return;
}
//--------------------------------------------------------------------------
void AliITSVertex::GetSigmaXYZ(Double_t sigma[3]) const {
//
// Return errors on vertex position in global frame
//
  Double_t cm[6];
  GetCovMatrix(cm);
  sigma[0] = TMath::Sqrt(cm[0]);
  sigma[1] = TMath::Sqrt(cm[2]);
  sigma[2] = TMath::Sqrt(cm[5]);

  return;
}
//--------------------------------------------------------------------------
void AliITSVertex::GetSigmaXYZ_ThrustFrame(Double_t sigma[3]) const {
//
// Return errors on vertex position in thrust frame
//
  sigma[0] = TMath::Sqrt(fCovXX);
  sigma[1] = TMath::Sqrt(fCovYY);
  sigma[2] = TMath::Sqrt(fCovZZ);

  return;
}
//--------------------------------------------------------------------------
void AliITSVertex::GetCovMatrix(Double_t covmatrix[6]) const {
//
// Return covariance matrix of the vertex
//
  Double_t cPhi = TMath::Cos(fPhi);
  Double_t sPhi = TMath::Sin(fPhi);

  covmatrix[0] = fCovXX*cPhi*cPhi-2.*fCovXY*cPhi*sPhi+fCovYY*sPhi*sPhi;
  covmatrix[1] = (fCovXX-fCovYY)*cPhi*sPhi+fCovXY*(cPhi*cPhi-sPhi*sPhi);
  covmatrix[2] = fCovXX*sPhi*sPhi+2.*fCovXY*cPhi*sPhi+fCovYY*cPhi*cPhi;
  covmatrix[3] = fCovXZ*cPhi-fCovYZ*sPhi;
  covmatrix[4] = fCovXZ*sPhi+fCovYZ*cPhi;
  covmatrix[5] = fCovZZ;

  return;
}
//--------------------------------------------------------------------------
void AliITSVertex::GetCovMatrix_ThrustFrame(Double_t covmatrix[6]) const {
//
// Return covariance matrix of the vertex
//
  covmatrix[0] = fCovXX;
  covmatrix[1] = fCovXY;
  covmatrix[2] = fCovYY;
  covmatrix[3] = fCovXZ;
  covmatrix[4] = fCovYZ;
  covmatrix[5] = fCovZZ;

  return;
}
//--------------------------------------------------------------------------
Double_t AliITSVertex::GetXv() const {
//
// Return global x
//
  return fPosition[0]*TMath::Cos(fPhi)-fPosition[1]*TMath::Sin(fPhi);
}
//--------------------------------------------------------------------------
Double_t AliITSVertex::GetYv() const {
//
// Return global y
//
  return fPosition[0]*TMath::Sin(fPhi)+fPosition[1]*TMath::Cos(fPhi);
}
//--------------------------------------------------------------------------
Double_t AliITSVertex::GetZv() const {
//
// Return global z
//
  return fPosition[2];
}
//--------------------------------------------------------------------------
Double_t AliITSVertex::GetXRes() const {
//
// Return error on global x
//
  Double_t cm[6];
  GetCovMatrix(cm);
  return TMath::Sqrt(cm[0]);
}
//--------------------------------------------------------------------------
Double_t AliITSVertex::GetYRes() const {
//
// Return error on global y
//
  Double_t cm[6];
  GetCovMatrix(cm);
  return TMath::Sqrt(cm[2]);
}
//--------------------------------------------------------------------------
Double_t AliITSVertex::GetZRes() const {
//
// Return error on global z
//
  Double_t cm[6];
  GetCovMatrix(cm);
  return TMath::Sqrt(cm[5]);
}
//--------------------------------------------------------------------------
void AliITSVertex::GetSNR(Double_t snr[3]) const {
//
// Return S/N ratios
//
  for(Int_t i=0;i<3;i++) snr[i] = fSNR[i];

  return;
}
//--------------------------------------------------------------------------
void AliITSVertex::PrintStatus() const {
//
// Print out information on all data members
//
  if(fPhi) {
    printf(" ! The vertex fitting has been done in the thrust frame !\n");
    printf("   The rotation angle is %f. Numbers are given in the rotated frame.\n",fPhi);
  }
  printf(" Vertex position:\n");
  printf("   x = %f +- %f\n",fPosition[0],TMath::Sqrt(fCovXX));
  printf("   y = %f +- %f\n",fPosition[1],TMath::Sqrt(fCovYY));
  printf("   z = %f +- %f\n",fPosition[2],TMath::Sqrt(fCovZZ));
  printf(" Covariance matrix:\n");
  printf(" %12.10f  %12.10f  %12.10f\n %12.10f  %12.10f  %12.10f\n %12.10f  %12.10f  %12.10f\n",fCovXX,fCovXY,fCovXZ,fCovXY,fCovYY,fCovYZ,fCovXZ,fCovYZ,fCovZZ);
  printf(" S/N = (%f, %f, %f)\n",fSNR[0],fSNR[1],fSNR[2]);
  printf(" chi2 = %f\n",fChi2);
  printf(" # tracks (or tracklets) = %d\n",fNContributors);

  printf(" True vertex position - for comparison: %12.10f  %12.10f  %12.10f\n ",fTruePos[0],fTruePos[1],fTruePos[2]);

  return;
}




