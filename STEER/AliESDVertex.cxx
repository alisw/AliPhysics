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
//           for the Event Data Summary class
//           This class contains the Primary Vertex
//           of the event coming from reconstruction
// Origin: A.Dainese, Padova, andrea.dainese@pd.infn.it
//-----------------------------------------------------------------

//---- standard headers ----
#include "Riostream.h"
//---- Root headers --------
#include <TMath.h>
#include <TROOT.h>
//---- AliRoot headers -----
#include "AliESDVertex.h"


ClassImp(AliESDVertex)

//--------------------------------------------------------------------------
AliESDVertex::AliESDVertex() {
//
// Default Constructor, set everything to 0
//
  SetToZero();
}
//--------------------------------------------------------------------------
AliESDVertex::AliESDVertex(Double_t positionZ,Double_t sigmaZ,
			   Int_t nContributors,const Char_t *vtxName) {
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
AliESDVertex::AliESDVertex(Double_t position[3],Double_t covmatrix[6],
			   Double_t chi2,Int_t nContributors, const Char_t *vtxName) {
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


    fChi2          = chi2;
    fNContributors = nContributors;

    SetName(vtxName);

}
//--------------------------------------------------------------------------
AliESDVertex::AliESDVertex(Double_t position[3],Double_t sigma[3],
			   const Char_t *vtxName) {
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
AliESDVertex::AliESDVertex(Double_t position[3],Double_t sigma[3],
			   Double_t snr[3], const Char_t *vtxName) {
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
void AliESDVertex::SetToZero() {
  //
  // Set some data members to 0. Used by constructors
  //
  for(Int_t i=0; i<3; i++){
    fPosition[i] = 0.;
    fTruePos[i] = 0;
    fSNR[i] = 0.;
  }
  fCovXX         = 0.005*0.005;
  fCovXY         = 0;
  fCovYY         = 0.005*0.005;
  fCovXZ         = 0;
  fCovYZ         = 0;
  fCovZZ         = 5.3*5.3;

  fChi2          = 0;
  fNContributors = 0;
}
//--------------------------------------------------------------------------
AliESDVertex::~AliESDVertex() {
//  
// Default Destructor
//

}
//--------------------------------------------------------------------------
void AliESDVertex::GetXYZ(Double_t position[3]) const {
//
// Return position of the vertex in global frame
//
  position[0] = fPosition[0];
  position[1] = fPosition[1];
  position[2] = fPosition[2];

  return;
}
//--------------------------------------------------------------------------
void AliESDVertex::GetSigmaXYZ(Double_t sigma[3]) const {
//
// Return errors on vertex position in thrust frame
//
  sigma[0] = TMath::Sqrt(fCovXX);
  sigma[1] = TMath::Sqrt(fCovYY);
  sigma[2] = TMath::Sqrt(fCovZZ);

  return;
}
//--------------------------------------------------------------------------
void AliESDVertex::GetCovMatrix(Double_t covmatrix[6]) const {
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
Double_t AliESDVertex::GetXv() const {
//
// Return global x
//
  return fPosition[0];
}
//--------------------------------------------------------------------------
Double_t AliESDVertex::GetYv() const {
//
// Return global y
//
  return fPosition[1];
}
//--------------------------------------------------------------------------
Double_t AliESDVertex::GetZv() const {
//
// Return global z
//
  return fPosition[2];
}
//--------------------------------------------------------------------------
Double_t AliESDVertex::GetXRes() const {
//
// Return error on global x
//
  return TMath::Sqrt(fCovXX);
}
//--------------------------------------------------------------------------
Double_t AliESDVertex::GetYRes() const {
//
// Return error on global y
//
  return TMath::Sqrt(fCovYY);
}
//--------------------------------------------------------------------------
Double_t AliESDVertex::GetZRes() const {
//
// Return error on global z
//
  return TMath::Sqrt(fCovZZ);
}
//--------------------------------------------------------------------------
void AliESDVertex::GetSNR(Double_t snr[3]) const {
//
// Return S/N ratios
//
  for(Int_t i=0;i<3;i++) snr[i] = fSNR[i];

  return;
}
//--------------------------------------------------------------------------
void AliESDVertex::Print(Option_t* /*option*/) const {
//
// Print out information on all data members
//
  printf("ESD vertex position:\n");
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




