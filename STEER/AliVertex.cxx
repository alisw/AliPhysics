/**************************************************************************
 * Copyright(c) 2006-2008, ALICE Experiment at CERN, All rights reserved. *
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
//           Implementation of the base Vertex class
//           This class contains the Secondary Vertex
//           of a set of tracks
//           And it is the base class for primary vertices
// Origin: F.Prino, Torino, prino@to.infn.it
//-----------------------------------------------------------------

#include "AliVertex.h"


ClassImp(AliVertex)

//--------------------------------------------------------------------------
AliVertex::AliVertex() : TNamed() {
//
// Default Constructor, set everything to 0
//
  for(Int_t k=0;k<3;k++) fPosition[k]   = 0;
  fSigma = 0;
  fNContributors=0;
}

//--------------------------------------------------------------------------
AliVertex::AliVertex(Double_t position[3],Double_t dispersion,
		Int_t nContributors): TNamed() {
  //
  // Standard Constructor
  //

  for(Int_t k=0;k<3;k++) fPosition[k]   = position[k];
  fSigma         = dispersion;
  fNContributors = nContributors;
  SetName("BaseVertex");

}

//--------------------------------------------------------------------------
AliVertex::AliVertex(const AliVertex &source): TNamed(source) {
  //
  // Copy constructor
  //
  for(Int_t i=0;i<3;i++)fPosition[i] = source.fPosition[i];
  fSigma = source.GetDispersion();
  fNContributors = source.GetNContributors();
}

//--------------------------------------------------------------------------
AliVertex &AliVertex::operator=(const AliVertex &source){
  //
  // assignment operator
  //
  if(&source == this) return *this;
  this->SetName(source.GetName());
  this->SetTitle(source.GetTitle());
  for(Int_t i=0;i<3;i++)fPosition[i] = source.fPosition[i];
  fSigma = source.GetDispersion();
  fNContributors = source.GetNContributors();
  return *this;
}


//--------------------------------------------------------------------------
AliVertex::~AliVertex() {
//  
// Default Destructor
//

}
//--------------------------------------------------------------------------
void AliVertex::GetXYZ(Double_t position[3]) const {
//
// Return position of the vertex in global frame
//
  position[0] = fPosition[0];
  position[1] = fPosition[1];
  position[2] = fPosition[2];

  return;
}
//--------------------------------------------------------------------------
void AliVertex::Print(Option_t* /*option*/) const {
//
// Print out information on all data members
//
  printf("Vertex position:\n");
  printf("   x = %f\n",fPosition[0]);
  printf("   y = %f\n",fPosition[1]);
  printf("   z = %f\n",fPosition[2]);
  printf(" Dispersion = %f\n",fSigma);
  printf(" # tracks = %d\n",fNContributors);

  return;
}




