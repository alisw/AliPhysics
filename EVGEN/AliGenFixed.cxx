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

// Simple particle gun. 
// Momentum, phi and theta of the partice as well as the particle type can be set.
// If fExplicit is true the user set momentum vector is used,
// otherwise it is calculated.
// andreas.morsch@cern.ch

#include "TPDGCode.h"

#include "AliGenFixed.h"
#include "AliRun.h"
  
ClassImp(AliGenFixed)

//_____________________________________________________________________________
AliGenFixed::AliGenFixed()
    :AliGenerator(), 
     fIpart(0),
     fExplicit(kFALSE)
{
  //
  // Default constructor
  //
}

//_____________________________________________________________________________
AliGenFixed::AliGenFixed(Int_t npart)
    :AliGenerator(npart),
     fIpart(kProton),
     fExplicit(kFALSE)
{
  //
  // Standard constructor
  //
  fName="Fixed";
  fTitle="Fixed Particle Generator";
}

//_____________________________________________________________________________
void AliGenFixed::Generate()
{
  //
  // Generate one trigger
  //
  Float_t polar[3]= {0,0,0};
  if(!fExplicit) {
    fP[0] = fPMin*TMath::Cos(fPhiMin)*TMath::Sin(fThetaMin);
    fP[1] = fPMin*TMath::Sin(fPhiMin)*TMath::Sin(fThetaMin);
    fP[2] = fPMin*TMath::Cos(fThetaMin);
  }
  Int_t i, nt;
  //
  Float_t o[3];
  o[0] = fOrigin[0];
  o[1] = fOrigin[1];
  o[2] = fOrigin[2];

  printf("Origin %f %f %f \n", o[0], o[1], o[2]);
  
  
  for(i=0;i<fNpart;i++) 
    PushTrack(fTrackIt,-1,fIpart,fP,o ,polar,0,kPPrimary,nt);
}
  
//_____________________________________________________________________________
void AliGenFixed::SetSigma(Float_t /*sx*/, Float_t /*sy*/, Float_t /*sz*/)
{
  //
  // Set the interaction point sigma
  //
  printf("Vertex smearing not implemented for fixed generator\n");
}
