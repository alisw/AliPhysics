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
    for (Int_t i = 0; i < 3; i++) fP[i] = 0.;
    
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
  for (Int_t i = 0; i < 3; i++) fP[i] = 0.;
}

//_____________________________________________________________________________
void AliGenFixed::Generate()
{
  //
  // Generate one trigger
  //
  Float_t polar[3]= {0,0,0};
  if(!fExplicit) {
    fP[0] = fPMin * TMath::Cos(fPhiMin) * TMath::Sin(fThetaMin);
    fP[1] = fPMin * TMath::Sin(fPhiMin) * TMath::Sin(fThetaMin);
    fP[2] = fPMin * TMath::Cos(fThetaMin);
  }
  Int_t i, j, nt;
  //
  Float_t o[3] = {0., 0., 0.}; 
  Float_t time = 0.;
  if(fVertexSmear == kPerEvent) {
      Vertex();
      for (j = 0;j < 3; j++) o[j] = fVertex[j];
      time = fTime;
  }
  
  for(i = 0; i < fNpart; i++) 
    PushTrack(fTrackIt, -1, fIpart, fP, o , polar, time, kPPrimary, nt);
}
  
