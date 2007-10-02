/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//     AOD class to store tracklets
//     Author: Jan Fiete Grosse-Oetringhaus, CERN
//     Class created from AliMultiplicity
//-------------------------------------------------------------------------

#include "AliAODTracklets.h"

ClassImp(AliAODTracklets)

AliAODTracklets::AliAODTracklets() : TNamed(), fNTracks(0), fTheta(0), fPhi(0), fDeltaPhi(0), fLabels(0)
{
  // default constructor
}

AliAODTracklets::AliAODTracklets(const char* name, const char* title) : TNamed(name, title), fNTracks(0), fTheta(0), fPhi(0), fDeltaPhi(0), fLabels(0)
{
  // TNamed constructor
}

void AliAODTracklets::CreateContainer(Int_t nTracks)
{
  // function that creates container to store tracklets

  DeleteContainer();
  
  fNTracks = nTracks;

  if (fNTracks <= 0)
    return;

  fTheta = new Double32_t[fNTracks];
  fPhi = new Double32_t[fNTracks];
  fDeltaPhi = new Double32_t[fNTracks];
  fLabels = new Int_t[fNTracks];
}

AliAODTracklets::~AliAODTracklets()
{
  // destructor

  DeleteContainer();
}

void AliAODTracklets::DeleteContainer()
{
  // deletes allocated memory

  if (fTheta)
  {
    delete[] fTheta;
    fTheta = 0;
  }

  if (fPhi)
  {
    delete[] fPhi;
    fPhi = 0;
  }

  if (fDeltaPhi)
  {
    delete[] fDeltaPhi;
    fDeltaPhi = 0;
  }

  if (fLabels)
  {
    delete[] fLabels;
    fLabels = 0;
  }

  fNTracks = 0;
}

Bool_t AliAODTracklets::SetTracklet(Int_t pos, Double32_t theta, Double32_t phi, Double32_t deltaPhi, Int_t label)
{
  // Sets a tracklet at the given position

  if (pos < 0 || pos >= fNTracks)
    return kFALSE;

  fTheta[pos] = theta;
  fPhi[pos] = phi;
  fDeltaPhi[pos] = deltaPhi;
  fLabels[pos] = label;

  return kTRUE;
}
