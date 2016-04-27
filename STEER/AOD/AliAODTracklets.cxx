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

#include <TString.h>
#include "AliAODTracklets.h"

ClassImp(AliAODTracklets)

AliAODTracklets::AliAODTracklets() 
: AliVMultiplicity(), fNTracks(0), fTheta(0), fPhi(0), fDeltaPhi(0), fLabels(0), fLabelsL2(0)
  ,fFastOrFiredChips(),fClusterFiredChips()
{
  fFiredChips[0] = fFiredChips[1] = 0;
  for (int i=6;i--;) fITSClusters[i] = 0;
  // default constructor
}

AliAODTracklets::AliAODTracklets(const char* name, const char* title)
: AliVMultiplicity(name, title), fNTracks(0), fTheta(0), fPhi(0), fDeltaPhi(0), fLabels(0), fLabelsL2(0)
, fFastOrFiredChips(),fClusterFiredChips()
{
  // Named constructor
  fFiredChips[0] = fFiredChips[1] = 0;
  for (int i=6;i--;) fITSClusters[i] = 0;
}

AliAODTracklets::AliAODTracklets(const AliAODTracklets& tracklet) :
    AliVMultiplicity(tracklet),
    fNTracks(tracklet.fNTracks),
    fTheta(0),
    fPhi(0),
    fDeltaPhi(0),
    fLabels(0), 
    fLabelsL2(0),
    fFastOrFiredChips(tracklet.fFastOrFiredChips),fClusterFiredChips(tracklet.fClusterFiredChips)
{
// Copy constructor
    fTheta = new Double32_t[fNTracks];
    fPhi = new Double32_t[fNTracks];
    fDeltaPhi = new Double32_t[fNTracks];
    fLabels = new Int_t[fNTracks];
    fLabelsL2 = new Int_t[fNTracks];
    for (Int_t i = 0; i < fNTracks; i++) {
	fTheta[i]    = tracklet.fTheta[i];
	fPhi[i]      = tracklet.fPhi[i];
	fDeltaPhi[i] = tracklet.fDeltaPhi[i];
	fLabels[i]   = tracklet.fLabels[i];
	fLabelsL2[i]   = tracklet.fLabelsL2[i];
    }
    fFiredChips[0] = tracklet.fFiredChips[0];
    fFiredChips[1] = tracklet.fFiredChips[1];
    for (int i=6;i--;) fITSClusters[i] = tracklet.fITSClusters[i];
}

AliAODTracklets& AliAODTracklets::operator=(const AliAODTracklets& tracklet)
{
// Assignment operator
    if(&tracklet == this) return *this;
    AliVMultiplicity::operator=(tracklet);
    if(fNTracks!=tracklet.fNTracks){
      fNTracks = tracklet.fNTracks;
      CreateContainer(fNTracks);
    }
    for (Int_t i = 0; i < fNTracks; i++) {
	fTheta[i]    = tracklet.fTheta[i];
	fPhi[i]      = tracklet.fPhi[i];
	fDeltaPhi[i] = tracklet.fDeltaPhi[i];
	fLabels[i]   = tracklet.fLabels[i];
	fLabelsL2[i]   = tracklet.fLabelsL2[i];
    }
    fFiredChips[0] = tracklet.fFiredChips[0];
    fFiredChips[1] = tracklet.fFiredChips[1];
    fFastOrFiredChips = tracklet.fFastOrFiredChips;
    fClusterFiredChips = tracklet.fClusterFiredChips;
    for (int i=6;i--;) fITSClusters[i] = tracklet.fITSClusters[i];
    return *this;
}

void AliAODTracklets::CreateContainer(Int_t nTracks)
{
  // function that creates container to store tracklets

  DeleteContainer();
  
  fNTracks = nTracks;

  if (fNTracks <= 0) {
    fNTracks = 0;
    return;
  }

  fTheta = new Double32_t[fNTracks];
  fPhi = new Double32_t[fNTracks];
  fDeltaPhi = new Double32_t[fNTracks];
  fLabels = new Int_t[fNTracks];
  fLabelsL2 = new Int_t[fNTracks];
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

  if (fLabelsL2)
  {
    delete[] fLabelsL2;
    fLabelsL2 = 0;
  }

  fNTracks = 0;
}

Bool_t AliAODTracklets::SetTracklet(Int_t pos, Double32_t theta, Double32_t phi, Double32_t deltaPhi, Int_t labelL1, Int_t labelL2)
{
  // Sets a tracklet at the given position

  if (pos < 0 || pos >= fNTracks)
    return kFALSE;

  fTheta[pos] = theta;
  fPhi[pos] = phi;
  fDeltaPhi[pos] = deltaPhi;
  fLabels[pos] = labelL1;
  fLabelsL2[pos] = labelL2;

  return kTRUE;
}

//______________________________________________________________________
void AliAODTracklets::Print(Option_t *opt) const
{
  // print
  printf("N.tracklets: %4d | ScaleDThtSin2T:%s\n",fNTracks,GetScaleDThetaBySin2T() ? "ON":"OFF");
  TString opts = opt; opts.ToLower();
  //
  if (opts.Contains("t")) {
    for (int i=0;i<fNTracks;i++) {
      printf("T#%3d| Eta:%+5.2f Th:%+6.3f Phi:%+6.3f DPhi:%+6.3f L1:%5d L2:%5d\n",
	     i,GetEta(i),fTheta[i],fPhi[i],fDeltaPhi[i],fLabels[i],fLabelsL2[i]);
    }
  }
  //
}

//________________________________________________________________
void AliAODTracklets::SetLabel(Int_t i, Int_t layer,Int_t label)  
{
  // set labels
  if (i>=0 && i<fNTracks) 
  {
    if(layer == 0) fLabels[i] = label;
    else fLabelsL2[i] = label;
  }
}

//________________________________________________________________
Int_t AliAODTracklets::GetLabel(Int_t i, Int_t layer) const 
{
  // access labels
  if (i>=0 && i<fNTracks) 
  {
    return (layer == 0) ? fLabels[i] : fLabelsL2[i];
  }
  else 
    Error("GetLabel","Invalid track number %d",i); return -9999;
}

//________________________________________________________________
Double_t AliAODTracklets::GetTheta(Int_t i) const 
{ 
  // access theta's
  if (i>=0 && i<fNTracks) 
  {
    return fTheta[i];
  }
  else 
    Error("GetTheta","Invalid track number %d",i); return -9999.;
}

//________________________________________________________________
Double_t AliAODTracklets::GetPhi(Int_t i) const 
{ 
  // access phi's
  if (i>=0 && i<fNTracks) 
  {
    return fPhi[i];
  }
  else 
    Error("GetPhi","Invalid track number %d",i); return -9999.;
}

//________________________________________________________________
Double_t AliAODTracklets::GetDeltaPhi(Int_t i) const 
{
  // access delta phi's
  if (i>=0 && i<fNTracks) 
  {
    return fDeltaPhi[i];
  }
  else 
    Error("GetDeltaPhi","Invalid track number %d",i); return -9999.;
}
