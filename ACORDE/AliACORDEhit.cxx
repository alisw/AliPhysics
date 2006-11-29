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

#include "AliACORDEhit.h"

#include <TMath.h>

#include "AliConst.h"

ClassImp(AliACORDEhit)

//____________________________________________________________________________
AliACORDEhit::AliACORDEhit()
  : AliHit(),
    fId(0),
    fPx(0),
    fPy(0),
    fPz(0),
    fEloss(0),
    fMedium(0)
{
  //
  // default ctor for AliACORDEhit object
  //
}

//_____________________________________________________________________________
AliACORDEhit::AliACORDEhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits)
  : AliHit(shunt, track),
    fId(hits[0]),
    fPx(hits[4]),
    fPy(hits[5]),
    fPz(hits[6]),
    fEloss(hits[7]),
    fMedium(vol[0])
{
  //
  // Constructor of hit object
  //
  fX = hits[1];
  fY = hits[2];
  fZ = hits[3];
}

//____________________________________________________________________________
AliACORDEhit::AliACORDEhit(const AliACORDEhit & hit)
  : AliHit(hit),
    fId(hit.fId),
    fPx(hit.fPx),
    fPy(hit.fPy),
    fPz(hit.fPz),
    fEloss(hit.fEloss),
    fMedium(hit.fMedium)
{
  //
  // copy ctor
  //
  fX      = hit.fX;
  fY      = hit.fY;
  fZ      = hit.fZ;
}

//_____________________________________________________________________________
AliACORDEhit::~AliACORDEhit()
{
  //
  // Default destructor.
  //
}

//_____________________________________________________________________________
AliACORDEhit& AliACORDEhit::operator=(const AliACORDEhit & hit)
{
  //
  // aisngment operator.
  //
  fId     = hit.fId;
  fX      = hit.fX;
  fY      = hit.fY;
  fZ      = hit.fZ;
  fPx     = hit.fPx;
  fPy     = hit.fPy;
  fPz     = hit.fPz;
  fEloss  = hit.fEloss;
  fMedium = hit.fMedium;
  return *this;
}

//_____________________________________________________________________________
Float_t AliACORDEhit::Energy() const
{
  //
  //
  //
  return TMath::Sqrt(fPx*fPx + fPy*fPy + fPz*fPz);
}

//_____________________________________________________________________________
Float_t AliACORDEhit::PolarAngle() const
{
  //
  //
  //
  return kRaddeg*TMath::ACos(-fPy/this->Energy());
}

//_____________________________________________________________________________
Float_t AliACORDEhit::AzimuthAngle() const
{
  //
  //
  //
  return kRaddeg*TMath::ATan2(-fPx, -fPz);
}

//_____________________________________________________________________________
Bool_t AliACORDEhit::operator==(const AliACORDEhit& hit)
{
  //
  //
  //
  Float_t energy = TMath::Sqrt(fPx*fPx + fPy*fPy + fPz*fPz);
  Float_t energy2=TMath::Sqrt(hit.fPx*hit.fPx+hit.fPy*hit.fPy+hit.fPz*hit.fPz);
  return (energy == energy2);
  //return (fTrack == hit.fTrack);
}

//_____________________________________________________________________________
Bool_t AliACORDEhit::operator<(const AliACORDEhit& hit)
{
  //
  //
  //
  Float_t energy = TMath::Sqrt(fPx*fPx + fPy*fPy + fPz*fPz);
  Float_t energy2=TMath::Sqrt(hit.fPx*hit.fPx+hit.fPy*hit.fPy+hit.fPz*hit.fPz);
  return (energy < energy2);
}
