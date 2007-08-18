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

// volume: 
//  [0] = module number 1-60 (1==>(0-0), 60 (5-9)
//  [1] = Plastic number: 0 (down) to 1 (up)
//
// hit
// [0] = PID
// [1-3] = x, y, z 
// [4] = time 
// [5-7] = px, py, pz
// [8] = energy 
// [9] = energy loss
// [10] = trak length in plastic

#include "AliACORDEhit.h"

#include <TMath.h>

#include "AliConst.h"

ClassImp(AliACORDEhit)

//____________________________________________________________________________
AliACORDEhit::AliACORDEhit()
  : AliHit(),
    fModule(0),
    fPlastic(0),
    fTrackId(0),
    fTime(0),
    fPx(0),
    fPy(0),
    fPz(0),
    fEloss(0),
    fEnergy(0),
    fTrkLength(0)
{
  //
  // default ctor for AliACORDEhit object
  //
}

//_____________________________________________________________________________
AliACORDEhit::AliACORDEhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits)
  : AliHit(shunt, track),
    fModule(vol[0]),
    fPlastic(vol[1]),
    fTrackId((Int_t) hits[0]),
    fTime(hits[4]),
    fPx(hits[5]),
    fPy(hits[6]),
    fPz(hits[7]),
    fEloss(hits[9]),
    fEnergy(hits[8]),
    fTrkLength(hits[10])
{
  //
  // Constructor of hit object
  //
  fX = hits[1];
  fY = hits[2];
  fZ = hits[3];
}


//_____________________________________________________________________________
AliACORDEhit::~AliACORDEhit()
{
  //
  // Default destructor.
  //
}


//_____________________________________________________________________________
Float_t AliACORDEhit::PolarAngle() const
{
  //
  //
  //
  //  return kRaddeg*TMath::ACos(-fPy/this->Energy());
  return kRaddeg*TMath::ACos(fPz/this->Energy());
}

//_____________________________________________________________________________
Float_t AliACORDEhit::AzimuthAngle() const
{
  //
  //
  //
  //  return kRaddeg*TMath::ATan2(-fPx, -fPz);
  return kRaddeg*TMath::ATan2(fPx, fPz);
}
