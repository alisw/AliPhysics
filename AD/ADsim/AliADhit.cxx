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

/* $Id: AliADhit.cxx  $ */

//_________________________________________________________________________
//
//      Hit class for AD detector
//
//_________________________________________________________________________

#include "AliADhit.h"

ClassImp(AliADhit);

//_____________________________________________________________________________
AliADhit::AliADhit()
  : AliHit()
  , fModule(0)
  , fEk(0.)
  , fPt(0.)
  , fPx(0.)
  , fPy(0.)
  , fPz(0.)
  , fTof(0.)
  , fTleng(0.)
  , fEloss(0.)
  , fNphot(0)
  , fCell(0)
  , fPrimary(0)
  , fPDG(0)
  , fPDGMother(0)
{
  // Default constructor
}

//_____________________________________________________________________________
AliADhit::AliADhit(Int_t shunt, Int_t track, Int_t* vol, Float_t* hits)
  : AliHit(shunt, track)
  , fModule(vol[4])
  , fEk(hits[3])
  , fPt(hits[4])
  , fPx(hits[5])
  , fPy(hits[6])
  , fPz(hits[7])
  , fTof(hits[8])
  , fTleng(hits[9])
  , fEloss(hits[10])
  , fNphot(vol[3])
  , fCell(vol[4])
  , fPrimary(vol[0])
  , fPDG(vol[1])
  , fPDGMother(vol[2])
{
  fX = hits[0];      // X position of hit
  fY = hits[1];      // Y position of hit
  fZ = hits[2];      // Z position of hit
}
//_____________________________________________________________________________
AliADhit::~AliADhit()
{
  //
  // Default destructor.
  //
}
