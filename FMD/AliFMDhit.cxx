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
#include "AliFMDhit.h"

ClassImp(AliFMDhit)

AliFMDhit::AliFMDhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track)
{
//Normal FMD  hit ctor
  fVolume = vol[0];
  fX=hits[0];
  fY=hits[1];
  fZ=hits[2];
  fPx=hits[3];
  fPy=hits[4];
  fPz=hits[5];
  fEdep=hits[6];
  fParticle=(Int_t)hits[7];
  fTime=hits[8];
}
 













