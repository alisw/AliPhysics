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

#include "AliPMDhit.h"
#include <TClonesArray.h>
#include "Riostream.h"
#include "Rtypes.h"

ClassImp(AliPMDhit)
  
//_____________________________________________________________________________
AliPMDhit::AliPMDhit():
  fEnergy(0.),
  fTime(0.)
{
  for (Int_t i=0; i<6; i++)
    {
      fVolume[i] = 0;
    }
}
//_____________________________________________________________________________
AliPMDhit::AliPMDhit(Int_t shunt,Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track),
  fEnergy(hits[3]),
  fTime(hits[4])
{
  //
  // Add a PMD hit
  //
  Int_t i;
  for (i=0; i<6; i++) fVolume[i] = vol[i];
  fX=hits[0];
  fY=hits[1];
  fZ=hits[2];
}
//_____________________________________________________________________________
AliPMDhit::AliPMDhit(AliPMDhit* oldhit):
  fEnergy(0.),
  fTime(0.)
{
  *this=*oldhit;
}

//_____________________________________________________________________________
int AliPMDhit::operator == (AliPMDhit &cell) const
{
  Int_t i;
  if(fTrack!=cell.GetTrack()) return 0;
  for (i=0; i<6; i++) if(fVolume[i]!=cell.GetVolume(i)) return 0;
  return 1;
}

void
AliPMDhit::Print(Option_t *) const {
   printf("PMD Cell %d %d %d %d %d %d\n   Primary %d -   Energy %f\n",
          fVolume[0],fVolume[1],fVolume[2],fVolume[3],
          fVolume[4],fVolume[5],fTrack,fEnergy);
}
