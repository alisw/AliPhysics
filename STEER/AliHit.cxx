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

/*
$Log$
Revision 1.4  2000/07/11 18:24:59  fca
Coding convention corrections + few minor bug fixes

Revision 1.3  1999/09/29 09:24:29  fca
Introduction of the Copyright and cvs Log

*/

#include "AliHit.h"
#include "TParticle.h"
#include "AliRun.h"
 
ClassImp(AliHit)

AliHit::AliHit()
{
  //
  // Default constructor
  //
  fTrack=0;	
}

AliHit::AliHit(Int_t shunt, Int_t track)
{
  //
  // Standard constructor
  //
  TClonesArray &particles = *(gAlice->Particles());
  if(shunt) {
    int primary = gAlice->GetPrimary(track);
    ((TParticle *)particles[primary])->SetBit(kKeepBit);
    fTrack=primary;
  } else {
    fTrack=track;
    gAlice->FlagTrack(fTrack);
  }
}

	 
