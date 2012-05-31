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

//-----------------------------------------------------------------------
//    Base Hit class for all detectors
//    Contains the coordinates of the hit (single energy deposition)
//    and the number of correspondent track
//    Author:
//-----------------------------------------------------------------------

#include "TParticle.h"

#include "AliHit.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliStack.h"

ClassImp(AliHit)

//_______________________________________________________________________
AliHit::AliHit():
  fTrack(0),
  fX(0),
  fY(0),
  fZ(0)
{
  //
  // Default constructor
  //
}

//_______________________________________________________________________
AliHit::AliHit(Int_t shunt, Int_t track):
  fTrack(0),
  fX(0),
  fY(0),
  fZ(0)
{
  //
  // Standard constructor
  //
  if(shunt == 1) {
    int primary = gAlice->GetMCApp()->GetPrimary(track);
    gAlice->GetMCApp()->Particle(primary)->SetBit(kKeepBit);
    fTrack=primary;
  } 

  else if (shunt == 2) {
    // the "primary" particle associated to the hit is
    // the last track that has been flagged in the StepManager
    // used by PHOS to associate the hit with the decay gamma
    // rather than with the original pi0 
    TParticle *part;
    Int_t current;
    Int_t parent=track;
    while (1) {
      current=parent;
      part = gAlice->GetMCApp()->Particle(current);
      parent=part->GetFirstMother();    
      if(parent<0 || part->TestBit(kKeepBit))
	break;
    }
    fTrack=current;   
  } else {
    fTrack=track;
    gAlice->GetMCApp()->FlagTrack(fTrack);
  }
}
