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
Revision 1.6  2001/01/26 19:58:48  hristov
Major upgrade of AliRoot code

Revision 1.5  2000/07/12 08:56:25  fca
Coding convention correction and warning removal

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
  if(shunt == 1) {
    int primary = gAlice->GetPrimary(track);
    gAlice->Particle(primary)->SetBit(kKeepBit);
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
      part = gAlice->Particle(current);
      parent=part->GetFirstMother();    
      if(parent<0 || part->TestBit(kKeepBit))
	break;
    }
    fTrack=current;   
  }

  else {
    fTrack=track;
    gAlice->FlagTrack(fTrack);
  }
}

	 
