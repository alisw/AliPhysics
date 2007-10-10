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
//     Realisation of AliVParticle for MC Particles
//     Implementation wraps a TParticle and delegates the methods
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TRefArray.h>

#include "AliMCParticle.h"


ClassImp(AliMCParticle)

AliMCParticle::AliMCParticle():
    AliVParticle(),
    fParticle(0),
    fTrackReferences(0),
    fNTrackRef(0)
{
    // Constructor
}

AliMCParticle::AliMCParticle(TParticle* part, TRefArray* rarray):
    AliVParticle(),
    fParticle(part),
    fTrackReferences(rarray),
    fNTrackRef(0)
{
    // Constructor
    if (rarray != 0) {
	fNTrackRef = fTrackReferences->GetEntriesFast();
    }
}
    
    
AliMCParticle::AliMCParticle(const AliMCParticle& mcPart) :
    AliVParticle(mcPart),
    fParticle(0)
{
// Copy constructor
}

AliMCParticle& AliMCParticle::operator=(const AliMCParticle& mcPart)
{ if (this!=&mcPart) { 
    AliVParticle::operator=(mcPart); 
  }
  
  return *this; 
}
