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

#include "AliMCParticle.h"


ClassImp(AliMCParticle)

AliMCParticle::AliMCParticle():
    AliVParticle(),
    fParticle(0)
{
    // Constructor
}

AliMCParticle::AliMCParticle(TParticle* part):
    AliVParticle(),
    fParticle(part)
{
    // Constructor
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
