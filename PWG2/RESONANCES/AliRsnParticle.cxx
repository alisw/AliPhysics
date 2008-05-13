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

// ***********************
// *** AliRsnParticle ****
// ***********************
// 
// This is a light-weight class which contains 
// informations from a MC particle, needed for resonance analysis.
//
// author: A. Pulvirenti --- email: alberto.pulvirenti@ct.infn.it


#include <TParticle.h>

#include "AliLog.h"

#include "AliRsnParticle.h"

ClassImp(AliRsnParticle)

//_____________________________________________________________________________
AliRsnParticle::AliRsnParticle() :
  fPDG(0),
  fMother(-1),
  fMotherPDG(0)
{
//=========================================================
// Default constructor.
// Initializes all data-members with meaningless values.
//=========================================================
}

//_____________________________________________________________________________
void AliRsnParticle::Adopt(TParticle *particle)
{
//=========================================================
// Copies data from a TParticle into "this":
//  - PDG code
//  - GEANT label of mother (if any, otherwise -1)
// If the argument is NULL, nothing is done, and an alert
// is given by the method.
//=========================================================
	
	if (!particle) {
	   AliError("NULL argument passed. Nothing done.");
	   return;
    }
    
    fPDG    = particle->GetPdgCode();
    fMother = particle->GetFirstMother();
    fMotherPDG = (Short_t)0;
}
