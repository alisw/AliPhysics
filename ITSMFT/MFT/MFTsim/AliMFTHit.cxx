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

//====================================================================================================================================================
//
//      Hit description for the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TLorentzVector.h"
#include "TParticle.h"
#include "AliHit.h"
#include "AliRun.h"
#include "AliMC.h"

#include "AliMFTHit.h"

ClassImp(AliMFTHit)

//====================================================================================================================================================

AliMFTHit::AliMFTHit():
  AliHit(),
  fStatus(0),
  fPlane(-1),
  fDetElemID(-1),
  fPx(0),
  fPy(0),
  fPz(0),
  fEloss(0),
  fTOF(0) 
{

  // default constructor 

}

//====================================================================================================================================================

TParticle* AliMFTHit::GetParticle() const {

    //   The TParticle of the track that created this hit.

    return gAlice->GetMCApp()->Particle(GetTrack());

}

//====================================================================================================================================================
