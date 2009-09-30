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
/************************************************************************
 *                                                                      *
 * Abstract PID base class for Detector PID classes                     *
 * Supplies detector PID classes with basic informations (i.e. Debug    *
 * Level)                                                               *
 *                                                                      *
 * Authors:                                                             *
 *   Markus Fasel <M.Fasel@gsi.de>                                      *
 *                                                                      *
 ************************************************************************/
#include <TParticle.h>

#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"

#include "AliHFEpidBase.h"

ClassImp(AliHFEpidBase)

//___________________________________________________________________
AliHFEpidBase::AliHFEpidBase(const Char_t *name):
  TNamed(name, ""),
  fMCEvent(0x0),
  fDebugLevel(0)
{
  //
  // Default constructor
  //
}

//___________________________________________________________________
AliHFEpidBase::AliHFEpidBase(const AliHFEpidBase &c):
  TNamed(),
  fMCEvent(0x0),
  fDebugLevel(0)
{
  //
  //Copy constructor
  //
  c.Copy(*this);
}

//___________________________________________________________________
AliHFEpidBase &AliHFEpidBase::operator=(const AliHFEpidBase &ref){
  //
  // Assignment operator
  //
  if(this != &ref){
    ref.Copy(*this);
  }

  return *this;
}

//___________________________________________________________________
void AliHFEpidBase::Copy(TObject &ref) const {
  AliHFEpidBase &target = dynamic_cast<AliHFEpidBase &>(ref);

  target.fMCEvent = fMCEvent;
  target.fDebugLevel = fDebugLevel;

  TNamed::Copy(ref);
}

//___________________________________________________________________
Int_t AliHFEpidBase::GetPdgCode(AliVParticle *track){
  //
  // returns the MC PDG code of the particle species
  //
  if(!fMCEvent) return 0;
  AliMCParticle *mctrack = 0x0;
  if(TString(track->IsA()->GetName()).CompareTo("AliESDtrack") == 0)
    mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs((dynamic_cast<AliESDtrack *>(track))->GetLabel())));
  else if(TString(track->IsA()->GetName()).CompareTo("AliMCParticle") == 0)
    mctrack = dynamic_cast<AliMCParticle *>(track);
  if(!mctrack) return 0;
  return mctrack->Particle()->GetPdgCode();
}
