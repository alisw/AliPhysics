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
//
// Debug event
// the tree is represented as reduced events
// 
// Authors:
//   M.Fasel <M.Fasel@gsi.de>
//

#include "TObjArray.h"
#include "AliHFEreducedTrack.h"
#include "AliHFEreducedMCParticle.h"
#include "AliHFEreducedEvent.h"

ClassImp(AliHFEreducedEvent)

//_______________________________________
AliHFEreducedEvent::AliHFEreducedEvent():
TObject(),
  fTracks(NULL),
  fMCparticles(NULL),
  fNtracks(0),
  fNmcparticles(0),
  fRunNumber(0),
  fTrigger(0),
  fVX(0.),
  fVY(0.),
  fVZ(0.),
  fNContrib(0),
  fSPDMultiplicity(0)
{
  //
  // Default constructor
  //
  fTracks = new TObjArray;
  fTracks->SetOwner();
  fMCparticles = new TObjArray;
  fMCparticles->SetOwner();
  memset(fCentrality, 0, sizeof(Float_t) * 6);
  memset(fV0Multiplicity, 0, sizeof(Float_t) * 2);
  memset(fZDCEnergy, 0, sizeof(Float_t) * 4);
}

//_______________________________________
AliHFEreducedEvent::AliHFEreducedEvent(const AliHFEreducedEvent &ref):
  TObject(ref),
  fTracks(NULL),
  fMCparticles(NULL),
  fNtracks(ref.fNtracks),
  fNmcparticles(ref.fNmcparticles),
  fRunNumber(ref.fRunNumber),
  fTrigger(ref.fTrigger),
  fVX(ref.fVX),
  fVY(ref.fVY),
  fVZ(ref.fVZ),
  fNContrib(ref.fNContrib),
  fSPDMultiplicity(ref.fSPDMultiplicity)
{
  //
  // Copy constructor
  //
  fTracks = new TObjArray;
  fTracks->SetOwner();
  for(int itrk = 0; itrk < ref.GetNumberOfTracks(); itrk++)
    fTracks->Add(new AliHFEreducedTrack(*(ref.GetTrack(itrk))));
  fMCparticles = new TObjArray;
  fMCparticles->SetOwner();
  for(int iprt = 0; iprt < ref.GetNumberOfMCParticles(); iprt++)
    fMCparticles->Add(new AliHFEreducedMCParticle(*(ref.GetMCParticle(iprt))));
  memcpy(fCentrality, ref.fCentrality, sizeof(Float_t) * 6);
  memcpy(fV0Multiplicity, ref.fV0Multiplicity, sizeof(Float_t) * 2);
  memcpy(fZDCEnergy, ref.fZDCEnergy, sizeof(Float_t) *4);
}

//_______________________________________
AliHFEreducedEvent &AliHFEreducedEvent::operator=(const AliHFEreducedEvent &ref){
  //
  // Assignment operator
  //
  if(&ref != this){
    TObject::operator=(ref);
    fTracks->SetOwner();
    fTracks->Clear();
    for(int itrk = 0; itrk < ref.GetNumberOfTracks(); itrk++)
      fTracks->Add(new AliHFEreducedTrack(*(ref.GetTrack(itrk))));
    fMCparticles->Clear();
    fMCparticles->SetOwner();
    for(int iprt = 0; iprt < ref.GetNumberOfMCParticles(); iprt++)
      fMCparticles->Add(new AliHFEreducedMCParticle(*(ref.GetMCParticle(iprt))));
    fNtracks = ref.fNtracks;
    fNmcparticles = ref.fNmcparticles;
    fRunNumber = ref.fRunNumber;
    fTrigger = ref.fTrigger;
    fVX = ref.fVX;
    fVY = ref.fVY;
    fVZ = ref.fVZ;
    fSPDMultiplicity = ref.fSPDMultiplicity;
    memcpy(fCentrality, ref.fCentrality, sizeof(Float_t) * 6);
    memcpy(fV0Multiplicity, ref.fV0Multiplicity, sizeof(Float_t) * 2);
    memcpy(fZDCEnergy, ref.fZDCEnergy, sizeof(Float_t) *4);
  }
  return *this;
}

//_______________________________________
AliHFEreducedEvent::~AliHFEreducedEvent(){
  //
  // Destructor: Clear tracks an MC particles
  //
  delete fTracks;
  delete fMCparticles;
}

//_______________________________________
void AliHFEreducedEvent::AddTrack(const AliHFEreducedTrack *track){
  //
  // Add track to the event
  //
  fTracks->Add(new AliHFEreducedTrack(*track));
  fNtracks++;
}

//_______________________________________
const AliHFEreducedTrack *AliHFEreducedEvent::GetTrack(Int_t itrk) const {
  //
  // Get Track
  //
  if(itrk < 0 || itrk >= fNtracks) return NULL;
  return dynamic_cast<const AliHFEreducedTrack *>(fTracks->At(itrk));
}

//_______________________________________
void AliHFEreducedEvent::AddMCParticle(const AliHFEreducedMCParticle *track){
  //
  // Add MC particle to the Event
  //
  fMCparticles->Add(new AliHFEreducedMCParticle(*track));
  fNmcparticles++;
}

//_______________________________________
const AliHFEreducedMCParticle *AliHFEreducedEvent::GetMCParticle(Int_t itrk) const {
  //
  // Get MC particle
  //
  if(itrk < 0 || itrk >= fNmcparticles) return NULL;
  return dynamic_cast<const AliHFEreducedMCParticle *>(fMCparticles->At(itrk));
}
