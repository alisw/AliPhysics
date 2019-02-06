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
// Originator:  M.Fasel <M.Fasel@gsi.de>
// Base class for AliHFEminiEventCreator.cxx/h
// Contributors: Nirbhay K. Behera, Jiyeon Kwon, Jonghan Park

#include "TObjArray.h"
#include "AliHFEreducedMCParticle.h"
#include "AliHFEminiTrack.h"
#include "AliHFEminiEvent.h"

ClassImp(AliHFEminiEvent)

//_______________________________________
AliHFEminiEvent::AliHFEminiEvent():
  TObject(),
  fTracks(NULL),
  fMCparticles(NULL),
  fNtracks(0),
  fNmcparticles(0),
  fVZ(0.),
  fSPDMultiplicity(0),
  fCentralityBin(0),
  fNprimaryNchMC(0)
  
{
  //
  // Default constructor
  //
  fTracks = new TObjArray;
  fTracks->SetOwner();
  fMCparticles = new TObjArray;
  fMCparticles->SetOwner();
  memset(fV0Multiplicity, 0, sizeof(Float_t) * 2);
}

//_______________________________________
AliHFEminiEvent::AliHFEminiEvent(const AliHFEminiEvent &ref):
  TObject(ref),
  fTracks(NULL),
  fMCparticles(NULL),
  fNtracks(ref.fNtracks),
  fNmcparticles(ref.fNmcparticles),
  fVZ(ref.fVZ),
  fSPDMultiplicity(ref.fSPDMultiplicity),
  fCentralityBin(ref.fCentralityBin),
  fNprimaryNchMC(ref.fNprimaryNchMC)
  
{
  //
  // Copy constructor
  //
  fTracks = new TObjArray;
  fTracks->SetOwner();
  for(int itrk = 0; itrk < ref.GetNumberOfTracks(); itrk++)
    fTracks->Add(new AliHFEminiTrack(*(ref.GetTrack(itrk))));
  fMCparticles = new TObjArray;
  fMCparticles->SetOwner();
  for(int iprt = 0; iprt < ref.GetNumberOfMCParticles(); iprt++)
    fMCparticles->Add(new AliHFEreducedMCParticle(*(ref.GetMCParticle(iprt))));
  memcpy(fV0Multiplicity, ref.fV0Multiplicity, sizeof(Float_t) * 2);
}

//_______________________________________
AliHFEminiEvent &AliHFEminiEvent::operator=(const AliHFEminiEvent &ref){
  //
  // Assignment operator
  //
  if(&ref != this){
    TObject::operator=(ref);
    fTracks->SetOwner();
    fTracks->Clear();
    for(int itrk = 0; itrk < ref.GetNumberOfTracks(); itrk++)
      fTracks->Add(new AliHFEminiTrack(*(ref.GetTrack(itrk))));
    fMCparticles->Clear();
    fMCparticles->SetOwner();
    for(int iprt = 0; iprt < ref.GetNumberOfMCParticles(); iprt++)
      fMCparticles->Add(new AliHFEreducedMCParticle(*(ref.GetMCParticle(iprt))));
    fNtracks = ref.fNtracks;
    fNmcparticles = ref.fNmcparticles;
    fVZ = ref.fVZ;
    fCentralityBin = ref.fCentralityBin;
    fSPDMultiplicity = ref.fSPDMultiplicity;
    fNprimaryNchMC = ref.fNprimaryNchMC;
    memcpy(fV0Multiplicity, ref.fV0Multiplicity, sizeof(Float_t) * 2);
  }
  return *this;
}

//_______________________________________
AliHFEminiEvent::~AliHFEminiEvent(){
  //
  // Destructor: Clear tracks an MC particles
  //
  delete fTracks;
  delete fMCparticles;
}

//_______________________________________
void AliHFEminiEvent::AddTrack(const AliHFEminiTrack *track){
  //
  // Add track to the event
  //
  fTracks->Add(new AliHFEminiTrack(*track));
  fNtracks++;
}

//_______________________________________
const AliHFEminiTrack *AliHFEminiEvent::GetTrack(Int_t itrk) const {
  //
  // Get Track
  //
  if(itrk < 0 || itrk >= fNtracks) return NULL;
  return dynamic_cast<const AliHFEminiTrack *>(fTracks->At(itrk));
}

//_______________________________________
void AliHFEminiEvent::AddMCParticle(const AliHFEreducedMCParticle *track){
  //
  // Add MC particle to the Event
  //
  fMCparticles->Add(new AliHFEreducedMCParticle(*track));
  fNmcparticles++;
}

//_______________________________________
const AliHFEreducedMCParticle *AliHFEminiEvent::GetMCParticle(Int_t itrk) const {
  //
  // Get MC particle
  //
  if(itrk < 0 || itrk >= fNmcparticles) return NULL;
  return dynamic_cast<const AliHFEreducedMCParticle *>(fMCparticles->At(itrk));
}
