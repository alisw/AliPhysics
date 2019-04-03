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
  fSPDMultiplicity(0),
  fPileupFlag(kFALSE),
  fV0PlanePhi(0.),
  fV0APlanePhi(0.),
  fV0CPlanePhi(0.),
  fTPCPlanePhi(0.),
  fV0PlanePhiCorrected(0.),
  fV0APlanePhiCorrected(0.),
  fV0CPlanePhiCorrected(0.),
  fMagneticField(0.)
{
  //
  // Default constructor
  //
  fTracks = new TObjArray;
  fTracks->SetOwner();
  fMCparticles = new TObjArray;
  fMCparticles->SetOwner();
  memset(fCentrality, 0, sizeof(Float_t) * kCentBuff);
  memset(fV0Multiplicity, 0, sizeof(Float_t) * 2);
  memset(fZDCEnergy, 0, sizeof(Float_t) * 4);
  memset(fVX, 0, sizeof(Float_t)*2);
  memset(fVY, 0, sizeof(Float_t)*2);
  memset(fVZ, 0, sizeof(Float_t)*2);
  memset(fVMC, 0, sizeof(Double_t)*3);
  memset(fNContrib, 0, sizeof(Int_t) * 2);
  fVertexResolution[1] = 999.;
  fVertexResolution[0] = fVertexResolution[1];
  fVertexDispersion[0] = 999.;
  fVertexDispersion[1] = fVertexDispersion[0];
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
  fSPDMultiplicity(ref.fSPDMultiplicity),
  fPileupFlag(ref.fPileupFlag),
  fV0PlanePhi(ref.fV0PlanePhi),
  fV0APlanePhi(ref.fV0APlanePhi),
  fV0CPlanePhi(ref.fV0CPlanePhi),
  fTPCPlanePhi(ref.fTPCPlanePhi),
  fMagneticField(ref.fMagneticField)
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
  memcpy(fCentrality, ref.fCentrality, sizeof(Float_t) * kCentBuff);
  memcpy(fV0Multiplicity, ref.fV0Multiplicity, sizeof(Float_t) * 2);
  memcpy(fZDCEnergy, ref.fZDCEnergy, sizeof(Float_t) *4);
  memcpy(fVX, ref.fVX, sizeof(Float_t)*2);
  memcpy(fVY, ref.fVY, sizeof(Float_t)*2);
  memcpy(fVZ, ref.fVZ, sizeof(Float_t)*2);
  memcpy(fVMC, ref.fVMC, sizeof(Double_t)*3);
  memcpy(fVertexResolution, ref.fVertexResolution, sizeof(Float_t)*2);
  memcpy(fVertexDispersion, ref.fVertexDispersion, sizeof(Float_t)*2);
  memcpy(fNContrib, ref.fNContrib, sizeof(Int_t) * 2);
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
    memcpy(fVX, ref.fVX, sizeof(Float_t)*2);
    memcpy(fVY, ref.fVY, sizeof(Float_t)*2);
    memcpy(fVZ, ref.fVZ, sizeof(Float_t)*2);
    memcpy(fVMC, ref.fVMC, sizeof(Double_t)*3);
    memcpy(fNContrib, ref.fNContrib, sizeof(Int_t) * 2);
    memcpy(fVertexResolution, ref.fVertexResolution, sizeof(Float_t)*2);
    memcpy(fVertexDispersion, ref.fVertexDispersion, sizeof(Float_t)*2);
    fSPDMultiplicity = ref.fSPDMultiplicity;
    fPileupFlag = ref.fPileupFlag;
    memcpy(fCentrality, ref.fCentrality, sizeof(Float_t) * kCentBuff);
    memcpy(fV0Multiplicity, ref.fV0Multiplicity, sizeof(Float_t) * 2);
    memcpy(fZDCEnergy, ref.fZDCEnergy, sizeof(Float_t) *4);
    fMagneticField = ref.fMagneticField;
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
