
#include "AliReducedHypTritEvent.h"
#include "TLorentzVector.h"
#include <TClonesArray.h>

ClassImp(AliReducedHypTritTrack)
ClassImp(AliReducedHypTritV0)
ClassImp(AliReducedHypTritEvent)

AliReducedHypTritTrack::AliReducedHypTritTrack() :
  TObject(),
  fP(),
  fPtrack(-1),
  fDca(-1),
  fDedx(-1),
  fDedxSigma(-9999),
  fDedxSigmaTriton(-9999),
  fEta(-999),
  fPhi(-999),
  fTpcNClusters(-1),
  fGeoLength(-1) {
}

AliReducedHypTritTrack::~AliReducedHypTritTrack() {
}

AliReducedHypTritV0::AliReducedHypTritV0() :
  TObject(),
  fPosition(),
  fPvect(),
  fPiTrack(0x0),
  fHeTrack(0x0),
  fP(-1),
  fPt(-1),
  fM(-1),
  fDcaV0(-999),
  fCosPointingAngle(0),
  fDecayLength(-1),
  fMcTruth(0),
  fRapidity(-999),
  fParticleSpecies(0),
  fCharge(-999),
  fOnFlyStatus() {
  if (!fPiTrack) fPiTrack = new AliReducedHypTritTrack();
  if (!fHeTrack) fHeTrack = new AliReducedHypTritTrack();
}

AliReducedHypTritV0::~AliReducedHypTritV0() {
}

AliReducedHypTritEvent::AliReducedHypTritEvent() :
  TObject(),
  fVertexPosition(),
  fV0s(0x0),
  fNumberV0s(0),
  fCentrality(-1),
  fRunNumber(0),
  fTrigger(),
  fTriggerClasses() {
  if (!fV0s) fV0s = new TClonesArray("AliReducedHypTritV0", 100000);
}

AliReducedHypTritEvent::~AliReducedHypTritEvent() {
}

void AliReducedHypTritEvent::ClearEvent() {
  if (fV0s) fV0s->Clear("C");
  fNumberV0s = 0;
}
