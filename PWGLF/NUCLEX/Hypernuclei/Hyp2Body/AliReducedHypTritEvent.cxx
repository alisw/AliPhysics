
#include "AliReducedHypTritEvent.h"
#include "TLorentzVector.h"
#include <TClonesArray.h>

ClassImp(AliReducedHypTritTrack)
ClassImp(AliReducedHypTritV0)
ClassImp(AliReducedHypTritEvent)

TClonesArray *AliReducedHypTritEvent::fgV0s = 0;

AliReducedHypTritTrack::AliReducedHypTritTrack() :
  TObject(),
  fP(),
  fDca(0),
  fDedx(0),
  fDedxSigma(0),
  fCharge(0),
  fEta(0),
  fPhi(0),
  fTpcNClusters(0) {
}

AliReducedHypTritTrack::~AliReducedHypTritTrack() {
}

AliReducedHypTritV0::AliReducedHypTritV0() :
  TObject(),
  fPosition(),
  fP(0),
  fPt(0),
  fM(0),
  fDCAv0(0),
  fSigmaD0(0),
  fCosPointingAngle(0),
  fDecayLength(0),
  fChi2(0),
  fMCTruth(0),
  fRapidity(-1),
  fCharge(0),
  fPiTrack(0x0),
  fHeTrack(0x0) {
  if (!fPiTrack) fPiTrack = new AliReducedHypTritTrack();
  if (!fHeTrack) fHeTrack = new AliReducedHypTritTrack();
}

AliReducedHypTritV0::~AliReducedHypTritV0() {
}

AliReducedHypTritEvent::AliReducedHypTritEvent() :
  TObject(),
  fEventId(0),
  fVertexPosition(),
  fNumberV0s(0),
  fCentrality(0),
  fRunNumber(0),
  fV0s(0x0),
  fTrigger() {
  if (!fgV0s) fgV0s = new TClonesArray("AliReducedHypTritV0", 100000);
  fV0s = fgV0s;
}

AliReducedHypTritEvent::~AliReducedHypTritEvent() {
}

void AliReducedHypTritEvent::ClearEvent() {
  if (fV0s) fV0s->Clear("C");
  fNumberV0s = 0;
  fCentrality = 0;
}
