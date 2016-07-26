
#include "AliReducedHypTritEvent.h"
#include "TLorentzVector.h"
#include <TClonesArray.h>

ClassImp(AliReducedHypTritTrack)
ClassImp(AliReducedHypTritV0)
ClassImp(AliReducedHypTritEvent)

TClonesArray *AliReducedHypTritEvent::fgV0s = 0;

AliReducedHypTritTrack::AliReducedHypTritTrack() :
  TObject(),
  fPt(0),
  fMomentum(),
  fDCAtoPrim(0),
  fP(0),
  fDedx(0),
  fSign(0),
  fEta(0),
  fPhi(0) {

}

AliReducedHypTritTrack::~AliReducedHypTritTrack() {

}

AliReducedHypTritV0::AliReducedHypTritV0() :
  TObject(),
  fSecVertexPos(),
  fDCAtoPrim(0),
  fMotherP(0),
  fMotherPt(0),
  fMotherInvMass(0),
  fDCAv0(0),
  fSigmaD0(0),
  fCosPointingAngle(0),
  fDecayRadius(0),
  fOnFlyStatus(0),
  fChi2(0),
  fMCTruth(0),
  fRapidity(-1),
  fAntiParticle(0),
  fPosTrack(0x0),
  fNegTrack(0x0) {
  if (!fPosTrack) fPosTrack = new AliReducedHypTritTrack();
  if (!fNegTrack) fNegTrack = new AliReducedHypTritTrack();
}

AliReducedHypTritV0::~AliReducedHypTritV0() {

}

AliReducedHypTritEvent::AliReducedHypTritEvent() :
  TObject(),
  fPrimVertexPos(),
  fNV0s(0),
  fCentrality(0),
  fRunNumber(0),
  fV0s(0x0),
  fTrigger() {
  if (!fgV0s) fgV0s = new TClonesArray("AliReducedHypTritV0", 100000);
  fV0s = fgV0s;
  for (Int_t i = 0; i < 3; i++) {
    fTrigger[i] = kFALSE;
  }

}


AliReducedHypTritEvent::~AliReducedHypTritEvent() {

}

void AliReducedHypTritEvent::ClearEvent() {
  if (fV0s) fV0s->Clear("C");
  fNV0s = 0;
  fCentrality = 0;
  for (Int_t i = 0; i < 3; i++) {
    fTrigger[i] = kFALSE;
  }
}
