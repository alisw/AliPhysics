#include "AliDielectronReducedTrack.h"
#include "AliVParticle.h"

AliDielectronReducedTrack::AliDielectronReducedTrack() :
fPx(0),
fPy(0),
fPz(0),
fCharge(0),
fIsTagged(0)
{}


AliDielectronReducedTrack::AliDielectronReducedTrack(Double_t px, Double_t py, Double_t pz, Short_t charge,Bool_t IsTagged ):
fPx(px),
fPy(py),
fPz(pz),
fCharge(charge),
fIsTagged(IsTagged)
{}

AliDielectronReducedTrack::~AliDielectronReducedTrack() {}
