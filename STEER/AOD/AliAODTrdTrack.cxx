#include "AliAODTrdTrack.h"
AliAODTrdTrack::AliAODTrdTrack() :
  AliVTrdTrack(),
  fGlobalStack(-1),
  fPID(0),
  fLayerMask(0),
  fA(0),
  fFlagsTiming(0),
  fTracklets(),
  fTrackMatch(0x0),
  fLabel(-1)
{
  /// default constructor

  fTracklets.SetClass("AliAODTrdTracklet", 6);
}

AliAODTrdTrack::AliAODTrdTrack(const AliVTrdTrack &rhs) :
  AliVTrdTrack(rhs),
  fGlobalStack(5*rhs.GetSector() + rhs.GetStack()),
  fPID(rhs.GetPID()),
  fLayerMask(rhs.GetLayerMask()),
  fA(rhs.GetA()),
  fFlagsTiming(rhs.GetFlagsTiming()),
  fTracklets(),
  fTrackMatch(rhs.GetTrackMatch()),
  fLabel(rhs.GetLabel())
{
  /// constructor from abstract base class

  fTracklets.SetClass("AliAODTrdTracklet", 6);

  // copy the contributing tracklets
  for (Int_t iTracklet = 0; iTracklet < 6; ++iTracklet) {
    const AliVTrdTracklet *trkl = rhs.GetTracklet(iTracklet);
    if (trkl)
      new (fTracklets[iTracklet]) AliAODTrdTracklet(*trkl);
    else
      new (fTracklets[iTracklet]) AliAODTrdTracklet();
  }
}

AliAODTrdTrack::AliAODTrdTrack(const AliAODTrdTrack& rhs) :
  AliVTrdTrack(rhs),
  fGlobalStack(rhs.fGlobalStack),
  fPID(rhs.fPID),
  fLayerMask(rhs.fLayerMask),
  fA(rhs.fA),
  fFlagsTiming(rhs.fFlagsTiming),
  fTracklets(),
  fTrackMatch(rhs.fTrackMatch),
  fLabel(rhs.fLabel)
{
  /// copy constructor

  fTracklets.SetClass("AliAODTrdTracklet", 6);

  // copy the contributing tracklets
  for (Int_t iTracklet = 0; iTracklet < 6; ++iTracklet) {
    const AliVTrdTracklet *trkl = rhs.GetTracklet(iTracklet);
    if (trkl)
      new (fTracklets[iTracklet]) AliAODTrdTracklet(*trkl);
    else
      new (fTracklets[iTracklet]) AliAODTrdTracklet();
  }
}

AliAODTrdTrack& AliAODTrdTrack::operator=(const AliAODTrdTrack& rhs)
{
  /// assignment operator

  if (&rhs != this)
    AliVTrdTrack::operator=(rhs);

  fGlobalStack = rhs.fGlobalStack;
  fPID         = rhs.fPID;
  fLayerMask   = rhs.fLayerMask;
  fA           = rhs.fA;
  fFlagsTiming = rhs.fFlagsTiming;
  fTrackMatch  = rhs.fTrackMatch;
  fLabel       = rhs.fLabel;

  // assign the contributing tracklets
  for (Int_t iTracklet = 0; iTracklet < 6; ++iTracklet) {
    const AliVTrdTracklet *trkl = rhs.GetTracklet(iTracklet);
    if (trkl)
      new (fTracklets[iTracklet]) AliAODTrdTracklet(*trkl);
    else
      new (fTracklets[iTracklet]) AliAODTrdTracklet();
  }

  return *this;
}

void AliAODTrdTrack::Copy(TObject &rhs) const
{
  /// copy

  AliVTrdTrack::Copy(rhs);
}

Int_t AliAODTrdTrack::GetPt() const
{
  /// calculate pt from a as done in hardware

  const Int_t maskIdLut[64] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,
    -1, -1, -1, -1, -1, -1, -1,  1, -1, -1, -1,  2, -1,  3,  4,  5,
    -1, -1, -1, -1, -1, -1, -1,  6, -1, -1, -1,  7, -1,  8,  9, 10,
    -1, -1, -1, 11, -1, 12, 13, 14, -1, 15, 16, 17, 18, 19, 20, 21
  };

  const Int_t c1Lut[32] = {
    -2371, -2474, -2474, -2474, -2563, -2448, -2578, -2578,
    -2578, -2670, -2557, -2578, -2578, -2670, -2557, -2578,
    -2670, -2557, -2763, -2557, -2644, -2523,    -1,    -1,
    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1
  };

  if (this->GetA() != 0) {
    Int_t layerMaskId = maskIdLut[this->GetLayerMask()];
    Int_t c1 = c1Lut[layerMaskId];
    Int_t c1Ext = c1 << 8;
    Int_t ptRawStage4 = c1Ext / ((this->GetA() >> 2) != 0 ? (this->GetA() >> 2) : 1 );
    Int_t ptRawComb4 = ptRawStage4;
    Int_t ptExtComb4 = (ptRawComb4 > 0) ? ptRawComb4 + 33 : ptRawComb4 - 30;

    return -ptExtComb4/2;
  }
  else
    return 0;
}
