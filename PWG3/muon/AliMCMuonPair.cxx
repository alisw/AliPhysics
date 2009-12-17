#include "AliAODMuonPair.h"
#include "AliMCMuonTrack.h"
#include "AliMCMuonPair.h"

ClassImp(AliAODMuonPair)

//-----------------------------------------------------------------------------
AliMCMuonPair::AliMCMuonPair() :
AliAODMuonPair(),
fIsFull(kFALSE),
fPGen(),
fSource(-1)
{
  //
  // default constructor
  //
  for (Int_t i=0; i<2; i++) fTrk[i] = 0;

}

//-----------------------------------------------------------------------------
AliMCMuonPair::AliMCMuonPair(AliMCMuonTrack *trk0, AliMCMuonTrack *trk1, Bool_t full) :
AliAODMuonPair(),
fIsFull(full),
fPGen(),
fSource(-1)
{
  //
  // default constructor
  //
  fTrk[0] = trk0;
  fTrk[1] = trk1;
  fPGen = trk0->GetPGen() + trk1->GetPGen();
  AliAODMuonPair::FillPairInfo();
  if (fIsFull) this->FindDimuonSourceFull();
  else this->FindDimuonSourceFast();
}

//-----------------------------------------------------------------------------
AliMCMuonPair::~AliMCMuonPair()
{
  //
  // destructor
  //
}

//-----------------------------------------------------------------------------
void AliMCMuonPair::FindDimuonSourceFast()
{
  AliMCMuonTrack *trk0 = (AliMCMuonTrack*)fTrk[0].GetObject();
  Int_t src0 = trk0->GetSource();
  if (src0<0 || src0==4 || src0==3) {
    fSource=5; return;
  }

  AliMCMuonTrack *trk1 = (AliMCMuonTrack*)fTrk[1].GetObject();
  Int_t src1 = trk1->GetSource();
  if (src1<0 || src1==4 || src1==3) {
    fSource=5; return;
  }

  // Drell-Yan is expected very small at LHC, we ingore it
  Int_t np0 = trk0->GetNParents() - 1;
  if (np0<0) {
    fSource=5; return;
  }

  Int_t np1 = trk1->GetNParents() - 1;
  if (np1<0) {
    fSource=5; return;
  }

  if (trk0->IsMotherAResonance(np0) && trk1->IsMotherAResonance(np1) &&
      (trk0->GetParentIndex(np0))==(trk1->GetParentIndex(np1))) {
    fSource=4; return;
  }

  if (src0==0 && src1==0) {
    if ((trk0->GetParentIndex(0))==(trk1->GetParentIndex(0)))
      fSource = 1;
    else
      fSource = 0;
    return;
  }

  if (src0==1 && src1==1) {
    if ((trk0->GetParentIndex(0))==(trk1->GetParentIndex(0)))
      fSource = 3;
    else
      fSource = 2;
    return;
  }

  fSource = 5;
  return;
}
