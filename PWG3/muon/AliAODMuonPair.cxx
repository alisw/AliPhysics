#include "AliAODMuonTrack.h"
#include "AliAODMuonPair.h"

ClassImp(AliAODMuonPair)

//-----------------------------------------------------------------------------
AliAODMuonPair::AliAODMuonPair() :
TObject(),
fP(),
fCharge(0)
{
  //
  // default constructor
  //
  for (Int_t i=0; i<2; i++) {
    fTrk[i] = 0;
    fTrigger[i] = 0;
    fDca[i] = 0.;
    fChi2[i] = 0.;
  }
}

//-----------------------------------------------------------------------------
AliAODMuonPair::AliAODMuonPair(AliAODMuonTrack *trk0, AliAODMuonTrack *trk1) :
TObject(),
fP(),
fCharge(0)
{
  //
  // default constructor
  //
  for (Int_t i=0; i<2; i++) {
    fTrigger[i]=0; fDca[i]=0.; fChi2[i]=0.;
  }

  fTrk[0] = trk0;
  fTrk[1] = trk1;
  FillPairInfo();
}

//-----------------------------------------------------------------------------
AliAODMuonPair::~AliAODMuonPair()
{
  //
  // destructor
  //
}

//-----------------------------------------------------------------------------
void AliAODMuonPair::FillPairInfo()
{
  AliAODMuonTrack *trk0 = (AliAODMuonTrack*)fTrk[0].GetObject();
  AliAODMuonTrack *trk1 = (AliAODMuonTrack*)fTrk[1].GetObject();

  fP = trk0->GetP() + trk1->GetP();
  fCharge = trk0->GetCharge() + trk1->GetCharge();

  fTrigger[0] = trk0->GetTrigger();
  fTrigger[1] = trk1->GetTrigger();

  fDca[0] = trk0->GetDCA();
  fDca[1] = trk1->GetDCA();

  fChi2[0] = trk0->GetChi2();
  fChi2[1] = trk1->GetChi2();

  return;
}
