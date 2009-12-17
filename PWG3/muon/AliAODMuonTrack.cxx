#include <TObject.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

#include "AliAODTrack.h"
#include "AliESDMuonTrack.h"
#include "AliMCMuonTrack.h"
#include "AliAODMuonTrack.h"

ClassImp(AliAODMuonTrack)

//-----------------------------------------------------------------------------
AliAODMuonTrack::AliAODMuonTrack() :
TObject(),
fP(),
fCharge(0),
fTrigger(0),
fDca(0.),
fChi2(0.),
fCentr(0.)
{
  //
  // default constructor
  //
}

//-----------------------------------------------------------------------------
AliAODMuonTrack::AliAODMuonTrack(AliAODTrack *trk) :
TObject(),
fP(),
fCharge(0),
fTrigger(0),
fDca(0.),
fChi2(0.),
fCentr(0.)
{
  //
  // default constructor
  //
 this->FillTrackInfo(trk);
}

//-----------------------------------------------------------------------------
AliAODMuonTrack::AliAODMuonTrack(AliESDMuonTrack *trk) :
TObject(),
fP(),
fCharge(0),
fTrigger(0),
fDca(0.),
fChi2(0.),
fCentr(0.)
{
  //
  // default constructor
  //
 this->FillTrackInfo(trk);
}

//-----------------------------------------------------------------------------
AliAODMuonTrack::~AliAODMuonTrack()
{
  //
  // destructor
  //
}

//-----------------------------------------------------------------------------
void AliAODMuonTrack::FillTrackInfo(AliAODTrack *trk)
{
  Double_t mMu    = TDatabasePDG::Instance()->GetParticle(13)->Mass();
  Double_t px     = trk->Px();
  Double_t py     = trk->Py();
  Double_t pz     = trk->Pz();
  Double_t energy = TMath::Sqrt(mMu*mMu + px*px + py*py + pz*pz);
  fP.SetPxPyPzE(px,py,pz,energy);
  fCharge = trk->Charge();
  fTrigger = trk->GetMatchTrigger();
  fDca = trk->DCA();
  fChi2 = trk->Chi2perNDF();
  return;
}

//-----------------------------------------------------------------------------
void AliAODMuonTrack::FillTrackInfo(AliESDMuonTrack *trk)
{
  trk->LorentzP(fP);
  fCharge = trk->Charge();
  fTrigger = trk->GetMatchTrigger();
  fDca = trk->GetDCA();
  fChi2 = trk->GetChi2()/(2.*trk->GetNHit()-5.);
  return;
}

//-----------------------------------------------------------------------------
Bool_t AliAODMuonTrack::SelectSingleMuon(Double_t cuts[10])
{
  TLorentzVector lorentzP = this->GetP();
  Double_t p = lorentzP.P();
  if (p<cuts[0] || (p>cuts[1])) return kFALSE;

  Double_t pt = lorentzP.Pt();
  if (pt<cuts[2] || pt>cuts[3]) return kFALSE;

  Double_t eta = lorentzP.Eta();
  if (eta<cuts[4] || eta>cuts[5]) return kFALSE;

  Double_t dca = this->GetDCA();
  if (dca<cuts[6] || dca>cuts[7]) return kFALSE;

  Int_t trigger = this->GetTrigger();
  if (trigger<cuts[8] || trigger>cuts[9]) return kFALSE;

  return kTRUE;
}
