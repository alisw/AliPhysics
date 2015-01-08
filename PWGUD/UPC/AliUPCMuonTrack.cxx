
//_____________________________________________________________________________
//    Class for UPC muon track
//    Author: Jaroslav Adam
//
//    contains parameters of the muon track relevant for UPC analysis
//_____________________________________________________________________________

#include "TLorentzVector.h"

#include "AliUPCMuonTrack.h"

ClassImp(AliUPCMuonTrack)

//_____________________________________________________________________________
AliUPCMuonTrack::AliUPCMuonTrack()
 :TObject(),
  fPt(0), fEta(0), fPhi(0), fCharge(0), fMatchTrigger(0), fRabs(0),
  fChi2perNDF(0), fDca(0), fPdca(0),
  fkMuonMass(0.105658)
{

  // Default constructor


}

void AliUPCMuonTrack::GetMomentum(TLorentzVector *v) const
{
  // get track 4-momentum
  v->SetPtEtaPhiM(fPt,fEta,fPhi,fkMuonMass);
}












































