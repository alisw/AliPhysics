// $Id$
//
// Emcal jet class.
//
// Author: C.Loizides

#include "AliEmcalJet.h"

ClassImp(AliEmcalJet)

//__________________________________________________________________________________________________
AliEmcalJet::AliEmcalJet() : 
  AliVParticle(), 
  fPt(0), 
  fEta(0), 
  fPhi(0), 
  fM(0), 
  fNEF(0),
  fArea(0),       
  fAreaEmc(-1), 
  fAxisInEmcal(0), 
  fMaxCPt(0), 
  fMaxNPt(0), 
  fMCPt(0),
  fNn(0), 
  fNch(0),        
  fPtEmc(0),
  fNEmc(0),
  fClusterIDs(),
  fTrackIDs(),
  fMatched(2),
  fPtSub(0)
{
  // Constructor.

  fClosestJets[0] = 0;
  fClosestJets[1] = 0; 
  fClosestJetsDist[0] = 999; 
  fClosestJetsDist[1] = 999; 
}

//__________________________________________________________________________________________________
AliEmcalJet::AliEmcalJet(Double_t px, Double_t py, Double_t pz) : 
  AliVParticle(), 
  fPt(TMath::Sqrt(px*px+py*py)), 
  fEta(TMath::ASinH(pz/fPt)),
  fPhi(0), 
  fM(0), 
  fNEF(0), 
  fArea(0), 
  fAreaEmc(-1), 
  fAxisInEmcal(0),
  fMaxCPt(0), 
  fMaxNPt(0), 
  fMCPt(0),
  fNn(0),
  fNch(0),
  fPtEmc(0),
  fNEmc(0),
  fClusterIDs(), 
  fTrackIDs(),
  fMatched(2),
  fPtSub(0)
{    
  // Constructor.

  if (fPt != 0) {
    fPhi = TMath::ATan2(py, px);
    if (fPhi<0.) 
      fPhi += 2. * TMath::Pi();
  }

  fClosestJets[0] = 0; 
  fClosestJets[1] = 0;
  fClosestJetsDist[0] = 999; 
  fClosestJetsDist[1] = 999;
}

//_________________________________________________________________________________________________
AliEmcalJet::AliEmcalJet(Double_t pt, Double_t eta, Double_t phi, Double_t m) :
  AliVParticle(), 
  fPt(pt), 
  fEta(eta), 
  fPhi(phi), 
  fM(m), 
  fNEF(0), 
  fArea(0), 
  fAreaEmc(-1), 
  fAxisInEmcal(0),
  fMaxCPt(0), 
  fMaxNPt(0),
  fMCPt(0),
  fNn(0),
  fNch(0), 
  fPtEmc(0),
  fNEmc(0),
  fClusterIDs(), 
  fTrackIDs(),
  fMatched(2),
  fPtSub(0)
{
  // Constructor.

  if (fPhi<0.) 
    fPhi += TMath::TwoPi();

  fClosestJets[0] = 0; 
  fClosestJets[1] = 0;
  fClosestJetsDist[0] = 999; 
  fClosestJetsDist[1] = 999;
}

//_________________________________________________________________________________________________
AliEmcalJet::AliEmcalJet(const AliEmcalJet &jet) :
  AliVParticle(jet),
  fPt(jet.fPt), 
  fEta(jet.fEta), 
  fPhi(jet.fPhi), 
  fM(jet.fM), 
  fNEF(jet.fNEF), 
  fArea(jet.fArea), 
  fAreaEmc(jet.fAreaEmc), 
  fAxisInEmcal(jet.fAxisInEmcal),
  fMaxCPt(jet.fMaxCPt), 
  fMaxNPt(jet.fMaxNPt), 
  fMCPt(jet.fMCPt),
  fNn(jet.fNn),
  fNch(jet.fNch),
  fPtEmc(jet.fPtEmc),
  fNEmc(jet.fNEmc),
  fClusterIDs(jet.fClusterIDs), 
  fTrackIDs(jet.fTrackIDs),
  fMatched(jet.fMatched),
  fPtSub(jet.fPtSub)
{
  // Copy constructor.

  fClosestJets[0]     = jet.fClosestJets[0]; 
  fClosestJets[1]     = jet.fClosestJets[1]; 
  fClosestJetsDist[0] = jet.fClosestJetsDist[0];  
  fClosestJetsDist[1] = jet.fClosestJetsDist[1]; 
}

//_________________________________________________________________________________________________
AliEmcalJet &AliEmcalJet::operator=(const AliEmcalJet &jet)
{
  // Assignment operator.

  if (this!=&jet) {
    AliVParticle::operator=(jet);
    fPt                 = jet.fPt;
    fEta                = jet.fEta;
    fPhi                = jet.fPhi;
    fM                  = jet.fM; 
    fNEF                = jet.fNEF;
    fArea               = jet.fArea; 
    fAreaEmc            = jet.fAreaEmc; 
    fAxisInEmcal        = jet.fAxisInEmcal; 
    fMaxCPt             = jet.fMaxCPt; 
    fMaxNPt             = jet.fMaxNPt;
    fMCPt               = jet.fMCPt;
    fNn                 = jet.fNn;
    fNch                = jet.fNch;
    fPtEmc              = jet.fPtEmc;
    fNEmc               = jet.fNEmc;
    fClusterIDs         = jet.fClusterIDs;
    fTrackIDs           = jet.fTrackIDs;
    fClosestJets[0]     = jet.fClosestJets[0]; 
    fClosestJets[1]     = jet.fClosestJets[1]; 
    fClosestJetsDist[0] = jet.fClosestJetsDist[0];  
    fClosestJetsDist[1] = jet.fClosestJetsDist[1]; 
    fMatched            = jet.fMatched;
    fPtSub              = jet.fPtSub;
  }

  return *this;
}

//__________________________________________________________________________________________________
void AliEmcalJet::GetMom(TLorentzVector &vec) const
{
  // Return momentum as four vector.

  Double_t p = fPt *TMath::CosH(fEta);
  vec.SetPtEtaPhiE(fPt,fEta,fPhi,TMath::Sqrt(p*p+fM*fM));
}

//__________________________________________________________________________________________________
void AliEmcalJet::Print(Option_t* /*option*/) const
{
  // Print jet information.

  printf("Jet pt=%.2f, eta=%.2f, phi=%.2f, area=%.2f, NEF=%.2f\n", fPt, fEta, fPhi, fArea, fNEF);
}

//__________________________________________________________________________________________________
void AliEmcalJet::SortConstituents()
{
  // Sort constituent by index (increasing).

  std::sort(fClusterIDs.GetArray(), fClusterIDs.GetArray() + fClusterIDs.GetSize());
  std::sort(fTrackIDs.GetArray(), fTrackIDs.GetArray() + fTrackIDs.GetSize());
}
