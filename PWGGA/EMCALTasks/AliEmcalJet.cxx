// $Id$

#include "AliEmcalJet.h"

ClassImp(AliEmcalJet)

//__________________________________________________________________________________________________
AliEmcalJet::AliEmcalJet(Double_t px, Double_t py, Double_t pz) 
  : AliVParticle(), 
    fPt(TMath::Sqrt(px*px+py*py)), 
    fEta(TMath::ASinH(pz/fPt)),
    fPhi(0), fM(0), fNEF(0), 
    fArea(0), fNch(0), fNn(0),
    fMaxCPt(0), fMaxNPt(0)
{    
  // Constructor.

  if (fPt != 0) {
    fPhi = TMath::ATan2(py, px);
    if (fPhi<0.) 
      fPhi += 2. * TMath::Pi();
  }
}

//_________________________________________________________________________________________________
AliEmcalJet::AliEmcalJet(Double_t pt, Double_t eta, Double_t phi, Double_t m) :
  AliVParticle(), 
  fPt(pt), fEta(eta), fPhi(phi), 
  fM(m), fNEF(0), fArea(0), 
  fNch(0), fNn(0), 
  fMaxCPt(0), fMaxNPt(0) 
{
  // Constructor.

 if (fPhi<0.) 
   fPhi += 2. * TMath::Pi();
}

//_________________________________________________________________________________________________
AliEmcalJet::AliEmcalJet(const AliEmcalJet &jet) :
  AliVParticle(jet),
  fPt(jet.fPt), fEta(jet.fEta), fPhi(jet.fPhi), 
  fM(jet.fM), fNEF(jet.fNEF), fArea(jet.fArea), 
  fNch(jet.fNch), fNn(jet.fNn),
  fMaxCPt(jet.fMaxCPt), fMaxNPt(jet.fMaxNPt)
{
  // Constructor.
}

//_________________________________________________________________________________________________
AliEmcalJet &AliEmcalJet::operator=(const AliEmcalJet &jet)
{
  // Assignment operator.

  if (this!=&jet) {
    AliVParticle::operator=(jet);
    fPt     = jet.fPt;
    fEta    = jet.fEta;
    fPhi    = jet.fPhi;
    fM      = jet.fM; 
    fNEF    = jet.fNEF;
    fArea   = jet.fArea; 
    fNch    = jet.fNch; 
    fNn     = jet.fNn;
    fMaxCPt = jet.fMaxCPt; 
    fMaxNPt = jet.fMaxNPt;
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
