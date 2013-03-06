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
  fAreaEta(0),       
  fAreaPhi(0),       
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
  fMatchingType(0),
  fPtSub(0),
  fPtVectSub(0)
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
  fAreaEta(0),       
  fAreaPhi(0),       
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
  fMatchingType(0),
  fPtSub(0),
  fPtVectSub(0)
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
  fAreaEta(0),       
  fAreaPhi(0),       
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
  fMatchingType(0),
  fPtSub(0),
  fPtVectSub(0)
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
  fAreaEta(jet.fAreaEta),       
  fAreaPhi(jet.fAreaPhi),       
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
  fMatchingType(jet.fMatchingType),
  fPtSub(jet.fPtSub),
  fPtVectSub(jet.fPtVectSub)
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
    fAreaEta            = jet.fAreaEta; 
    fAreaPhi            = jet.fAreaPhi; 
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
    fPtVectSub          = jet.fPtVectSub;
  }

  return *this;
}

//_________________________________________________________________________________________________
Int_t AliEmcalJet::Compare(const TObject* obj) const
{
  //Return -1 if this is smaller than obj, 0 if objects are equal and 1 if this is larger than obj.

  const AliEmcalJet *jet = static_cast<const AliEmcalJet *>(obj);
  if (!obj)
    return 0;
  if (Pt()>jet->Pt())
    return -1;
  return 1;
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
Double_t AliEmcalJet::PtSubVect(Double_t rho) const
{
  // Return vectorial subtracted transverse momentum.

  Double_t dx = Px() - rho * fArea * TMath::Cos(fAreaPhi);
  Double_t dy = Py() - rho * fArea * TMath::Sin(fAreaPhi);
  //Double_t dz = Pz() - rho * fArea * TMath::SinH(fAreaEta);
  return TMath::Sqrt(dx*dx+dy*dy);
}

//__________________________________________________________________________________________________
void AliEmcalJet::SortConstituents()
{
  // Sort constituent by index (increasing).

  std::sort(fClusterIDs.GetArray(), fClusterIDs.GetArray() + fClusterIDs.GetSize());
  std::sort(fTrackIDs.GetArray(), fTrackIDs.GetArray() + fTrackIDs.GetSize());
}

//__________________________________________________________________________________________________
AliVParticle* AliEmcalJet::GetLeadingTrack(TClonesArray *tracks) const
{
  AliVParticle* maxTrack = 0;
  for (Int_t i = 0; i < GetNumberOfTracks(); i++) {
    AliVParticle *track = TrackAt(i, tracks);
    if (!maxTrack || track->Pt() > maxTrack->Pt()) 
      maxTrack = track;
  }

  return maxTrack;
}

//__________________________________________________________________________________________________
AliVCluster* AliEmcalJet::GetLeadingCluster(TClonesArray *clusters) const
{
  AliVCluster* maxCluster = 0;
  for (Int_t i = 0; i < GetNumberOfClusters(); i++) {
    AliVCluster *cluster = ClusterAt(i, clusters);
    if (!maxCluster || cluster->E() > maxCluster->E()) 
      maxCluster = cluster;
  }

  return maxCluster;
}

//__________________________________________________________________________________________________
void AliEmcalJet::ResetMatching()
{
  fClosestJets[0] = 0;
  fClosestJets[1] = 0; 
  fClosestJetsDist[0] = 999; 
  fClosestJetsDist[1] = 999; 
  fMatched = 2;
}
