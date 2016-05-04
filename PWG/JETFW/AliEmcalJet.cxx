//
// Emcal jet class.
//
// Author: C.Loizides

#include "AliEmcalJet.h"

#include "AliLog.h"
#include "Riostream.h"

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
  fAreaE(0),
  fAreaEmc(-1),
  fAxisInEmcal(0),
  fFlavourTagging(0),
  fFlavourTracks(0),
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
  fTaggedJet(0x0),
  fTagStatus(-1),
  fPtSub(0),
  fPtSubVect(0),
  fTriggers(0),
  fLabel(-1),
  fHasGhost(kFALSE),
  fGhosts(),
  fJetShapeProperties(0)
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
  fPt(TMath::Sqrt(px * px + py* py)),
  fEta(TMath::ASinH(pz / fPt)),
  fPhi(0),
  fM(0),
  fNEF(0),
  fArea(0),
  fAreaEta(0),
  fAreaPhi(0),
  fAreaE(0),
  fAreaEmc(-1),
  fAxisInEmcal(0),
  fFlavourTagging(0),
  fFlavourTracks(0),
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
  fTaggedJet(0x0),
  fTagStatus(-1),
  fPtSub(0),
  fPtSubVect(0),
  fTriggers(0),
  fLabel(-1),
  fHasGhost(kFALSE),
  fGhosts(),
  fJetShapeProperties(0)
{
  // Constructor.

  if(fPt != 0) {
    fPhi = TVector2::Phi_0_2pi(TMath::ATan2(py, px));
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
  fAreaE(0),
  fAreaEmc(-1),
  fAxisInEmcal(0),
  fFlavourTagging(0),
  fFlavourTracks(0),
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
  fTaggedJet(0x0),
  fTagStatus(-1),
  fPtSub(0),
  fPtSubVect(0),
  fTriggers(0),
  fLabel(-1),
  fHasGhost(kFALSE),
  fGhosts(),
  fJetShapeProperties(0)
{
  // Constructor.

  fPhi = TVector2::Phi_0_2pi(fPhi);

  fClosestJets[0] = 0;
  fClosestJets[1] = 0;
  fClosestJetsDist[0] = 999;
  fClosestJetsDist[1] = 999;
}

//_________________________________________________________________________________________________
AliEmcalJet::AliEmcalJet(const AliEmcalJet& jet) :
  AliVParticle(jet),
  fPt(jet.fPt),
  fEta(jet.fEta),
  fPhi(jet.fPhi),
  fM(jet.fM),
  fNEF(jet.fNEF),
  fArea(jet.fArea),
  fAreaEta(jet.fAreaEta),
  fAreaPhi(jet.fAreaPhi),
  fAreaE(jet.fAreaE),
  fAreaEmc(jet.fAreaEmc),
  fAxisInEmcal(jet.fAxisInEmcal),
  fFlavourTagging(jet.fFlavourTagging),
  fFlavourTracks((jet.fFlavourTracks) ? new TObjArray(*(jet.fFlavourTracks)) : 0),
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
  fTaggedJet(jet.fTaggedJet),
  fTagStatus(jet.fTagStatus),
  fPtSub(jet.fPtSub),
  fPtSubVect(jet.fPtSubVect),
  fTriggers(jet.fTriggers),
  fLabel(jet.fLabel),
  fHasGhost(jet.fHasGhost),
  fGhosts(jet.fGhosts),
  fJetShapeProperties(0)
{
  // Copy constructor.
  fClosestJets[0]     = jet.fClosestJets[0];
  fClosestJets[1]     = jet.fClosestJets[1];
  fClosestJetsDist[0] = jet.fClosestJetsDist[0];
  fClosestJetsDist[1] = jet.fClosestJetsDist[1];

  if (jet.fJetShapeProperties) {
    fJetShapeProperties = new AliEmcalJetShapeProperties(*(jet.fJetShapeProperties));
  }
}

//_________________________________________________________________________________________________
AliEmcalJet::~AliEmcalJet()
{
  if (fJetShapeProperties) delete fJetShapeProperties;
}

//_________________________________________________________________________________________________
AliEmcalJet& AliEmcalJet::operator=(const AliEmcalJet& jet)
{
  // Assignment operator.

  if(this != &jet) {
    AliVParticle::operator=(jet);
    fPt                 = jet.fPt;
    fEta                = jet.fEta;
    fPhi                = jet.fPhi;
    fM                  = jet.fM;
    fNEF                = jet.fNEF;
    fArea               = jet.fArea;
    fAreaEta            = jet.fAreaEta;
    fAreaPhi            = jet.fAreaPhi;
    fAreaE              = jet.fAreaE;
    fAreaEmc            = jet.fAreaEmc;
    fAxisInEmcal        = jet.fAxisInEmcal;
    fFlavourTagging     = jet.fFlavourTagging;
    fFlavourTracks      = (jet.fFlavourTracks) ? new TObjArray(*(jet.fFlavourTracks)) : 0;
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
    fTaggedJet          = jet.fTaggedJet;
    fTagStatus          = jet.fTagStatus;
    fPtSub              = jet.fPtSub;
    fPtSubVect          = jet.fPtSubVect;
    fTriggers           = jet.fTriggers;
    fLabel              = jet.fLabel;
    fHasGhost = jet.fHasGhost;
    fGhosts   = jet.fGhosts;
    if (jet.fJetShapeProperties) {
      fJetShapeProperties = new AliEmcalJetShapeProperties(*(jet.fJetShapeProperties));
    }
  }

  return *this;
}

//_________________________________________________________________________________________________
Int_t AliEmcalJet::Compare(const TObject* obj) const
{
  //Return -1 if this is smaller than obj, 0 if objects are equal and 1 if this is larger than obj.

  const AliEmcalJet* jet = static_cast<const AliEmcalJet*>(obj);
  if(!obj)
    return 0;
  if(Pt() > jet->Pt())
    return -1;
  return 1;
}

//__________________________________________________________________________________________________
void AliEmcalJet::GetMomentum(TLorentzVector& vec) const
{
  // Return momentum as four-vector.

  vec.SetPtEtaPhiE(fPt, fEta, fPhi, E());
}

//__________________________________________________________________________________________________
Double_t AliEmcalJet::PtSub(Double_t rho, Bool_t save)
{
  // Return transverse momentum after scalar subtraction. Save the result if required.
  // Result can be negative.

  Double_t ptcorr = fPt - rho * fArea;
  if(save)
    fPtSub = ptcorr;
  return ptcorr;
}

//__________________________________________________________________________________________________
Double_t AliEmcalJet::PtSubVect(Double_t rho, Bool_t save)
{
  // Return transverse momentum after vectorial subtraction. Save the result if required.
  // Result cannot be negative.

  Double_t dx = Px() - rho * fArea * TMath::Cos(fAreaPhi);
  Double_t dy = Py() - rho * fArea * TMath::Sin(fAreaPhi);
  //Double_t dz = Pz() - rho * fArea * TMath::SinH(fAreaEta);
  Double_t ptcorr = TMath::Sqrt(dx * dx + dy * dy);
  if(save)
    fPtSubVect = ptcorr;
  return ptcorr;
}

//__________________________________________________________________________________________________
TLorentzVector AliEmcalJet::SubtractRhoVect(Double_t rho, Bool_t save)
{
  // Return four-momentum after vectorial subtraction. Save pt if required.
  // Saved value of pt is negative if the corrected momentum is pointing to the opposite half-plane in the x-y plane w.r.t. the raw momentum.

  TLorentzVector vecCorr;
  GetMom(vecCorr);
  TLorentzVector vecBg;
  vecBg.SetPtEtaPhiE(fArea, fAreaEta, fAreaPhi, fAreaE);
  vecBg *= rho;
  vecCorr -= vecBg;
  if(save)
  {
    Double_t dPhi = TMath::Abs(TVector2::Phi_mpi_pi(Phi() - vecCorr.Phi()));
    Int_t signum = dPhi <= TMath::PiOver2() ? 1 : -1;
    fPtSubVect = signum * vecCorr.Pt();
  }
  return vecCorr;
}

//__________________________________________________________________________________________________
void AliEmcalJet::SortConstituents()
{
  // Sort constituent by index (increasing).

  std::sort(fClusterIDs.GetArray(), fClusterIDs.GetArray() + fClusterIDs.GetSize());
  std::sort(fTrackIDs.GetArray(), fTrackIDs.GetArray() + fTrackIDs.GetSize());
}

//__________________________________________________________________________________________________
Double_t AliEmcalJet::DeltaR(const AliVParticle* part) const
{
  // Helper function to calculate the distance between two jets or a jet and a particle

  Double_t dPhi = Phi() - part->Phi();
  Double_t dEta = Eta() - part->Eta();
  dPhi = TVector2::Phi_mpi_pi(dPhi);
  return TMath::Sqrt(dPhi * dPhi + dEta * dEta);
}


//__________________________________________________________________________________________________
std::vector<int> AliEmcalJet::SortConstituentsPt(TClonesArray* tracks) const
{
  // Sorting by p_T (decreasing) jet constituents

  typedef std::pair<Double_t, Int_t> ptidx_pair;

  // Create vector for Pt sorting
  std::vector<ptidx_pair> pair_list ;

  for(Int_t i_entry = 0; i_entry < GetNumberOfTracks(); i_entry++)
  {
    AliVParticle* track = TrackAt(i_entry, tracks);
    if(!track)
    {
      AliError(Form("Unable to find jet track %d in collection %s (pos in collection %d, max %d)", i_entry, tracks->GetName(), TrackAt(i_entry), tracks->GetEntriesFast()));
      continue;
    }

    pair_list.push_back(std::make_pair(track->Pt(), i_entry));
  }

  std::stable_sort(pair_list.begin() , pair_list.end() , sort_descend());

  // return an vector of indexes of constituents (sorted descending by pt)
  std::vector <int> index_sorted_list;

  for(std::vector< std::pair<Double_t, Int_t> >::iterator it = pair_list.begin(); it != pair_list.end(); ++it)
  { index_sorted_list.push_back((*it).second); }   // populating the return object with indexes of sorted tracks

  return index_sorted_list;
}

//________________________________________________________________________
Double_t AliEmcalJet::GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz) const
{
  // Get the z of a constituent inside of a jet

  Double_t pJetSq = P();
  pJetSq *= pJetSq;

  if(pJetSq > 1e-6)
  { return (trkPx * Px() + trkPy * Py() + trkPz * Pz()) / pJetSq ; }
  else
  { AliWarning(Form("%s: strange, pjet*pjet seems to be zero pJetSq: %f", GetName(), pJetSq)); return -1; }
}

//________________________________________________________________________
Double_t AliEmcalJet::GetZ(const AliVParticle* trk) const
{
  // Get Z of constituent trk

  return GetZ(trk->Px(), trk->Py(), trk->Pz());
}

//__________________________________________________________________________________________________
AliVParticle* AliEmcalJet::GetLeadingTrack(TClonesArray* tracks) const
{
  AliVParticle* maxTrack = 0;
  for(Int_t i = 0; i < GetNumberOfTracks(); i++) {
    AliVParticle* track = TrackAt(i, tracks);
    if(!track) {
      AliError(Form("Unable to find jet track %d in collection %s (pos in collection %d, max %d)",
                    i, tracks->GetName(), TrackAt(i), tracks->GetEntriesFast()));
      continue;
    }
    if(!maxTrack || track->Pt() > maxTrack->Pt())
      maxTrack = track;
  }

  return maxTrack;
}

//__________________________________________________________________________________________________
AliVCluster* AliEmcalJet::GetLeadingCluster(TClonesArray* clusters) const
{
  AliVCluster* maxCluster = 0;
  for(Int_t i = 0; i < GetNumberOfClusters(); i++) {
    AliVCluster* cluster = ClusterAt(i, clusters);
    if(!cluster) {
      AliError(Form("Unable to find jet cluster %d in collection %s (pos in collection %d, max %d)",
                    i, clusters->GetName(), ClusterAt(i), clusters->GetEntriesFast()));
      continue;
    }
    if(!maxCluster || cluster->E() > maxCluster->E())
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

//__________________________________________________________________________________________________
Int_t AliEmcalJet::ContainsTrack(Int_t it) const
{
  for (Int_t i = 0; i < fTrackIDs.GetSize(); i++) {
    if (it == fTrackIDs[i]) return i;
  }
  return -1;
}

//__________________________________________________________________________________________________
Int_t AliEmcalJet::ContainsCluster(Int_t ic) const
{
  for (Int_t i = 0; i < fClusterIDs.GetSize(); i++) {
    if (ic == fClusterIDs[i]) return i;
  }
  return -1;
}

/**
 * Create string representation of the Jet. The string respresentation contains
 * - \f$ p_{t} \f$ of the jet
 * - \f$ \eta \f$ of the jet
 * - \f$ \phi \f$ of the jet
 * - \f$ p_{t} \f$ of the maximum charged and neutral particle
 * - Number of constituent tracks
 * - Number of constituent clusters
 * - Jet area
 * - Neutral energy fraction
 * @return String representation of the jet
 */
TString AliEmcalJet::toString() const {
  return TString::Format("Jet pT = %.2f, eta = %.2f, phi = %.2f, max charged pT = %.2f, max neutral pT = %.2f, N tracks = %d, N clusters = %d, Area = %.2f, NEF = %.2f",
         Pt(), Eta(), Phi(), MaxChargedPt(), MaxNeutralPt(), GetNumberOfTracks(), GetNumberOfClusters(), Area(), NEF());
}

/**
 * Print basic jet information using the string representation provided by
 * AliEmcalJet::toString
 *
 * @param unused
 */
void AliEmcalJet::Print(Option_t* /*opt*/) const
{
  Printf("%s\n", toString().Data());
}

/**
 * Print basic jet information on an output stream using the string representation provided by
 * AliEmcalJet::toString. Used by operator<<
 * @param in output stream stream
 * @return reference to the output stream
 */
std::ostream &AliEmcalJet::Print(std::ostream &in) const {
  in << toString().Data();
  return in;
}

//________________________________________________________________________
void AliEmcalJet::PrintConstituents(TClonesArray* tracks, TClonesArray* clusters) const
{
  if (tracks) {
    for (Int_t i = 0; i < GetNumberOfTracks(); i++) {
      AliVParticle* part = TrackAt(i, tracks);
      if (part) {
        Printf("Track %d (index = %d) pT = %.2f, eta = %.2f, phi = %.2f, PDG code = %d", i, TrackAt(i), part->Pt(), part->Eta(), part->Phi(), part->PdgCode());
      }
    }
  }

  if (clusters) {
    for (Int_t i = 0; i < GetNumberOfClusters(); i++) {
      AliVCluster* clus = ClusterAt(i, clusters);
      if (clus) {
        Printf("Cluster %d (index = %d) E = %.2f", i, ClusterAt(i), clus->E());
      }
    }
  }
}

//________________________________________________________________________
Double_t AliEmcalJet::GetFlavourTrackZ(Int_t i)  const
{
  if (P() < 1e-6) return 0.;
  AliVParticle* hftrack = GetFlavourTrack(i);
  return hftrack != 0 ? hftrack->P() / P() : 0.;
}

/**
 * Implementation of the output stream operator for AliEmcalJet. Printing
 * basic jet information provided by function toString
 * @param in output stream
 * @param myjet Jet which will be printed
 * @return Reference to the output stream
 */
std::ostream &operator<<(std::ostream &in, const AliEmcalJet &myjet) {
  std::ostream &result = myjet.Print(in);
  return result;
}
