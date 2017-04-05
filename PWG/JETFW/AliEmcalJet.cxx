/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliEmcalJet.h"

#include "AliLog.h"
#include "Riostream.h"

#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

/// \cond CLASSIMP
ClassImp(AliEmcalJet);
/// \endcond

/**
 * Default constructor
 */
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
  fJetShapeProperties(0),
  fJetAcceptanceType(0)
{
  fClosestJets[0] = 0;
  fClosestJets[1] = 0;
  fClosestJetsDist[0] = 999;
  fClosestJetsDist[1] = 999;
}

/**
 * Constructor that uses the 3-momentum to define the jet axis.
 * It assumes zero mass for the jet.
 * @param px First transverse component of the jet momentum
 * @param px Second transverse component of the jet momentum
 * @param pz Longitudinal component of the jet momentum
 */
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
  fJetShapeProperties(0),
  fJetAcceptanceType(0)
{
  if (fPt != 0) {
    fPhi = TVector2::Phi_0_2pi(TMath::ATan2(py, px));
  }

  fClosestJets[0] = 0;
  fClosestJets[1] = 0;
  fClosestJetsDist[0] = 999;
  fClosestJetsDist[1] = 999;
}

/**
 * Constructor that uses the 4-momentum to define the jet axis.
 * Coordinates are given in cylindrical system plus the mass.
 * @param pt Transverse component of the jet momentum
 * @param eta Pseudo-rapidity of the jet
 * @param phi Azimuthal angle of the jet axis
 * @param m Mass of the jet
 */
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
  fJetShapeProperties(0),
  fJetAcceptanceType(0)
{
  fPhi = TVector2::Phi_0_2pi(fPhi);

  fClosestJets[0] = 0;
  fClosestJets[1] = 0;
  fClosestJetsDist[0] = 999;
  fClosestJetsDist[1] = 999;
}

/**
 * Copy constructor.
 * @param jet Constant reference to copy from
 */
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
  fJetShapeProperties(0),
  fJetAcceptanceType(jet.fJetAcceptanceType)
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

/**
 * Destructor.
 */
AliEmcalJet::~AliEmcalJet()
{
  if (fJetShapeProperties) delete fJetShapeProperties;
}

/**
 * Assignment operator
 * @param jet Constant reference to copy from
 * @return A reference to this
 */
AliEmcalJet& AliEmcalJet::operator=(const AliEmcalJet& jet)
{
  if (this != &jet) {
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
    fJetAcceptanceType  = jet.fJetAcceptanceType;
  }

  return *this;
}

/**
 * Compares two instances of AliEmcalJet, ordering them based on their transverse momentum.
 * @param obj Pointer to another instance of AliEmcalJet
 * @return -1 if this is smaller than obj, 1 if this is larger than obj, 0 if objects are equal or if obj is NULL or not an instance of AliEmcalJet
 */
Int_t AliEmcalJet::Compare(const TObject* obj) const
{
  //Return -1 if this is smaller than obj, 0 if objects are equal and 1 if this is larger than obj.

  if (obj == this) return 0;

  const AliEmcalJet* jet = dynamic_cast<const AliEmcalJet*>(obj);
  if (!jet) return 0;

  if (Pt() > jet->Pt()) return -1;
  else if (Pt() < jet->Pt()) return 1;
  else return 0;
}

/**
 * Builds a 4-momentum object using information contained in this instance of AliEmcalJet
 * @param[out] vec Reference to TLorentzVector where the 4-momentum is returned
 */
void AliEmcalJet::GetMomentum(TLorentzVector& vec) const
{
  vec.SetPtEtaPhiE(fPt, fEta, fPhi, E());
}

/**
 * Calculates transverse momentum after scalar average background subtraction.
 * The result can be negative. It saves the result if requested.
 * @param rho Scalar average background
 * @param save If kTRUE, stores the result in a class field
 * @return The subtracted transverse momentum
 */
Double_t AliEmcalJet::PtSub(Double_t rho, Bool_t save)
{
  Double_t ptcorr = fPt - rho * fArea;
  if (save) fPtSub = ptcorr;
  return ptcorr;
}

/**
 * Calculates transverse momentum after vectorial average background subtraction.
 * The result can be negative. It saves the result if requested.
 * @param rho Vectorial average background
 * @param save If kTRUE, stores the result in a class field
 * @return The subtracted transverse momentum
 */
Double_t AliEmcalJet::PtSubVect(Double_t rho, Bool_t save)
{
  Double_t dx = Px() - rho * fArea * TMath::Cos(fAreaPhi);
  Double_t dy = Py() - rho * fArea * TMath::Sin(fAreaPhi);
  //Double_t dz = Pz() - rho * fArea * TMath::SinH(fAreaEta);
  Double_t ptcorr = TMath::Sqrt(dx * dx + dy * dy);
  if (save) fPtSubVect = ptcorr;
  return ptcorr;
}

/**
 * Calculates 4-momentum after vectorial average background subtraction.
 * It saves the result if requested.
 * @param rho Vectorial average background
 * @param save If kTRUE, stores the result in a class field
 * @return The subtracted 4-momentum
 */
TLorentzVector AliEmcalJet::SubtractRhoVect(Double_t rho, Bool_t save)
{
  TLorentzVector vecCorr;
  GetMomentum(vecCorr);
  TLorentzVector vecBg;
  vecBg.SetPtEtaPhiE(fArea, fAreaEta, fAreaPhi, fAreaE);
  vecBg *= rho;
  vecCorr -= vecBg;
  if (save) {
    Double_t dPhi = TMath::Abs(TVector2::Phi_mpi_pi(Phi() - vecCorr.Phi()));
    Int_t signum = dPhi <= TMath::PiOver2() ? 1 : -1;
    fPtSubVect = signum * vecCorr.Pt();
  }
  return vecCorr;
}

/**
 *  Sort constituent by index (increasing).
 *
 */
void AliEmcalJet::SortConstituents()
{
  std::sort(fClusterIDs.GetArray(), fClusterIDs.GetArray() + fClusterIDs.GetSize());
  std::sort(fTrackIDs.GetArray(), fTrackIDs.GetArray() + fTrackIDs.GetSize());
}

/**
 * Helper function to calculate the distance between two jets or a jet and a particle
 * @param part Constant pointer to another particle
 * @return Distance in the eta-phi phase space
 */
Double_t AliEmcalJet::DeltaR(const AliVParticle* part) const
{
  Double_t dPhi = Phi() - part->Phi();
  Double_t dEta = Eta() - part->Eta();
  dPhi = TVector2::Phi_mpi_pi(dPhi);
  return TMath::Sqrt(dPhi * dPhi + dEta * dEta);
}

/**
 * Sorting jet constituents by pT (decreasing)  
 * It returns a standard vector with the indexes of the constituents relative to fTrackIDs.
 * To retrieve the track do:
 * ~~~{.cxx}
 * TClonesArray* fTracksContArray = jetCont->GetParticleContainer()->GetArray();
 * std::vector< int > index_sorted_list = jet->GetPtSortedTrackConstituentIndexes(fTracksContArray);
 * for (std::size_t i = 0; i < jet->GetNumberOfTracks(); i++ ) {
 * track = jet->TrackAt ( index_sorted_list.at (i), fTracksContArray );
 * // use track;
 * }
 * ~~~
 * @param tracks Array containing pointers to the tracks from which jet constituents are drawn
 * @return Standard vector with the list of constituent indexes (relative to fTrackIDs)
 */
std::vector<int> AliEmcalJet::GetPtSortedTrackConstituentIndexes(TClonesArray* tracks) const
{
  typedef std::pair<Double_t, Int_t> ptidx_pair;

  // Create vector for Pt sorting
  std::vector<ptidx_pair> pair_list;

  for (Int_t i_entry = 0; i_entry < GetNumberOfTracks(); i_entry++) {
    AliVParticle* track = TrackAt(i_entry, tracks);
    if (!track) {
      AliError(Form("Unable to find jet track %d in collection %s (pos in collection %d, max %d)", i_entry, tracks->GetName(), TrackAt(i_entry), tracks->GetEntriesFast()));
      continue;
    }
    pair_list.push_back(std::make_pair(track->Pt(), i_entry));
  }

  std::stable_sort(pair_list.begin() , pair_list.end() , sort_descend());

  // return a vector of indexes of constituents (sorted descending by pt)
  std::vector <int> index_sorted_list;

  // populating the return object with indexes of sorted tracks
  for (auto it : pair_list) index_sorted_list.push_back(it.second);

  return index_sorted_list;
}

/**
 * Get the momentum fraction of a jet constituent
 * @param trkPx First transverse component of the momentum of the jet constituent
 * @param trkPy Second transverse component of the momentum of the jet constituent
 * @param trkPz Longitudinal component of the momentum of the jet constituent
 * @return Momentum fraction
 */
Double_t AliEmcalJet::GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz) const
{
  Double_t pJetSq = P();
  pJetSq *= pJetSq;

  if (pJetSq > 1e-6) {
    return (trkPx * Px() + trkPy * Py() + trkPz * Pz()) / pJetSq;
  }
  else {
    AliWarning(Form("%s: strange, pjet*pjet seems to be zero pJetSq: %f", GetName(), pJetSq));
    return -1;
  }
}

/**
 * Get the momentum fraction of a jet constituent
 * @param trk Jet constituent
 * @return Momentum fraction
 */
Double_t AliEmcalJet::GetZ(const AliVParticle* trk) const
{
  return GetZ(trk->Px(), trk->Py(), trk->Pz());
}

/**
 * Find the leading track constituent of the jet.
 * @param tracks Array containing the pointers to the tracks from which jet constituents are drawn
 * @return Pointer to the leading track of the jet
 */
AliVParticle* AliEmcalJet::GetLeadingTrack(TClonesArray* tracks) const
{
  AliVParticle* maxTrack = 0;
  for (Int_t i = 0; i < GetNumberOfTracks(); i++) {
    AliVParticle* track = TrackAt(i, tracks);
    if (!track) {
      AliError(Form("Unable to find jet track %d in collection %s (pos in collection %d, max %d)",
          i, tracks->GetName(), TrackAt(i), tracks->GetEntriesFast()));
      continue;
    }
    if (!maxTrack || track->Pt() > maxTrack->Pt())
      maxTrack = track;
  }

  return maxTrack;
}

/**
 * Find the leading cluster constituent of the jet.
 * @param tracks Array containing the pointers to the clusters from which jet constituents are drawn
 * @return Pointer to the leading cluster of the jet
 */
AliVCluster* AliEmcalJet::GetLeadingCluster(TClonesArray* clusters) const
{
  AliVCluster* maxCluster = 0;
  for (Int_t i = 0; i < GetNumberOfClusters(); i++) {
    AliVCluster* cluster = ClusterAt(i, clusters);
    if (!cluster) {
      AliError(Form("Unable to find jet cluster %d in collection %s (pos in collection %d, max %d)",
          i, clusters->GetName(), ClusterAt(i), clusters->GetEntriesFast()));
      continue;
    }
    if (!maxCluster || cluster->E() > maxCluster->E())
      maxCluster = cluster;
  }

  return maxCluster;
}

/**
 * Reset jet matching information.
 */
void AliEmcalJet::ResetMatching()
{
  fClosestJets[0] = 0;
  fClosestJets[1] = 0;
  fClosestJetsDist[0] = 999;
  fClosestJetsDist[1] = 999;
  fMatched = 2;
}

/**
 * Checks whether a certain track is among the jet constituents by looking for its index
 * @param Index of the track to search
 * @return The position of the track in the jet constituent array, if the track is found; -1 if the track is not a jet constituent
 */
Int_t AliEmcalJet::ContainsTrack(Int_t it) const
{
  for (Int_t i = 0; i < fTrackIDs.GetSize(); i++) {
    if (it == fTrackIDs[i]) return i;
  }
  return -1;
}

/**
 * Checks whether a certain cluster is among the jet constituents by looking for its index
 * @param Index of the cluster to search
 * @return The position of the cluster in the jet constituent array, if the cluster is found; -1 if the cluster is not a jet constituent
 */
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

/**
 * Prints the list of constituents in the standard output
 * @param tracks Array containing the pointers to tracks
 * @param clusters Array containing the pointers to the clusters
 */
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

/**
 * Calculates the momentum fraction carried by the "flavor" track
 * @param i Position of the flavor track in the flavor track constituent array
 * @return The momentum fraction carried by the flavor track
 */
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

/**
 * Add a ghost particle to the ghost particle array.
 * This function should be called by the jet finder to fill the ghost particle array.
 * @param dPx First component of the transverse momentum of the particle
 * @param dPy First component of the transverse momentum of the particle
 * @param dPz Longitudinal component of the momentum of the particle
 * @param dE Energy of the particle
 */
void AliEmcalJet::AddGhost(const Double_t dPx, const Double_t dPy, const Double_t dPz, const Double_t dE)
{
  TLorentzVector ghost(dPx, dPy, dPz, dE);
  fGhosts.push_back(ghost);
  if (!fHasGhost) fHasGhost = kTRUE;
  return;
}

/**
 * Clear this object: remove matching information, jet constituents, ghosts
 */
void AliEmcalJet::Clear(Option_t */*option*/)
{
  fClusterIDs.Set(0);
  fTrackIDs.Set(0);
  fClosestJets[0] = 0;
  fClosestJets[1] = 0;
  fClosestJetsDist[0] = 0;
  fClosestJetsDist[1] = 0;
  fMatched = 0;
  fPtSub = 0;
  fGhosts.clear();
  fHasGhost = kFALSE;
}

/**
 * Retrieve the track constituent corresponding to the index found at a certain position.
 * Automatically retrieves the particle from the proper TClonesArray. This function is preferred to
 * TrackAt(Int_t, TClonesArray), which is only kept for backwards compatibility.
 *
 * @param idx Position of the track constituent
 * @return Pointer to the track constituent requested (if found)
 */
AliVParticle* AliEmcalJet::Track(Int_t idx) const
{
  return AliParticleContainer::GetEmcalContainerIndexMap().GetObjectFromGlobalIndex(TrackAt(idx));
}

/**
 * Finds the track constituent corresponding to the index found at a certain position.
 * @param idx Position of the track constituent
 * @param ta Array with pointers to the tracks from which jet constituents are drawn
 * @return Pointer to the track constituent requested (if found)
 */
AliVParticle* AliEmcalJet::TrackAt(Int_t idx, TClonesArray *ta) const
{
  if (!ta) return 0;
  auto res =  AliParticleContainer::GetEmcalContainerIndexMap().LocalIndexFromGlobalIndex(TrackAt(idx));
  if (res.second != ta) {
    AliWarning(Form("TClonesArray %s that was passed does not correspond to the passed index! The index belongs to a different TClonesArray named %s! Returning the object corresponding to the index (not the passed TClonesArray)! Consider fixing by updating to jet->Track(index).", ta->GetName(), res.second->GetName()));
  }
  return dynamic_cast<AliVParticle*>(res.second->At(res.first));
}

/**
 * Checks whether a given track is among the jet constituents
 * @param track Pointer to the track to be searched
 * @param tracks Array with pointers to the tracks from which jet constituents are drawn
 * @return Position of the track among the jet constituents, if the track is found; -1 otherwise
 */
Int_t AliEmcalJet::ContainsTrack(AliVParticle* track, TClonesArray* tracks) const
{
  if (!tracks || !track) return 0;
  return ContainsTrack(tracks->IndexOf(track));
}

/**
 * Retrieve the cluster constituent corresponding to the index found at a certain position.
 * Automatically retrieves the particle from the proper TClonesArray. This function is preferred to
 * ClusterAt(Int_t, TClonesArray), which is only kept for backwards compatibility.
 *
 * @param idx Position of the cluster constituent
 * @return Pointer to the cluster constituent requested (if found)
 */
AliVCluster* AliEmcalJet::Cluster(Int_t idx) const
{
  return AliClusterContainer::GetEmcalContainerIndexMap().GetObjectFromGlobalIndex(ClusterAt(idx));
}

/**
 * Finds the cluster constituent corresponding to the index found at a certain position.
 * @param idx Position of the cluster constituent
 * @param ta Array with pointers to the clusters from which jet constituents are drawn
 * @return Pointer to the cluster constituent requested (if found)
 */
AliVCluster* AliEmcalJet::ClusterAt(Int_t idx, TClonesArray *ca) const
{
  if (!ca) return 0;

  auto res =  AliClusterContainer::GetEmcalContainerIndexMap().LocalIndexFromGlobalIndex(ClusterAt(idx));
  if (res.second != ca) {
    AliWarning(Form("TClonesArray %s that was passed does not correspond to the passed index! The index belongs to a different TClonesArray named %s! Returning the object corresponding to the index (not the passed TClonesArray)! Consider fixing by updating to jet->Cluster(index).", ca->GetName(), res.second->GetName()));
  }
  return dynamic_cast<AliVCluster*>(res.second->At(res.first));
}

/**
 * Checks whether a given cluster is among the jet constituents
 * @param cluster Pointer to the track to be searched
 * @param clusters Array with pointers to the clusters from which jet constituents are drawn
 * @return Position of the cluster among the jet constituents, if the cluster is found; -1 otherwise
 */
Int_t AliEmcalJet::ContainsCluster(AliVCluster* cluster, TClonesArray* clusters) const
{
  if (!clusters || !cluster) return 0;
  return ContainsCluster(clusters->IndexOf(cluster));
}

/**
 * Get Xi = Log(1 / z) of constituent track
 * @param trk Pointer to a constituent track
 * @return Xi of the constituent
 */
Double_t AliEmcalJet::GetXi(const AliVParticle* trk) const
{
  return TMath::Log(1 / GetZ(trk));
}

/**
 * Get Xi = Log(1 / z) of constituent track
 * @param trkPx First transverse component of the momentum of the jet constituent
 * @param trkPy Second transverse component of the momentum of the jet constituent
 * @param trkPz Longitudinal component of the momentum of the jet constituent
 * @return Xi of the constituent
 */
Double_t AliEmcalJet::GetXi( const Double_t trkPx, const Double_t trkPy, const Double_t trkPz ) const
{
  return TMath::Log(1 / GetZ(trkPx, trkPy, trkPz));
}

/**
 * Add a track to the list of flavor tagging tracks
 * @param Pointer to the flavor track
 */
void AliEmcalJet::AddFlavourTrack(AliVParticle* hftrack)
{
  if (!fFlavourTracks) fFlavourTracks = new TObjArray();
  fFlavourTracks->Add(hftrack);
}

/**
 * Finds the flavor track at a given array position
 * @param i Position of the flavor track in the array
 * @return Pointer to the flavor track
 */
AliVParticle* AliEmcalJet::GetFlavourTrack(Int_t i) const
{
  if (!fFlavourTracks || i < 0 || i >= fFlavourTracks->GetEntriesFast()) return 0;

  return  static_cast<AliVParticle*>(fFlavourTracks->At(i));
}

/**
 * Finds the flavor track at a given array position and removes it from the array
 * @param i Position of the flavor track in the array
 * @return Pointer to the flavor track that was removed
 */
AliVParticle* AliEmcalJet::RemoveFlavourTrack(Int_t i)
{
  if (!fFlavourTracks || i < 0 || i >= fFlavourTracks->GetEntriesFast()) return 0;

  return static_cast<AliVParticle*>(fFlavourTracks->RemoveAt(i));
}
