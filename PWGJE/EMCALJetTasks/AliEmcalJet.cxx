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
  fAreaEmc(-1), 
  fAxisInEmcal(0), 
  fFlavourTagging(0),
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
  fPtVectSub(0),
  fTriggers(0),
  fJetShapeMassFirstDer(0),
  fJetShapeMassSecondDer(0),
  fJetShapeMassFirstSub(0),
  fJetShapeMassSecondSub(0),
  fLabel(-1),
  fGRNumerator(0),
  fGRDenominator(0),
  fGRNumeratorSub(0),
  fGRDenominatorSub(0),
  fJetShapeAngularityFirstDer(0),
  fJetShapeAngularitySecondDer(0),
  fJetShapeAngularityFirstSub(0),
  fJetShapeAngularitySecondSub(0),
  fJetShapepTDFirstDer(0),
  fJetShapepTDSecondDer(0),
  fJetShapepTDFirstSub(0),
  fJetShapepTDSecondSub(0),
  fJetShapeCircularityFirstDer(0),
  fJetShapeCircularitySecondDer(0),
  fJetShapeCircularityFirstSub(0),
  fJetShapeCircularitySecondSub(0),
  fJetShapeConstituentFirstDer(0),
  fJetShapeConstituentSecondDer(0),
  fJetShapeConstituentFirstSub(0),
  fJetShapeConstituentSecondSub(0),
  fJetShapeLeSubFirstDer(0),
  fJetShapeLeSubSecondDer(0),
  fJetShapeLeSubFirstSub(0),
  fJetShapeLeSubSecondSub(0)
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
  fFlavourTagging(0),
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
  fPtVectSub(0),
  fTriggers(0),
  fJetShapeMassFirstDer(0),
  fJetShapeMassSecondDer(0),
  fJetShapeMassFirstSub(0),
  fJetShapeMassSecondSub(0),
  fLabel(-1),
  fGRNumerator(0),
  fGRDenominator(0),
  fGRNumeratorSub(0),
  fGRDenominatorSub(0),
    fJetShapeAngularityFirstDer(0),
  fJetShapeAngularitySecondDer(0),
  fJetShapeAngularityFirstSub(0),
  fJetShapeAngularitySecondSub(0),
   fJetShapepTDFirstDer(0),
  fJetShapepTDSecondDer(0),
  fJetShapepTDFirstSub(0),
  fJetShapepTDSecondSub(0),
  fJetShapeCircularityFirstDer(0),
  fJetShapeCircularitySecondDer(0),
  fJetShapeCircularityFirstSub(0),
  fJetShapeCircularitySecondSub(0),
  fJetShapeConstituentFirstDer(0),
  fJetShapeConstituentSecondDer(0),
  fJetShapeConstituentFirstSub(0),
  fJetShapeConstituentSecondSub(0),
   fJetShapeLeSubFirstDer(0),
  fJetShapeLeSubSecondDer(0),
  fJetShapeLeSubFirstSub(0),
  fJetShapeLeSubSecondSub(0)
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
  fFlavourTagging(0),
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
  fPtVectSub(0),
  fTriggers(0),
  fJetShapeMassFirstDer(0),
  fJetShapeMassSecondDer(0),
  fJetShapeMassFirstSub(0),
  fJetShapeMassSecondSub(0),
  fLabel(-1),
  fGRNumerator(0),
  fGRDenominator(0),
  fGRNumeratorSub(0),
  fGRDenominatorSub(0),
  fJetShapeAngularityFirstDer(0),
  fJetShapeAngularitySecondDer(0),
  fJetShapeAngularityFirstSub(0),
  fJetShapeAngularitySecondSub(0),
  fJetShapepTDFirstDer(0),
  fJetShapepTDSecondDer(0),
  fJetShapepTDFirstSub(0),
  fJetShapepTDSecondSub(0),
  fJetShapeCircularityFirstDer(0),
  fJetShapeCircularitySecondDer(0),
  fJetShapeCircularityFirstSub(0),
  fJetShapeCircularitySecondSub(0),
  fJetShapeConstituentFirstDer(0),
  fJetShapeConstituentSecondDer(0),
  fJetShapeConstituentFirstSub(0),
  fJetShapeConstituentSecondSub(0),
  fJetShapeLeSubFirstDer(0),
  fJetShapeLeSubSecondDer(0),
  fJetShapeLeSubFirstSub(0),
  fJetShapeLeSubSecondSub(0)

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
  fFlavourTagging(jet.fFlavourTagging),
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
  fPtVectSub(jet.fPtVectSub),
  fTriggers(jet.fTriggers),
  fJetShapeMassFirstDer(jet.fJetShapeMassFirstDer),
  fJetShapeMassSecondDer(jet.fJetShapeMassSecondDer),
  fJetShapeMassFirstSub(jet.fJetShapeMassFirstSub),
  fJetShapeMassSecondSub(jet.fJetShapeMassSecondSub),
  fLabel(jet.fLabel),
  fGRNumerator(jet.fGRNumerator),
  fGRDenominator(jet.fGRDenominator),
  fGRNumeratorSub(jet.fGRNumeratorSub),
  fGRDenominatorSub(jet.fGRDenominatorSub),
   fJetShapeAngularityFirstDer(jet.fJetShapeAngularityFirstDer),
  fJetShapeAngularitySecondDer(jet.fJetShapeAngularitySecondDer),
  fJetShapeAngularityFirstSub(jet.fJetShapeAngularityFirstSub),
  fJetShapeAngularitySecondSub(jet.fJetShapeAngularitySecondSub),

  fJetShapepTDFirstDer(jet.fJetShapepTDFirstDer),
  fJetShapepTDSecondDer(jet.fJetShapepTDSecondDer),
  fJetShapepTDFirstSub(jet.fJetShapepTDFirstSub),
  fJetShapepTDSecondSub(jet.fJetShapepTDSecondSub),

  fJetShapeCircularityFirstDer(jet.fJetShapeCircularityFirstDer),
  fJetShapeCircularitySecondDer(jet.fJetShapeCircularitySecondDer),
  fJetShapeCircularityFirstSub(jet.fJetShapeCircularityFirstSub),
  fJetShapeCircularitySecondSub(jet.fJetShapeCircularitySecondSub),

  fJetShapeConstituentFirstDer(jet.fJetShapeConstituentFirstDer),
  fJetShapeConstituentSecondDer(jet.fJetShapeConstituentSecondDer),
  fJetShapeConstituentFirstSub(jet.fJetShapeConstituentFirstSub),
  fJetShapeConstituentSecondSub(jet.fJetShapeConstituentSecondSub),

  fJetShapeLeSubFirstDer(jet.fJetShapeLeSubFirstDer),
  fJetShapeLeSubSecondDer(jet.fJetShapeLeSubSecondDer),
  fJetShapeLeSubFirstSub(jet.fJetShapeLeSubFirstSub),
  fJetShapeLeSubSecondSub(jet.fJetShapeLeSubSecondSub)
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
    fFlavourTagging     = jet.fFlavourTagging;
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
    fPtVectSub          = jet.fPtVectSub;
    fTriggers           = jet.fTriggers;
    fJetShapeMassFirstDer  = jet.fJetShapeMassFirstDer;
    fJetShapeMassSecondDer = jet.fJetShapeMassSecondDer;
    fJetShapeMassFirstSub  = jet.fJetShapeMassFirstSub;
    fJetShapeMassSecondSub = jet.fJetShapeMassSecondSub;
    fLabel              = jet.fLabel;
    fGRNumerator        = jet.fGRNumerator;
    fGRDenominator      = jet.fGRDenominator;
    fGRNumeratorSub     = jet.fGRNumeratorSub;
    fGRDenominatorSub   = jet.fGRDenominatorSub;
           fJetShapeAngularityFirstDer  = jet.fJetShapeAngularityFirstDer;
    fJetShapeAngularitySecondDer = jet.fJetShapeAngularitySecondDer;
    fJetShapeAngularityFirstSub  = jet.fJetShapeAngularityFirstSub;
    fJetShapeAngularitySecondSub = jet.fJetShapeAngularitySecondSub;

    fJetShapepTDFirstDer  = jet.fJetShapepTDFirstDer;
    fJetShapepTDSecondDer = jet.fJetShapepTDSecondDer;
    fJetShapepTDFirstSub  = jet.fJetShapepTDFirstSub;
    fJetShapepTDSecondSub = jet.fJetShapepTDSecondSub;

    fJetShapeCircularityFirstDer  = jet.fJetShapeCircularityFirstDer;
    fJetShapeCircularitySecondDer = jet.fJetShapeCircularitySecondDer;
    fJetShapeCircularityFirstSub  = jet.fJetShapeCircularityFirstSub;
    fJetShapeCircularitySecondSub = jet.fJetShapeCircularitySecondSub;

 
    fJetShapeConstituentFirstDer  = jet.fJetShapeConstituentFirstDer;
    fJetShapeConstituentSecondDer = jet.fJetShapeConstituentSecondDer;
    fJetShapeConstituentFirstSub  = jet.fJetShapeConstituentFirstSub;
    fJetShapeConstituentSecondSub = jet.fJetShapeConstituentSecondSub;
    fJetShapeLeSubFirstDer  = jet.fJetShapeLeSubFirstDer;
    fJetShapeLeSubSecondDer = jet.fJetShapeLeSubSecondDer;
    fJetShapeLeSubFirstSub  = jet.fJetShapeLeSubFirstSub;
    fJetShapeLeSubSecondSub = jet.fJetShapeLeSubSecondSub;
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

  Printf("Jet pt=%.2f, eta=%.2f, phi=%.2f, area=%.2f, NEF=%.2f", fPt, fEta, fPhi, fArea, fNEF);
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
Double_t AliEmcalJet::DeltaR (const AliVParticle* part) const
{ // Helper function to calculate the distance between two jets or a jet and a particle
    Double_t dPhi = this->Phi() - part->Phi();
    Double_t dEta = this->Eta() - part->Eta();
    dPhi = TVector2::Phi_mpi_pi ( dPhi );

    return TMath::Sqrt ( dPhi * dPhi + dEta * dEta );
}


//__________________________________________________________________________________________________
std::vector<int> AliEmcalJet::SortConstituentsPt( TClonesArray *tracks ) const
{   //___________________________________________
    // Sorting by p_T (decreasing) jet constituents
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    typedef std::pair<Double_t, Int_t> ptidx_pair;

    // Create vector for Pt sorting
    std::vector<ptidx_pair> pair_list ;

    for ( Int_t i_entry = 0; i_entry < GetNumberOfTracks(); i_entry++ )
        {
        AliVParticle *track = TrackAt(i_entry, tracks);
        if (!track)
            {
            AliError(Form("Unable to find jet track %d in collection %s (pos in collection %d, max %d)", i_entry, tracks->GetName(), TrackAt(i_entry), tracks->GetEntriesFast()));
            continue;
            }

        pair_list.push_back( std::make_pair ( track->Pt(), i_entry ) );
        }

    std::stable_sort( pair_list.begin() , pair_list.end() , sort_descend() );

    // return an vector of indexes of constituents (sorted descending by pt)
    std::vector <int> index_sorted_list;

    for ( std::vector< std::pair<Double_t,Int_t> >::iterator it = pair_list.begin(); it != pair_list.end(); ++it)
        { index_sorted_list.push_back( (*it).second ); } // populating the return object with indexes of sorted tracks

    return index_sorted_list;
}

//__________________________________________________________________________________________________
AliVParticle* AliEmcalJet::GetLeadingTrack(TClonesArray *tracks) const
{
  AliVParticle* maxTrack = 0;
  for (Int_t i = 0; i < GetNumberOfTracks(); i++) {
    AliVParticle *track = TrackAt(i, tracks);
    if (!track) {
      AliError(Form("Unable to find jet track %d in collection %s (pos in collection %d, max %d)",
		    i,tracks->GetName(),TrackAt(i),tracks->GetEntriesFast()));
      continue;
    }
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
    if (!cluster) {
      AliError(Form("Unable to find jet cluster %d in collection %s (pos in collection %d, max %d)",
		    i,clusters->GetName(),ClusterAt(i),clusters->GetEntriesFast()));
      continue;
    }
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

//__________________________________________________________________________________________________
void AliEmcalJet::PrintGR() {
  for(Int_t i = 0; i<fGRNumerator.GetSize(); i++) {
    Printf("num[%d] = %f",i,fGRNumerator.At(i));
  }
}
