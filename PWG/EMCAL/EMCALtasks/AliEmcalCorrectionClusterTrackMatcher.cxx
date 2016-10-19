// AliEmcalCorrectionClusterTrackMatcher
//

#include "AliEmcalCorrectionClusterTrackMatcher.h"

#include <TH1.h>
#include <TList.h>

#include "AliClusterContainer.h"
#include "AliParticleContainer.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliVParticle.h"
#include "AliEmcalParticle.h"
#include "AliEMCALGeometry.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionClusterTrackMatcher);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionClusterTrackMatcher> AliEmcalCorrectionClusterTrackMatcher::reg("AliEmcalCorrectionClusterTrackMatcher");

//________________________________________________________________________
AliEmcalCorrectionClusterTrackMatcher::AliEmcalCorrectionClusterTrackMatcher() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionClusterTrackMatcher"),
  fPropDist(440),
  fDoPropagation(kFALSE),
  fAttemptProp(kTRUE),
  fAttemptPropMatch(kFALSE),
  fMaxDistance(0.1),
  fUpdateTracks(kTRUE),
  fUpdateClusters(kTRUE),
  fEmcalTracks(0),
  fEmcalClusters(0),
  fNEmcalTracks(0),
  fNEmcalClusters(0),
  fHistMatchEtaAll(0),
  fHistMatchPhiAll(0)
{
  // Default constructor
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  
  for(Int_t icent=0; icent<8; ++icent) {
    for(Int_t ipt=0; ipt<9; ++ipt) {
      for(Int_t ieta=0; ieta<2; ++ieta) {
        fHistMatchEta[icent][ipt][ieta] = 0;
        fHistMatchPhi[icent][ipt][ieta] = 0;
      }
    }
  }

}

//________________________________________________________________________
AliEmcalCorrectionClusterTrackMatcher::~AliEmcalCorrectionClusterTrackMatcher()
{
  // Destructor
}

//________________________________________________________________________
Bool_t AliEmcalCorrectionClusterTrackMatcher::Initialize()
{
  // Initialization
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::Initialize();
  // Do base class initializations and if it fails -> bail out
  //AliAnalysisTaskEmcal::ExecOnce();
  //if (!fInitialized) return;
  
  GetProperty("createHistos", fCreateHisto);
  
  GetProperty("maxDist", fMaxDistance);
  GetProperty("updateClusters", fUpdateClusters);
  GetProperty("updateTracks", fUpdateTracks);
  fDoPropagation = fEsdMode;
  
  return kTRUE;
}

//________________________________________________________________________
void AliEmcalCorrectionClusterTrackMatcher::UserCreateOutputObjects()
{   
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::UserCreateOutputObjects();

  fEmcalTracks = new TClonesArray("AliEmcalParticle");
  fEmcalTracks->SetName(Form("EmcalTracks_%s", fPartCont->GetArrayName().Data()));
  fEmcalClusters = new TClonesArray("AliEmcalParticle");
  fEmcalClusters->SetName(Form("EmcalClusters_%s", fClusCont->GetArrayName().Data()));
 
  // Create my user objects.
  if (fCreateHisto){
    fHistMatchEtaAll = new TH1F("fHistMatchEtaAll", "fHistMatchEtaAll", 400, -0.2, 0.2);
    fHistMatchPhiAll = new TH1F("fHistMatchPhiAll", "fHistMatchPhiAll", 400, -0.2, 0.2);
    fOutput->Add(fHistMatchEtaAll);
    fOutput->Add(fHistMatchPhiAll);
    
    const Int_t nCentChBins = fNcentBins * 2;
    for(Int_t icent=0; icent<nCentChBins; ++icent) {
      for(Int_t ipt=0; ipt<9; ++ipt) {
        for(Int_t ieta=0; ieta<2; ++ieta) {
          TString nameEta(Form("fHistMatchEta_%i_%i_%i",icent,ipt,ieta));
          fHistMatchEta[icent][ipt][ieta] = new TH1F(nameEta, nameEta, 400, -0.2, 0.2);
          fHistMatchEta[icent][ipt][ieta]->SetXTitle("#Delta#eta");
          TString namePhi(Form("fHistMatchPhi_%i_%i_%i",icent,ipt,ieta));
          fHistMatchPhi[icent][ipt][ieta] = new TH1F(namePhi, namePhi, 400, -0.2, 0.2);
          fHistMatchPhi[icent][ipt][ieta]->SetXTitle("#Delta#phi");
          fOutput->Add(fHistMatchEta[icent][ipt][ieta]);
          fOutput->Add(fHistMatchPhi[icent][ipt][ieta]);
        }
      }
    }
    fOutput->SetOwner(kTRUE);
  }
}

//________________________________________________________________________
Bool_t AliEmcalCorrectionClusterTrackMatcher::Run()
{
  // Run
  AliDebug(3, Form("%s", __PRETTY_FUNCTION__));
  AliEmcalCorrectionComponent::Run();

  // Run the matching.
  GenerateEmcalParticles();
  DoMatching();
  if (fUpdateTracks) UpdateTracks();
  if (fUpdateClusters) UpdateClusters();
  
  return kTRUE;
}

//________________________________________________________________________
Int_t AliEmcalCorrectionClusterTrackMatcher::GetMomBin(Double_t p) const
{
  // Get momenum bin.
  
  Int_t pbin=-1;
  if (p<0.5)
    pbin=0;
  else if (p>=0.5 && p<1.0)
    pbin=1;
  else if (p>=1.0 && p<1.5)
    pbin=2;
  else if (p>=1.5 && p<2.)
    pbin=3;
  else if (p>=2. && p<3.)
    pbin=4;
  else if (p>=3. && p<4.)
    pbin=5;
  else if (p>=4. && p<5.)
    pbin=6;
  else if (p>=5. && p<8.)
    pbin=7;
  else if (p>=8.)
    pbin=8;
  
  return pbin;
}

//________________________________________________________________________
void AliEmcalCorrectionClusterTrackMatcher::GenerateEmcalParticles()
{
  // Create AliEmcalParticle collections to handle the matching efficiently.
  // At the same time propagates tracks, if requested.

  fEmcalTracks->Delete();
  fEmcalClusters->Delete();
  
  fNEmcalTracks = 0;
  fNEmcalClusters = 0;
  
  AliVCluster* cluster = 0;
  AliVTrack* track = 0;
  
  fClusCont->ResetCurrentID();
  while ((cluster = static_cast<AliVCluster*>(fClusCont->GetNextAcceptCluster()))) {
    
    // Clears the matching info
    cluster->SetEmcCpvDistance(-1);
    cluster->SetTrackDistance(1024, 1024);
    AliAODCaloCluster *ac = dynamic_cast<AliAODCaloCluster*>(cluster);
    AliESDCaloCluster *ec = 0;
    if (ac) {
      const Int_t N = ac->GetNTracksMatched();
      AliDebug(2, TString::Format("Number of matched tracks: %d", N));
      for (Int_t i = N - 1; i >= 0; i--) {
        TObject *ptr = ac->GetTrackMatched(i);
        ac->RemoveTrackMatched(ptr);
        AliDebug(2, TString::Format("N tracks matched: %i of %i", ac->GetNTracksMatched(), N));
      }
    }
    else {
      ec = dynamic_cast<AliESDCaloCluster*>(cluster);
      TArrayI *arr = ec->GetTracksMatched();
      if (arr) arr->Set(0);
    }
    
    // Create AliEmcalParticle objects to handle the matching
    AliEmcalParticle* emcalCluster = new ((*fEmcalClusters)[fNEmcalClusters])
    AliEmcalParticle(cluster, fClusCont->GetCurrentID(), fVertex[0], fVertex[1], fVertex[2], AliVCluster::kNonLinCorr);
    emcalCluster->SetMatchedPtr(fEmcalTracks);
    
    fNEmcalClusters++;
  }
  
  fPartCont->ResetCurrentID();
  while ((track = static_cast<AliVTrack*>(fPartCont->GetNextAcceptParticle()))) {

    // Clears the matching info
    track->ResetStatus(AliVTrack::kEMCALmatch);
    track->SetEMCALcluster(-1);

    // Propagate tracks if requested
    Bool_t propthistrack = kFALSE;
    if (fDoPropagation) {
      propthistrack = kTRUE;
    }
    else if (!track->IsExtrapolatedToEMCAL()) {
      if (fAttemptProp) {
        propthistrack = kTRUE;
      }
      else if (fAttemptPropMatch && IsTrackInEmcalAcceptance(track)) {
        propthistrack = kTRUE;
      }
    }
    if (propthistrack) AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(track, fPropDist);
    
    // Create AliEmcalParticle objects to handle the matching
    AliEmcalParticle* emcalTrack = new ((*fEmcalTracks)[fNEmcalTracks]) AliEmcalParticle(track, fPartCont->GetCurrentID());
    emcalTrack->SetMatchedPtr(fEmcalClusters);
    
    AliDebug(2, Form("Now adding track (pT = %.3f, eta = %.3f, phi = %.3f)"
                     "Phi, Eta on EMCal = %.3f, %.3f",
                     emcalTrack->Pt(), emcalTrack->Eta(), emcalTrack->Phi(),
                     track->GetTrackEtaOnEMCal(), track->GetTrackPhiOnEMCal()));
    
    fNEmcalTracks++;
  }
}

//________________________________________________________________________
void AliEmcalCorrectionClusterTrackMatcher::DoMatching()
{
  // Set the links between tracks and clusters.
  const Double_t maxd2 = fMaxDistance*fMaxDistance;

  for (Int_t itrack = 0; itrack < fNEmcalTracks; itrack++) {
    AliEmcalParticle* emcalTrack = static_cast<AliEmcalParticle*>(fEmcalTracks->At(itrack));
    AliVTrack* track = emcalTrack->GetTrack();

    for (Int_t icluster = 0; icluster < fNEmcalClusters; icluster++) {
      AliEmcalParticle* emcalCluster = static_cast<AliEmcalParticle*>(fEmcalClusters->At(icluster));
      AliVCluster* cluster = emcalCluster->GetCluster();
      
      Double_t deta = 999;
      Double_t dphi = 999;
      GetEtaPhiDiff(track, cluster, dphi, deta);
      Double_t d2 = deta * deta + dphi * dphi;

      if (d2 > maxd2) continue;
      
      Double_t d = TMath::Sqrt(d2);
      emcalCluster->AddMatchedObj(itrack, d);
      emcalTrack->AddMatchedObj(icluster, d);
      AliDebug(2, Form("Now matching cluster E = %.3f, pT = %.3f, eta = %.3f, phi = %.3f "
                       "with track pT = %.3f, eta = %.3f, phi = %.3f"
                       "Track eta, phi on EMCal = %.3f, %.3f, d = %.3f",
                       cluster->GetNonLinCorrEnergy(), emcalCluster->Pt(), emcalCluster->Eta(), emcalCluster->Phi(),
                       emcalTrack->Pt(), emcalTrack->Eta(), emcalTrack->Phi(),
                       track->GetTrackEtaOnEMCal(), track->GetTrackPhiOnEMCal(), d));
      
      if (fCreateHisto) {
        Int_t mombin = GetMomBin(track->P());
        Int_t centbinch = fCentBin;
        if (track->Charge() < 0) centbinch += fNcentBins;
        Int_t etabin = 0;
        if(track->Eta() > 0) etabin = 1;

        fHistMatchEta[centbinch][mombin][etabin]->Fill(deta);
        fHistMatchPhi[centbinch][mombin][etabin]->Fill(dphi);
        fHistMatchEtaAll->Fill(deta);
        fHistMatchPhiAll->Fill(dphi);
      }
    }
  }
}

//________________________________________________________________________
void AliEmcalCorrectionClusterTrackMatcher::UpdateClusters()
{
  // Update clusters with matching info.
  
  for (Int_t icluster = 0; icluster < fNEmcalClusters; icluster++) {
    AliEmcalParticle* emcalCluster = static_cast<AliEmcalParticle*>(fEmcalClusters->At(icluster));
    const Int_t N = emcalCluster->GetNumberOfMatchedObj();
    AliVCluster* cluster = emcalCluster->GetCluster();
    AliDebug(3, Form("Cluster E = %.2f, eta = %.2f, phi = %.2f, Nmatch = %d", cluster->GetNonLinCorrEnergy(), emcalCluster->Eta(), emcalCluster->Phi(), N));
    
    if (N <= 0) continue;
    
    // Set the first match distance
    const UInt_t firstMatchId = emcalCluster->GetMatchedObjId();
    AliEmcalParticle* emcalTrackFirstMatch = static_cast<AliEmcalParticle*>(fEmcalTracks->At(firstMatchId));
    AliVTrack* trackFirstMatch = emcalTrackFirstMatch->GetTrack();
    Double_t deta = 999;
    Double_t dphi = 999;
    GetEtaPhiDiff(trackFirstMatch, cluster, dphi, deta);
    cluster->SetTrackDistance(dphi, deta);
    
    // Cast into ESD/AOD objects
    AliAODCaloCluster *ac = dynamic_cast<AliAODCaloCluster*>(cluster);
    AliESDCaloCluster *ec = 0;
    if (!ac) ec = dynamic_cast<AliESDCaloCluster*>(cluster);
    
    // Copy the matched tracks in the cluster. Note: different methods for ESD/AOD
    if (ac) {
      // TEMP
      /*std::cout << "ProcessID for current process: " << TProcessID::GetPID()->GetName() << "/" << TProcessID::GetPID()->GetTitle() << ". Memory Address: " << TProcessID::GetPID() << std::endl;
      std::cout << "N tracks matched: " << ac->GetNTracksMatched() << std::endl;
      std::cout << "ProcessID for cluster: " << TProcessID::GetProcessWithUID(ac)->GetName() << "/" << TProcessID::GetProcessWithUID(ac)->GetTitle() << ". Memory Address: " << TProcessID::GetProcessWithUID(ac) << std::endl;*/
      for (Int_t i=0; i < N; ++i) {
        Int_t id = emcalCluster->GetMatchedObjId(i);
        AliEmcalParticle* emcalTrack = static_cast<AliEmcalParticle*>(fEmcalTracks->At(id));
        
        AliDebug(3, Form("Pt = %.2f, eta = %.2f, phi = %.2f", emcalTrack->Pt(), emcalTrack->Eta(), emcalTrack->Phi()));
        // TEMP
        //std::cout << "ProcessID for emcal track: " << TProcessID::GetProcessWithUID(emcalTrack)->GetName() << "/" << TProcessID::GetProcessWithUID(emcalTrack)->GetTitle() << ". Memory Address: " << TProcessID::GetProcessWithUID(emcalTrack) << std::endl;
        
        TObject *obj = emcalTrack->GetTrack();
        // TEMP
        // Superceded by assigning a new process ID earlier
        //obj->SetBit(TObject::kIsReferenced);
        //std::cout << "ProcessID for track: " << TProcessID::GetProcessWithUID(obj)->GetName() << "/" << TProcessID::GetProcessWithUID(obj)->GetTitle() << ". Memory Address: " << TProcessID::GetProcessWithUID(obj) << std::endl;
        ac->AddTrackMatched(obj);
      }
    }
    else {
      TArrayI arr(N);
      for (Int_t i = 0; i < N; ++i) {
        Int_t id = emcalCluster->GetMatchedObjId(i);
        AliEmcalParticle* emcalTrack = static_cast<AliEmcalParticle*>(fEmcalTracks->At(id));
        arr.AddAt(emcalTrack->IdInCollection(), i);
      }
      ec->AddTracksMatched(arr);
    }
  }
}

//________________________________________________________________________
void AliEmcalCorrectionClusterTrackMatcher::UpdateTracks()
{
  // Update tracks with matching info.
  
  for (Int_t itrack = 0; itrack < fNEmcalTracks; itrack++) {
    AliEmcalParticle* emcalTrack = static_cast<AliEmcalParticle*>(fEmcalTracks->At(itrack));

    if (emcalTrack->GetNumberOfMatchedObj() <= 0) continue;
    
    AliEmcalParticle* emcalCluster = static_cast<AliEmcalParticle*>(fEmcalClusters->At(emcalTrack->GetMatchedObjId()));
    
    AliVTrack* track = emcalTrack->GetTrack();
    track->SetEMCALcluster(emcalCluster->IdInCollection());
    track->SetStatus(AliVTrack::kEMCALmatch);
  }
}

/**
 * Determines if a track is inside the EMCal acceptance, using \f$\eta\f$/\f$\phi\f$ at the vertex (no propagation).
 * Includes +/- edges. Useful to determine whether track propagation should be attempted.
 * @param[in] part Particle to check
 * @param[in] edges Size of the edges in \f$\phi\f$ excluded from the EMCAL acceptance
 * @return True if a particle is inside the EMCAL acceptance, false otherwise
 */
Bool_t AliEmcalCorrectionClusterTrackMatcher::IsTrackInEmcalAcceptance(AliVParticle* part, Double_t edges) const
{
  
  if (!fGeom) {
    AliWarning(Form("%s - AliAnalysisTaskEmcal::IsTrackInEmcalAcceptance - Geometry is not available!", GetName()));
    return kFALSE;
  }
  
  Double_t minPhi = fGeom->GetArm1PhiMin() - edges;
  Double_t maxPhi = fGeom->GetArm1PhiMax() + edges;
  
  if (part->Phi() > minPhi && part->Phi() < maxPhi) {
    return kTRUE;
  }
  else {
    return kFALSE;
  }
}
