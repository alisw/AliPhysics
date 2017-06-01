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
#include "AliMCEvent.h"

/// \cond CLASSIMP
ClassImp(AliEmcalCorrectionClusterTrackMatcher);
/// \endcond

// Actually registers the class with the base class
RegisterCorrectionComponent<AliEmcalCorrectionClusterTrackMatcher> AliEmcalCorrectionClusterTrackMatcher::reg("AliEmcalCorrectionClusterTrackMatcher");

/**
 * Default constructor
 */
AliEmcalCorrectionClusterTrackMatcher::AliEmcalCorrectionClusterTrackMatcher() :
  AliEmcalCorrectionComponent("AliEmcalCorrectionClusterTrackMatcher"),
  fPropDist(440),
  fDoPropagation(kFALSE),
  fAttemptProp(kTRUE),
  fAttemptPropMatch(kFALSE),
  fMaxDistance(0.1),
  fUsePIDmass(kTRUE),
  fUseDCA(kTRUE),
  fUpdateTracks(kTRUE),
  fUpdateClusters(kTRUE),
  fClusterContainerIndexMap(),
  fParticleContainerIndexMap(),
  fEmcalTracks(0),
  fEmcalClusters(0),
  fNEmcalTracks(0),
  fNEmcalClusters(0),
  fHistMatchEtaAll(0),
  fHistMatchPhiAll(0),
  fMCGenerToAcceptForTrack(1),
  fNMCGenerToAccept(0)
{
  for(Int_t icent=0; icent<8; ++icent) {
    for(Int_t ipt=0; ipt<9; ++ipt) {
      for(Int_t ieta=0; ieta<2; ++ieta) {
        fHistMatchEta[icent][ipt][ieta] = 0;
        fHistMatchPhi[icent][ipt][ieta] = 0;
      }
    }
  }
  
  for(Int_t j = 0; j <  5;    j++)  fMCGenerToAccept[j] =  "";
}

/**
 * Destructor
 */
AliEmcalCorrectionClusterTrackMatcher::~AliEmcalCorrectionClusterTrackMatcher()
{
}

/**
 * Initialize and configure the component.
 */
Bool_t AliEmcalCorrectionClusterTrackMatcher::Initialize()
{
  // Initialization
  AliEmcalCorrectionComponent::Initialize();
  
  GetProperty("createHistos", fCreateHisto);
  GetProperty("usePIDmass", fUsePIDmass);
  GetProperty("useDCA", fUseDCA);
  GetProperty("maxDist", fMaxDistance);
  GetProperty("updateClusters", fUpdateClusters);
  GetProperty("updateTracks", fUpdateTracks);
  fDoPropagation = fEsdMode;
  
  Bool_t enableFracEMCRecalc = kFALSE;
  GetProperty("enableFracEMCRecalc", enableFracEMCRecalc);
  Int_t removeNMCGenerators = 0;
  GetProperty("removeNMCGenerators", removeNMCGenerators);
  GetProperty("enableMCGenRemovTrack", fMCGenerToAcceptForTrack);
  std::string removeMcGen1 = "";
  GetProperty("removeMCGen1", removeMcGen1);
  TString removeMCGen1 = removeMcGen1.c_str();
  std::string removeMcGen2 = "";
  GetProperty("removeMCGen2", removeMcGen2);
  TString removeMCGen2 = removeMcGen2.c_str();
  
  if(enableFracEMCRecalc){
    if(removeNMCGenerators > 0)
    {
      printf("\t gen1 <%s>, gen2 <%s>, remove tracks %d\n",removeMCGen1.Data(),removeMCGen2.Data(),fMCGenerToAcceptForTrack);
      SetNumberOfMCGeneratorsToAccept(removeNMCGenerators) ;
      SetNameOfMCGeneratorsToAccept(0,removeMCGen1);
      SetNameOfMCGeneratorsToAccept(1,removeMCGen2);
    }
  }

  return kTRUE;
}

/**
 * Create run-independent objects for output. Called before running over events.
 */
void AliEmcalCorrectionClusterTrackMatcher::UserCreateOutputObjects()
{   
  AliEmcalCorrectionComponent::UserCreateOutputObjects();

  // Determine all particle container array names for naming the AliEmcalParticle array
  std::string particleContainerNames = "";
  bool firstLoop = true;
  AliParticleContainer * partCont = 0;
  TIter nextPartCont(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(nextPartCont()))) {
    if (firstLoop != true) {
      particleContainerNames += "_";
    }
    else {
      firstLoop = false;
    }
    particleContainerNames += partCont->GetArrayName();
  }
  // Determine all cluster container array names for naming the AliEmcalParticle array
  std::string clusterContainerNames = "";
  firstLoop = true;
  AliClusterContainer * clusCont = 0;
  TIter nextClusCont(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(nextClusCont()))) {
    if (firstLoop != true) {
      clusterContainerNames += "_";
    }
    else {
      firstLoop = false;
    }
    clusterContainerNames += clusCont->GetArrayName();
  }

  fEmcalTracks = new TClonesArray("AliEmcalParticle");
  fEmcalTracks->SetName(Form("EmcalTracks_%s", particleContainerNames.c_str()));
  fEmcalClusters = new TClonesArray("AliEmcalParticle");
  fEmcalClusters->SetName(Form("EmcalClusters_%s", clusterContainerNames.c_str()));
 
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

/**
 * Called before the first event to initialize the correction.
 */
void AliEmcalCorrectionClusterTrackMatcher::ExecOnce()
{
  fClusterContainerIndexMap.CopyMappingFrom(AliClusterContainer::GetEmcalContainerIndexMap(), fClusterCollArray);
  fParticleContainerIndexMap.CopyMappingFrom(AliParticleContainer::GetEmcalContainerIndexMap(), fParticleCollArray);
}

/**
 * Called for each event to process the event data.
 */
Bool_t AliEmcalCorrectionClusterTrackMatcher::Run()
{
  AliEmcalCorrectionComponent::Run();
  
  CheckIfRunChanged();

  // Run the matching.
  GenerateEmcalParticles();
  DoMatching();
  if (fUpdateTracks) UpdateTracks();
  if (fUpdateClusters) UpdateClusters();
  
  return kTRUE;
}

/**
 * Get momentum bin.
 */
Int_t AliEmcalCorrectionClusterTrackMatcher::GetMomBin(Double_t p) const
{
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

/**
 * Create AliEmcalParticle collections to handle the matching efficiently.
 * At the same time propagates tracks, if requested.
 */
void AliEmcalCorrectionClusterTrackMatcher::GenerateEmcalParticles()
{
  fEmcalTracks->Delete();
  fEmcalClusters->Delete();

  fNEmcalTracks = 0;
  fNEmcalClusters = 0;

  AliVCluster* cluster = 0;
  AliVTrack* track = 0;
  
  AliClusterContainer * clusCont = 0;
  TIter nextClusCont(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(nextClusCont()))) {
    auto clusItCont = clusCont->accepted_momentum();
    for (AliClusterIterableMomentumContainer::iterator clusIterator = clusItCont.begin(); clusIterator != clusItCont.end(); ++clusIterator) {
      cluster = static_cast<AliVCluster *>(clusIterator->second);

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
      AliEmcalParticle(cluster, fClusterContainerIndexMap.GlobalIndexFromLocalIndex(clusCont, clusIterator.current_index()), fVertex[0], fVertex[1], fVertex[2], AliVCluster::kNonLinCorr);
      emcalCluster->SetMatchedPtr(fEmcalTracks);

      fNEmcalClusters++;
    }
  }
  
  Double_t mass;
  if (fUsePIDmass) {
    // use PID-based mass, and fUseDCA to determine starting point of propagation
    mass = -1;
  }
  else {
    // use pion mass, and fUseDCA to determine starting point of propagation
    mass = 0.1396;
  }

  AliParticleContainer * partCont = 0;
  TIter nextPartCont(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(nextPartCont()))) {
    auto partItCont = partCont->accepted_momentum();
    for (AliParticleIterableMomentumContainer::iterator partIterator = partItCont.begin(); partIterator != partItCont.end(); ++partIterator) {
      track = static_cast<AliVTrack *>(partIterator->second);

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
      if (propthistrack) {
        
        // Check if track comes from a particular MC generator, do not include it if it is not a selected one
        Int_t mcLabel = TMath::Abs(track->GetLabel());
        TString genName;
        if( fMCEvent && fMCGenerToAcceptForTrack && fNMCGenerToAccept > 0 )
        {
          fMCEvent->GetCocktailGenerator(mcLabel,genName);
          
          Bool_t generOK = kFALSE;
          for(Int_t ig = 0; ig < fNMCGenerToAccept; ig++)
          {
            if ( genName.Contains(fMCGenerToAccept[ig]) ) generOK = kTRUE;
          }
          
          if ( !generOK ) continue;
        }
        
        // Propagate the track
        AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(track, fPropDist, mass, 20, 0.35, kFALSE, fUseDCA);
      }

      // Reset properties of the track to fix TRefArray errors which occur when AddTrackMatched(obj) is called.
      // This particular combination is from AODTrackFilterTask
      // Resetting just the TProcessID of the track is not sufficient!
      // It is included here because multiple track collections can come from different files with different
      // Process IDs (for example, when embedding). As long as the updated tracks are not stored in a file,
      // (normally they are not) then resetting these values won't adversely change the behavior of the code
      // (except for fixing the TRefArray errors).
      track->SetUniqueID(0);
      track->ResetBit(TObject::kHasUUID);
      track->ResetBit(TObject::kIsReferenced);

      // Create AliEmcalParticle objects to handle the matching
      AliEmcalParticle* emcalTrack = new ((*fEmcalTracks)[fNEmcalTracks]) AliEmcalParticle(track, fParticleContainerIndexMap.GlobalIndexFromLocalIndex(partCont, partIterator.current_index()));
      emcalTrack->SetMatchedPtr(fEmcalClusters);
      
      AliDebug(2, Form("Now adding track %i (pT = %.3f, eta = %.3f, phi = %.3f)"
                       "Phi, Eta on EMCal = %.3f, %.3f",
                       emcalTrack->IdInCollection(), emcalTrack->Pt(), emcalTrack->Eta(), emcalTrack->Phi(),
                       track->GetTrackEtaOnEMCal(), track->GetTrackPhiOnEMCal()));
      
      fNEmcalTracks++;
    }
  }
}

/**
 * Set the links between tracks and clusters.
 */
void AliEmcalCorrectionClusterTrackMatcher::DoMatching()
{
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

/**
 * Update clusters with matching info.
 */
void AliEmcalCorrectionClusterTrackMatcher::UpdateClusters()
{
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
      for (Int_t i=0; i < N; ++i) {
        Int_t id = emcalCluster->GetMatchedObjId(i);
        AliEmcalParticle* emcalTrack = static_cast<AliEmcalParticle*>(fEmcalTracks->At(id));
        
        AliDebug(3, Form("Pt = %.2f, eta = %.2f, phi = %.2f", emcalTrack->Pt(), emcalTrack->Eta(), emcalTrack->Phi()));
        
        TObject *obj = emcalTrack->GetTrack();
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

/**
 * Update tracks with matching info.
 */
void AliEmcalCorrectionClusterTrackMatcher::UpdateTracks()
{
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
