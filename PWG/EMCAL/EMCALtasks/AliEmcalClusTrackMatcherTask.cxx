//
// Track/cluster matcher
// 
// Author: C.Loizides, S.Aiola

#include "AliEmcalClusTrackMatcherTask.h"

#include <TClonesArray.h>
#include <TClass.h>

#include <AliAODCaloCluster.h>
#include <AliESDCaloCluster.h>
#include <AliLog.h>
#include <AliVCluster.h>
#include <AliVTrack.h>
#include <AliEMCALRecoUtils.h>

#include "AliEmcalParticle.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

ClassImp(AliEmcalClusTrackMatcherTask)

//________________________________________________________________________
AliEmcalClusTrackMatcherTask::AliEmcalClusTrackMatcherTask() : 
  AliAnalysisTaskEmcal("AliEmcalClusTrackMatcherTask", kFALSE),
  fPropDist(440),
  fDoPropagation(kFALSE),
  fAttemptProp(kTRUE),
  fAttemptPropMatch(kFALSE),
  fMaxDistance(0.1),
  fAttachEmcalParticles(kFALSE),
  fUpdateTracks(kTRUE),
  fUpdateClusters(kTRUE),
  fEmcalTracks(0),
  fEmcalClusters(0),
  fNEmcalTracks(0),
  fNEmcalClusters(0),
  fHistMatchEtaAll(0),
  fHistMatchPhiAll(0)
{
  // Constructor.

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
AliEmcalClusTrackMatcherTask::AliEmcalClusTrackMatcherTask(const char *name, Bool_t histo) : 
  AliAnalysisTaskEmcal(name, histo),
  fPropDist(440),
  fDoPropagation(kFALSE),
  fAttemptProp(kTRUE),
  fAttemptPropMatch(kFALSE),
  fMaxDistance(0.1),
  fAttachEmcalParticles(kFALSE),
  fUpdateTracks(kTRUE),
  fUpdateClusters(kTRUE),
  fEmcalTracks(0),
  fEmcalClusters(0),
  fNEmcalTracks(0),
  fNEmcalClusters(0),
  fHistMatchEtaAll(0),
  fHistMatchPhiAll(0)
{
  // Standard constructor.

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
AliEmcalClusTrackMatcherTask::~AliEmcalClusTrackMatcherTask()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalClusTrackMatcherTask::ExecOnce()
{
  // Initialize the analysis.

  AliParticleContainer* tracks = GetParticleContainer(0);
  if (tracks) {
    TClass trackClass(tracks->GetClassName());
    if (!trackClass.InheritsFrom("AliVTrack")) {
      tracks->SetClassName("AliVTrack"); // enforce only AliVTrack and derived classes
    }
  }

  AliClusterContainer* clusters = GetClusterContainer(0);

  AliAnalysisTaskEmcal::ExecOnce();
  if (!fInitialized) return;

  TString emcalTracksName(Form("EmcalTracks_%s", tracks->GetArrayName().Data()));
  TString emcalClustersName(Form("EmcalClusters_%s", clusters->GetArrayName().Data()));

  fEmcalTracks = new TClonesArray("AliEmcalParticle");
  fEmcalTracks->SetName(emcalTracksName);
  fEmcalClusters = new TClonesArray("AliEmcalParticle");
  fEmcalClusters->SetName(emcalClustersName);

  if (fAttachEmcalParticles) {
    AddObjectToEvent(fEmcalTracks);
    AddObjectToEvent(fEmcalClusters);
  }
}

//________________________________________________________________________
void AliEmcalClusTrackMatcherTask::UserCreateOutputObjects()
{
  // Create my user objects.

  if (!fCreateHisto) return;

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  const Int_t nCentChBins = fNcentBins * 2;

  fHistMatchEtaAll = new TH1F("fHistMatchEtaAll", "fHistMatchEtaAll", 400, -0.2, 0.2);
  fHistMatchPhiAll = new TH1F("fHistMatchPhiAll", "fHistMatchPhiAll", 400, -0.2, 0.2);
  fOutput->Add(fHistMatchEtaAll);
  fOutput->Add(fHistMatchPhiAll);

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

  PostData(1, fOutput);
}

//________________________________________________________________________
Int_t AliEmcalClusTrackMatcherTask::GetMomBin(Double_t p) const
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
Bool_t AliEmcalClusTrackMatcherTask::Run() 
{
  // Run the matching.

  GenerateEmcalParticles();
  DoMatching();
  if (fUpdateTracks) UpdateTracks();
  if (fUpdateClusters) UpdateClusters();

  return kTRUE;
}

//________________________________________________________________________
void AliEmcalClusTrackMatcherTask::GenerateEmcalParticles()
{
  // Create AliEmcalParticle collections to handle the matching efficiently.
  // At the same time propagates tracks, if requested.

  AliParticleContainer *tracks = GetParticleContainer(0);
  AliClusterContainer *clusters = GetClusterContainer(0);

  fEmcalTracks->Delete();
  fEmcalClusters->Delete();

  fNEmcalTracks = 0;
  fNEmcalClusters = 0;

  AliVCluster* cluster = 0;
  AliVTrack* track = 0;

  clusters->ResetCurrentID();
  while ((cluster = static_cast<AliVCluster*>(clusters->GetNextAcceptCluster()))) {

    // Clears the matching info
    cluster->SetEmcCpvDistance(-1);
    cluster->SetTrackDistance(1024, 1024);
    AliAODCaloCluster *ac = dynamic_cast<AliAODCaloCluster*>(cluster);
    AliESDCaloCluster *ec = 0;
    if (ac) {
      const Int_t N = ac->GetNTracksMatched();
      for (Int_t i = N - 1; i >= 0; i--) {
        TObject *ptr = ac->GetTrackMatched(i);
        ac->RemoveTrackMatched(ptr);
      }
    }
    else {
      ec = dynamic_cast<AliESDCaloCluster*>(cluster);
      TArrayI *arr = ec->GetTracksMatched(); 
      if (arr) arr->Set(0);
    }

    // Create AliEmcalParticle objects to handle the matching
    AliEmcalParticle* emcalCluster = new ((*fEmcalClusters)[fNEmcalClusters])
          AliEmcalParticle(cluster, clusters->GetCurrentID(), fVertex[0], fVertex[1], fVertex[2], AliVCluster::kNonLinCorr);
    emcalCluster->SetMatchedPtr(fEmcalTracks);

    fNEmcalClusters++;
  }

  tracks->ResetCurrentID();
  while ((track = static_cast<AliVTrack*>(tracks->GetNextAcceptParticle()))) {

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
    AliEmcalParticle* emcalTrack = new ((*fEmcalTracks)[fNEmcalTracks]) AliEmcalParticle(track, tracks->GetCurrentID());
    emcalTrack->SetMatchedPtr(fEmcalClusters);

    AliDebug(2, Form("Now adding track (pT = %.3f, eta = %.3f, phi = %.3f)"
        "Phi, Eta on EMCal = %.3f, %.3f",
        emcalTrack->Pt(), emcalTrack->Eta(), emcalTrack->Phi(),
        track->GetTrackEtaOnEMCal(), track->GetTrackPhiOnEMCal()));

    fNEmcalTracks++;
  }
}

//________________________________________________________________________
void AliEmcalClusTrackMatcherTask::DoMatching() 
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
void AliEmcalClusTrackMatcherTask::UpdateClusters() 
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

//________________________________________________________________________
void AliEmcalClusTrackMatcherTask::UpdateTracks() 
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
