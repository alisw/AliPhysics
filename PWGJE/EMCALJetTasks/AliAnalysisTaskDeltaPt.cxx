// $Id$
//
// Jet deltaPt task.
//
// Author: S.Aiola

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"

#include "AliAnalysisTaskDeltaPt.h"

ClassImp(AliAnalysisTaskDeltaPt)

//________________________________________________________________________
AliAnalysisTaskDeltaPt::AliAnalysisTaskDeltaPt() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskDeltaPt", kTRUE),
  fMCJetPtThreshold(1),
  fMinRC2LJ(-1),
  fEmbJetsName(""),
  fEmbTracksName(""),
  fEmbCaloName(""),
  fRandTracksName("TracksRandomized"),
  fRandCaloName("CaloClustersRandomized"),
  fRCperEvent(-1),
  fEmbJets(0),
  fEmbTracks(0),
  fEmbCaloClusters(0),
  fRandTracks(0),
  fRandCaloClusters(0),
  fEmbeddedClusterNIds(0),
  fEmbeddedTrackNIds(0),
  fHistRCPhiEta(0),
  fHistRCPtExLJVSDPhiLJ(0),
  fHistEmbJetPhiEta(0),
  fHistEmbPartPhiEta(0)
{
  // Default constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistRCPt[i] = 0;
    fHistRCPtExLJ[i] = 0;
    fHistRCPtRand[i] = 0; 
    fHistRhoVSRCPt[i] = 0;
    fHistDeltaPtRC[i] = 0;
    fHistDeltaPtRCExLJ[i] = 0;
    fHistDeltaPtRCRand[i] = 0;
    fHistEmbNotFoundPhiEta[i] = 0;
    fHistEmbJetsPtArea[i] = 0;
    fHistEmbJetsCorrPtArea[i] = 0;
    fHistEmbPartPt[i] = 0;
    fHistDistEmbPartJetAxis[i] = 0;
    fHistEmbBkgArea[i] = 0;
    fHistRhoVSEmbBkg[i] = 0;
    fHistDeltaPtEmbArea[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskDeltaPt::AliAnalysisTaskDeltaPt(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fMCJetPtThreshold(1),
  fMinRC2LJ(-1),
  fEmbJetsName(""),
  fEmbTracksName(""),
  fEmbCaloName(""),
  fRandTracksName("TracksRandomized"),
  fRandCaloName("CaloClustersRandomized"),
  fRCperEvent(-1),
  fEmbJets(0),
  fEmbTracks(0),
  fEmbCaloClusters(0),
  fRandTracks(0),
  fRandCaloClusters(0),
  fEmbeddedClusterNIds(0),
  fEmbeddedTrackNIds(0),
  fHistRCPhiEta(0),
  fHistRCPtExLJVSDPhiLJ(0),
  fHistEmbJetPhiEta(0),
  fHistEmbPartPhiEta(0)
{
  // Standard constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistRCPt[i] = 0;
    fHistRCPtExLJ[i] = 0;
    fHistRCPtRand[i] = 0; 
    fHistRhoVSRCPt[i] = 0;
    fHistDeltaPtRC[i] = 0;
    fHistDeltaPtRCExLJ[i] = 0;
    fHistDeltaPtRCRand[i] = 0;
    fHistEmbNotFoundPhiEta[i] = 0;
    fHistEmbJetsPtArea[i] = 0;
    fHistEmbJetsCorrPtArea[i] = 0;
    fHistEmbPartPt[i] = 0;
    fHistDistEmbPartJetAxis[i] = 0;
    fHistEmbBkgArea[i] = 0;
    fHistRhoVSEmbBkg[i] = 0;
    fHistDeltaPtEmbArea[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskDeltaPt::~AliAnalysisTaskDeltaPt()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskDeltaPt::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fHistRCPhiEta = new TH2F("fHistRCPhiEta","Phi-Eta distribution of rigid cones", 50, -1, 1, 101, 0, TMath::Pi() * 2.02);
  fHistRCPhiEta->GetXaxis()->SetTitle("#eta");
  fHistRCPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistRCPhiEta);

  fHistRCPtExLJVSDPhiLJ = new TH2F("fHistRCPtExLJVSDPhiLJ","fHistRCPtExLJVSDPhiLJ", fNbins, fMinBinPt, fMaxBinPt, 128, -1.6, 4.8);
  fHistRCPtExLJVSDPhiLJ->GetXaxis()->SetTitle("rigid cone #it{p}_{T} (GeV/#it{c})");
  fHistRCPtExLJVSDPhiLJ->GetYaxis()->SetTitle("#Delta#phi");
  fOutput->Add(fHistRCPtExLJVSDPhiLJ);

  if (!fJetsName.IsNull()) {
    fHistEmbJetPhiEta = new TH2F("fHistEmbJetPhiEta","Phi-Eta distribution of embedded jets", 50, -1, 1, 101, 0, TMath::Pi() * 2.02);
    fHistEmbJetPhiEta->GetXaxis()->SetTitle("#eta");
    fHistEmbJetPhiEta->GetYaxis()->SetTitle("#phi");
    fOutput->Add(fHistEmbJetPhiEta);
    
    fHistEmbPartPhiEta = new TH2F("fHistEmbPartPhiEta","Phi-Eta distribution of embedded particles", 50, -1, 1, 101, 0, TMath::Pi() * 2.02);
    fHistEmbPartPhiEta->GetXaxis()->SetTitle("#eta");
    fHistEmbPartPhiEta->GetYaxis()->SetTitle("#phi");
    fOutput->Add(fHistEmbPartPhiEta);
  }

  TString histname;

  for (Int_t i = 0; i < 4; i++) {
    histname = "fHistRCPt_";
    histname += i;
    fHistRCPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
    fHistRCPt[i]->GetXaxis()->SetTitle("rigid cone #it{p}_{T} (GeV/#it{c})");
    fHistRCPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistRCPt[i]);

    histname = "fHistRCPtExLJ_";
    histname += i;
    fHistRCPtExLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
    fHistRCPtExLJ[i]->GetXaxis()->SetTitle("rigid cone #it{p}_{T}^{RC} (GeV/#it{c})");
    fHistRCPtExLJ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistRCPtExLJ[i]);

    histname = "fHistRCPtRand_";
    histname += i;
    fHistRCPtRand[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
    fHistRCPtRand[i]->GetXaxis()->SetTitle("rigid cone #it{p}_{T}^{RC} (GeV/#it{c})");
    fHistRCPtRand[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistRCPtRand[i]);

    histname = "fHistRhoVSRCPt_";
    histname += i;
    fHistRhoVSRCPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
    fHistRhoVSRCPt[i]->GetXaxis()->SetTitle("A#rho (GeV/#it{c})");
    fHistRhoVSRCPt[i]->GetYaxis()->SetTitle("rigid cone #it{p}_{T} (GeV/#it{c})");
    fOutput->Add(fHistRhoVSRCPt[i]);

    histname = "fHistDeltaPtRC_";
    histname += i;
    fHistDeltaPtRC[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
    fHistDeltaPtRC[i]->GetXaxis()->SetTitle("#delta#it{p}_{T}^{RC} (GeV/#it{c})");
    fHistDeltaPtRC[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtRC[i]);

    histname = "fHistDeltaPtRCExLJ_";
    histname += i;
    fHistDeltaPtRCExLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
    fHistDeltaPtRCExLJ[i]->GetXaxis()->SetTitle("#delta#it{p}_{T}^{RC} (GeV/#it{c})");
    fHistDeltaPtRCExLJ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtRCExLJ[i]);

    histname = "fHistDeltaPtRCRand_";
    histname += i;
    fHistDeltaPtRCRand[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
    fHistDeltaPtRCRand[i]->GetXaxis()->SetTitle("#delta#it{p}_{T}^{RC} (GeV/#it{c})");
    fHistDeltaPtRCRand[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtRCRand[i]);

    if (!fEmbJetsName.IsNull()) {
      histname = "fHistEmbJetsPtArea_";
      histname += i;
      fHistEmbJetsPtArea[i] = new TH2F(histname.Data(), histname.Data(), 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
      fHistEmbJetsPtArea[i]->GetXaxis()->SetTitle("area");
      fHistEmbJetsPtArea[i]->GetYaxis()->SetTitle("embedded jet #it{p}_{T}^{raw} (GeV/#it{c})");
      fOutput->Add(fHistEmbJetsPtArea[i]);

      histname = "fHistEmbJetsCorrPtArea_";
      histname += i;
      fHistEmbJetsCorrPtArea[i] = new TH2F(histname.Data(), histname.Data(), 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins * 2, -fMaxBinPt, fMaxBinPt);
      fHistEmbJetsCorrPtArea[i]->GetXaxis()->SetTitle("area");
      fHistEmbJetsCorrPtArea[i]->GetYaxis()->SetTitle("embedded jet #it{p}_{T}^{corr} (GeV/#it{c})");
      fOutput->Add(fHistEmbJetsCorrPtArea[i]);

      histname = "fHistEmbPartPt_";
      histname += i;
      fHistEmbPartPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
      fHistEmbPartPt[i]->GetXaxis()->SetTitle("embedded particle #it{p}_{T}^{emb} (GeV/#it{c})");
      fHistEmbPartPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistEmbPartPt[i]);

      histname = "fHistDistEmbPartJetAxis_";
      histname += i;
      fHistDistEmbPartJetAxis[i] = new TH1F(histname.Data(), histname.Data(), 50, 0, 0.5);
      fHistDistEmbPartJetAxis[i]->GetXaxis()->SetTitle("distance");
      fHistDistEmbPartJetAxis[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistDistEmbPartJetAxis[i]);

      histname = "fHistEmbNotFoundPhiEta_";
      histname += i;
      fHistEmbNotFoundPhiEta[i] = new TH2F(histname.Data(), histname.Data(), 40, -2, 2, 64, 0, 6.4);
      fHistEmbNotFoundPhiEta[i]->GetXaxis()->SetTitle("#eta");
      fHistEmbNotFoundPhiEta[i]->GetYaxis()->SetTitle("#phi");
      fOutput->Add(fHistEmbNotFoundPhiEta[i]);

      histname = "fHistEmbBkgArea_";
      histname += i;
      fHistEmbBkgArea[i] = new TH2F(histname.Data(), histname.Data(), 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
      fHistEmbBkgArea[i]->GetXaxis()->SetTitle("area");
      fHistEmbBkgArea[i]->GetYaxis()->SetTitle("background of embedded track (GeV/#it{c})");
      fOutput->Add(fHistEmbBkgArea[i]);

      histname = "fHistRhoVSEmbBkg_";
      histname += i;
      fHistRhoVSEmbBkg[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
      fHistRhoVSEmbBkg[i]->GetXaxis()->SetTitle("A#rho (GeV/#it{c})");
      fHistRhoVSEmbBkg[i]->GetYaxis()->SetTitle("background of embedded track (GeV/#it{c})");
      fOutput->Add(fHistRhoVSEmbBkg[i]);
      
      histname = "fHistDeltaPtEmbArea_";
      histname += i;
      fHistDeltaPtEmbArea[i] = new TH2F(histname.Data(), histname.Data(), 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins * 2, -fMaxBinPt, fMaxBinPt);
      fHistDeltaPtEmbArea[i]->GetXaxis()->SetTitle("area");
      fHistDeltaPtEmbArea[i]->GetYaxis()->SetTitle("#delta#it{p}_{T}^{emb} (GeV/#it{c})");
      fOutput->Add(fHistDeltaPtEmbArea[i]);
    }
  }

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
Bool_t AliAnalysisTaskDeltaPt::FillHistograms()
{
  // Fill histograms.

  Int_t *sortedJets = GetSortedArray(fJets);
  
  AliEmcalJet* jet = 0;

  if (sortedJets && sortedJets[0] > 0) 
    jet = static_cast<AliEmcalJet*>(fJets->At(sortedJets[0]));

  // ************
  // Random cones
  // _________________________________
  
  const Float_t rcArea = fJetRadius * fJetRadius * TMath::Pi();

  for (Int_t i = 0; i < fRCperEvent; i++) {
    // Simple random cones
    Float_t RCpt = 0;
    Float_t RCeta = 0;
    Float_t RCphi = 0;
    GetRandomCone(RCpt, RCeta, RCphi, 0);
    if (RCpt > 0) {
      fHistRCPhiEta->Fill(RCeta, RCphi);
      fHistRhoVSRCPt[fCentBin]->Fill(fRhoVal * rcArea, RCpt);

      fHistRCPt[fCentBin]->Fill(RCpt);
      fHistDeltaPtRC[fCentBin]->Fill(RCpt - rcArea * fRhoVal);
    }
  
    // Random cones far from leading jet
    RCpt = 0;
    RCeta = 0;
    RCphi = 0;
    GetRandomCone(RCpt, RCeta, RCphi, jet);
    if (RCpt > 0) {
      if (jet) {
	Float_t dphi = RCphi - jet->Phi();
	if (dphi > 4.8) dphi -= TMath::Pi() * 2;
	if (dphi < -1.6) dphi += TMath::Pi() * 2; 
	fHistRCPtExLJVSDPhiLJ->Fill(RCpt, dphi);
      }
      fHistRCPtExLJ[fCentBin]->Fill(RCpt);
      fHistDeltaPtRCExLJ[fCentBin]->Fill(RCpt - rcArea * fRhoVal);
    }
    
    // Random cones with randomized particles
    RCpt = 0;
    RCeta = 0;
    RCphi = 0;
    GetRandomCone(RCpt, RCeta, RCphi, 0, fRandTracks, fRandCaloClusters);
    if (RCpt > 0) {
      fHistRCPtRand[fCentBin]->Fill(RCpt);
      fHistDeltaPtRCRand[fCentBin]->Fill(RCpt - rcArea * fRhoVal);
    }  
  }
  
  // ************
  // Embedding
  // _________________________________

  if (!fEmbJets)
    return kTRUE;

  AliEmcalJet *embJet = NextEmbeddedJet(0);

  Int_t countEmbJets = 0;

  while (embJet != 0) {

    countEmbJets++;

    Double_t maxClusterPt = 0;
    Double_t maxClusterEta = 0;
    Double_t maxClusterPhi = 0;

    Double_t maxTrackPt = 0;
    Double_t maxTrackEta = 0;
    Double_t maxTrackPhi = 0;

    Double_t probePt = 0;
    Double_t probeEta = 0;
    Double_t probePhi = 0;

    AliVCluster *cluster = embJet->GetLeadingCluster(fEmbCaloClusters);
    if (cluster) {
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);

      maxClusterEta = nPart.Eta();
      maxClusterPhi = nPart.Phi();
      maxClusterPt = nPart.Pt();
    }

    AliVParticle *track = embJet->GetLeadingTrack(fEmbTracks);
    if (track) {
      maxTrackEta = track->Eta();
      maxTrackPhi = track->Phi();
      maxTrackPt = track->Pt();
    }

    if (!track && !cluster) {
      AliWarning(Form("%s - Embedded jet found but no leading particle was found (?) !", GetName()));
      return kTRUE;
    }

    if (maxTrackPt > maxClusterPt) {
      probePt = maxTrackPt;
      probeEta = maxTrackEta;
      probePhi = maxTrackPhi;
    }
    else {
      probePt = maxClusterPt;
      probeEta = maxClusterEta;
      probePhi = maxClusterPhi;
    }

    Double_t distProbeJet = TMath::Sqrt((embJet->Eta() - probeEta) * (embJet->Eta() - probeEta) + (embJet->Phi() - probePhi) * (embJet->Phi() - probePhi));

    fHistEmbPartPt[fCentBin]->Fill(probePt);
    fHistEmbPartPhiEta->Fill(probeEta, probePhi);
    fHistDistEmbPartJetAxis[fCentBin]->Fill(distProbeJet);

    fHistEmbJetsPtArea[fCentBin]->Fill(embJet->Area(), embJet->Pt());
    fHistEmbJetsCorrPtArea[fCentBin]->Fill(embJet->Area(), embJet->Pt() - fRhoVal * embJet->Area());
    fHistEmbJetPhiEta->Fill(embJet->Eta(), embJet->Phi());

    fHistEmbBkgArea[fCentBin]->Fill(embJet->Area(), embJet->Pt() - embJet->MCPt());
    fHistRhoVSEmbBkg[fCentBin]->Fill(fRhoVal * embJet->Area(), embJet->Pt() - embJet->MCPt());
    fHistDeltaPtEmbArea[fCentBin]->Fill(embJet->Area(), embJet->Pt() - embJet->Area() * fRhoVal - embJet->MCPt());

    embJet = NextEmbeddedJet();
  }

  if (countEmbJets==0) {
    if (fEmbTracks) {
      DoEmbTrackLoop();
      AliVParticle* maxTrack = 0;
      for (Int_t i = 0; i < fEmbeddedTrackNIds; i++) {
	AliVParticle *track2 = static_cast<AliVParticle*>(fEmbTracks->At(fEmbeddedTrackIds[i]));
	if (!maxTrack || track2->Pt() > maxTrack->Pt())
	  maxTrack = track2;
      }
      if (maxTrack)
	fHistEmbNotFoundPhiEta[fCentBin]->Fill(maxTrack->Eta(), maxTrack->Phi());
    }

    if (fEmbCaloClusters) {
      DoEmbClusterLoop();
      TLorentzVector maxCluster;
      for (Int_t i = 0; i < fEmbeddedClusterNIds; i++) {
	AliVCluster *cluster2 = static_cast<AliVCluster*>(fEmbCaloClusters->At(fEmbeddedClusterIds[i]));
	TLorentzVector nPart;
	cluster2->GetMomentum(nPart, fVertex);
	if (nPart.Pt() > maxCluster.Pt())
	  maxCluster = nPart;
      }
      if (maxCluster.Pt() > 0)
	fHistEmbNotFoundPhiEta[fCentBin]->Fill(maxCluster.Eta(), maxCluster.Phi());
    }
  }

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskDeltaPt::DoEmbTrackLoop()
{
  // Do track loop.

  fEmbeddedTrackNIds = 0;

  if (!fEmbTracks)
    return;

  Int_t ntracks = fEmbTracks->GetEntriesFast();

  for (Int_t i = 0; i < ntracks; i++) {

    AliVParticle* track = static_cast<AliVParticle*>(fEmbTracks->At(i)); // pointer to reconstructed to track  

    if (!track) {
      AliError(Form("Could not retrieve track %d",i)); 
      continue; 
    }

    AliVTrack* vtrack = dynamic_cast<AliVTrack*>(track); 
    
    if (vtrack && !AcceptTrack(vtrack)) 
      continue;

    if (track->GetLabel() >= 0) {
      fEmbeddedTrackIds[fEmbeddedTrackNIds] = i;
      fEmbeddedTrackNIds++;
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskDeltaPt::DoEmbClusterLoop()
{
  // Do cluster loop.

  fEmbeddedClusterNIds = 0;

  if (!fEmbCaloClusters)
    return;

  Int_t nclusters =  fEmbCaloClusters->GetEntriesFast();

  for (Int_t iClusters = 0; iClusters < nclusters; iClusters++) {
    AliVCluster* cluster = static_cast<AliVCluster*>(fEmbCaloClusters->At(iClusters));
    if (!cluster) {
      AliError(Form("Could not receive cluster %d", iClusters));
      continue;
    }  

    if (!AcceptCluster(cluster)) 
      continue;

    if (cluster->GetLabel() >= 0) {
      fEmbeddedClusterIds[fEmbeddedClusterNIds] = iClusters;
      fEmbeddedClusterNIds++;
    }
  }
}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskDeltaPt::NextEmbeddedJet(Int_t i)
{
  // Do the embedded jet loop.

  static Int_t iJet = 0;

  if (i >= 0)
    iJet = i;

  if (!fEmbJets)
    return 0;

  TLorentzVector maxClusVect;

  const Int_t nembjets = fEmbJets->GetEntriesFast();

  for (; iJet < nembjets; iJet++) {
      
    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fEmbJets->At(iJet));
      
    if (!jet) {
      AliError(Form("Could not receive jet %d", iJet));
      continue;
    } 
      
    if (!AcceptJet(jet))
      continue;

    if (jet->MCPt() >= fMCJetPtThreshold)
      return jet;
  }

  return 0;
}

//________________________________________________________________________
void AliAnalysisTaskDeltaPt::GetRandomCone(Float_t &pt, Float_t &eta, Float_t &phi,
					AliEmcalJet *jet, TClonesArray* tracks, TClonesArray* clusters) const
{
  // Get rigid cone.

  if (!tracks)
    tracks = fTracks;

  if (!clusters)
    clusters = fCaloClusters;

  eta = 0;
  phi = 0;
  pt = 0;

  if (!tracks && !clusters)
    return;

  Float_t LJeta = 999;
  Float_t LJphi = 999;

  if (jet) {
    LJeta = jet->Eta();
    LJphi = jet->Phi();
  }

  Float_t maxEta = fJetMaxEta;
  Float_t minEta = fJetMinEta;
  Float_t maxPhi = fJetMaxPhi;
  Float_t minPhi = fJetMinPhi;

  if (maxPhi > TMath::Pi() * 2) maxPhi = TMath::Pi() * 2;
  if (minPhi < 0) minPhi = 0;
  
  Float_t dLJ = 0;
  Int_t repeats = 0;
  do {
    eta = gRandom->Rndm() * (maxEta - minEta) + minEta;
    phi = gRandom->Rndm() * (maxPhi - minPhi) + minPhi;
    dLJ = TMath::Sqrt((LJeta - eta) * (LJeta - eta) + (LJphi - phi) * (LJphi - phi));
    repeats++;
  } while (dLJ < fMinRC2LJ && repeats < 999);

  if (repeats == 999) {
    AliWarning(Form("%s: Could not get random cone!", GetName()));
    return;
  }

  if (clusters) {
    Int_t nclusters =  clusters->GetEntriesFast();
    for (Int_t iClusters = 0; iClusters < nclusters; iClusters++) {
      AliVCluster* cluster = static_cast<AliVCluster*>(clusters->At(iClusters));
      if (!cluster) {
	AliError(Form("Could not receive cluster %d", iClusters));
	continue;
      }  
      
      if (!AcceptCluster(cluster))
	continue;
      
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, const_cast<Double_t*>(fVertex));
     
      Float_t d = TMath::Sqrt((nPart.Eta() - eta) * (nPart.Eta() - eta) + (nPart.Phi() - phi) * (nPart.Phi() - phi));

      if (d <= fJetRadius) 
	pt += nPart.Pt();
    }
  }

  if (tracks) {
    Int_t ntracks = tracks->GetEntriesFast();
    for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
      AliVTrack* track = static_cast<AliVTrack*>(tracks->At(iTracks));         
      if(!track) {
	AliError(Form("Could not retrieve track %d",iTracks)); 
	continue; 
      }
      
      if (!AcceptTrack(track)) 
	continue;
      
      Float_t tracketa = track->Eta();
      Float_t trackphi = track->Phi();
      
      if (TMath::Abs(trackphi - phi) > TMath::Abs(trackphi - phi + 2 * TMath::Pi()))
	trackphi += 2 * TMath::Pi();
      if (TMath::Abs(trackphi - phi) > TMath::Abs(trackphi - phi - 2 * TMath::Pi()))
	trackphi -= 2 * TMath::Pi();
      
      Float_t d = TMath::Sqrt((tracketa - eta) * (tracketa - eta) + (trackphi - phi) * (trackphi - phi));
      if (d <= fJetRadius)
	pt += track->Pt();
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskDeltaPt::ExecOnce()
{
  // Initialize the analysis.

  if (!fEmbJetsName.IsNull() && !fEmbJets) {
    fEmbJets =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fEmbJetsName));
    if (!fEmbJets) {
      AliError(Form("%s: Could not retrieve embedded jets %s!", GetName(), fEmbJetsName.Data()));
      return;
    }
    else if (!fEmbJets->GetClass()->GetBaseClass("AliEmcalJet")) {
      AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fEmbJetsName.Data())); 
      fEmbJets = 0;
      return;
    }
  }

  if (!fEmbCaloName.IsNull() && !fEmbCaloClusters) {
    fEmbCaloClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fEmbCaloName));
    if (!fEmbCaloClusters) {
      AliError(Form("%s: Could not retrieve embedded clusters %s!", GetName(), fEmbCaloName.Data()));
      return;
    }
    else if (!fEmbCaloClusters->GetClass()->GetBaseClass("AliVCluster") && !fEmbCaloClusters->GetClass()->GetBaseClass("AliEmcalParticle")) {
      AliError(Form("%s: Collection %s does not contain AliVCluster nor AliEmcalParticle objects!", GetName(), fEmbCaloName.Data())); 
      fEmbCaloClusters = 0;
      return;
    }
  }

  if (!fEmbTracksName.IsNull() && !fEmbTracks) {
    fEmbTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fEmbTracksName));
    if (!fEmbTracks) {
      AliError(Form("%s: Could not retrieve embedded tracks %s!", GetName(), fEmbTracksName.Data()));
      return;
    }
    else if (!fEmbTracks->GetClass()->GetBaseClass("AliVParticle") && !fEmbTracks->GetClass()->GetBaseClass("AliEmcalParticle")) {
      AliError(Form("%s: Collection %s does not contain AliVParticle nor AliEmcalParticle objects!", GetName(), fEmbTracksName.Data())); 
      fEmbTracks = 0;
      return;
    }
  }

  if (!fRandCaloName.IsNull() && !fRandCaloClusters) {
    fRandCaloClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fRandCaloName));
    if (!fRandCaloClusters) {
      AliError(Form("%s: Could not retrieve randomized clusters %s!", GetName(), fRandCaloName.Data()));
      return;
    }
    else if (!fRandCaloClusters->GetClass()->GetBaseClass("AliVCluster") && !fRandCaloClusters->GetClass()->GetBaseClass("AliEmcalParticle")) {
      AliError(Form("%s: Collection %s does not contain AliVCluster nor AliEmcalParticle objects!", GetName(), fRandCaloName.Data())); 
      fRandCaloClusters = 0;
      return;
    }
  }

  if (!fRandTracksName.IsNull() && !fRandTracks) {
    fRandTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fRandTracksName));
    if (!fRandTracks) {
      AliError(Form("%s: Could not retrieve randomized tracks %s!", GetName(), fRandTracksName.Data()));
      return;
    }
    else if (!fRandTracks->GetClass()->GetBaseClass("AliVParticle") && !fRandTracks->GetClass()->GetBaseClass("AliEmcalParticle")) {
      AliError(Form("%s: Collection %s does not contain AliVParticle nor AliEmcalParticle objects!", GetName(), fRandTracksName.Data())); 
      fRandTracks = 0;
      return;
    }
  }

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fRCperEvent < 0) {
    Double_t area = (fJetMaxEta - fJetMinEta) * (fJetMaxPhi - fJetMinPhi);
    Double_t jetArea = TMath::Pi() * fJetRadius * fJetRadius;
    fRCperEvent = TMath::FloorNint(area / jetArea - 0.5);
    if (fRCperEvent == 0)
      fRCperEvent = 1;
  }

  if (fMinRC2LJ < 0)
    fMinRC2LJ = fJetRadius * 1.5;

  const Float_t maxDist = TMath::Max(fJetMaxPhi - fJetMinPhi, fJetMaxEta - fJetMinEta) / 2;
  if (fMinRC2LJ > maxDist) {
    AliWarning(Form("The parameter fMinRC2LJ = %f is too large for the considered acceptance. "
                    "Will use fMinRC2LJ = %f", fMinRC2LJ, maxDist));
    fMinRC2LJ = maxDist;
  }
}

//________________________________________________________________________
void AliAnalysisTaskDeltaPt::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
