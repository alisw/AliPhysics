// $Id$
//
// Jet deltaPt task.
//
// Author: S.Aiola

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
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
  fHistEmbJetsPhiEta(0),
  fHistLeadPartPhiEta(0)
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
    fHistEmbRejectedJetsPhiEta[i] = 0;
    fHistEmbRejectedJetsPtArea[i] = 0;
    fHistEmbNotFoundPt[i] = 0;
    fHistEmbNotFoundPhiEta[i] = 0;
    fHistEmbJetsPtArea[i] = 0;
    fHistEmbJetsCorrPtArea[i] = 0;
    fHistEmbPartPtvsJetPt[i] = 0;
    fHistEmbPartPtvsJetCorrPt[i] = 0;
    fHistJetPtvsJetCorrPt[i] = 0;
    fHistDistLeadPart2JetAxis[i] = 0;
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
  fHistEmbJetsPhiEta(0),
  fHistLeadPartPhiEta(0)
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
    fHistEmbRejectedJetsPhiEta[i] = 0;
    fHistEmbRejectedJetsPtArea[i] = 0;
    fHistEmbNotFoundPt[i] = 0;
    fHistEmbNotFoundPhiEta[i] = 0;
    fHistEmbJetsPtArea[i] = 0;
    fHistEmbJetsCorrPtArea[i] = 0;
    fHistEmbPartPtvsJetPt[i] = 0;
    fHistEmbPartPtvsJetCorrPt[i] = 0;
    fHistJetPtvsJetCorrPt[i] = 0;
    fHistDistLeadPart2JetAxis[i] = 0;
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

  if (!fTracksName.IsNull() || !fCaloName.IsNull()) {
    fHistRCPhiEta = new TH2F("fHistRCPhiEta","fHistRCPhiEta", 100, -1, 1, 201, 0, TMath::Pi() * 2.01);
    fHistRCPhiEta->GetXaxis()->SetTitle("#eta");
    fHistRCPhiEta->GetYaxis()->SetTitle("#phi");
    fOutput->Add(fHistRCPhiEta);

    if (!fJetsName.IsNull()) {
      fHistRCPtExLJVSDPhiLJ = new TH2F("fHistRCPtExLJVSDPhiLJ","fHistRCPtExLJVSDPhiLJ", fNbins, fMinBinPt, fMaxBinPt, 128, -1.6, 4.8);
      fHistRCPtExLJVSDPhiLJ->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fHistRCPtExLJVSDPhiLJ->GetYaxis()->SetTitle("#Delta#phi");
      fOutput->Add(fHistRCPtExLJVSDPhiLJ);
    }
  }

  if (!fEmbJetsName.IsNull()) {
    fHistEmbJetsPhiEta = new TH2F("fHistEmbJetsPhiEta","fHistEmbJetsPhiEta", 100, -1, 1, 201, 0, TMath::Pi() * 2.01);
    fHistEmbJetsPhiEta->GetXaxis()->SetTitle("#eta");
    fHistEmbJetsPhiEta->GetYaxis()->SetTitle("#phi");
    fOutput->Add(fHistEmbJetsPhiEta);
    
    fHistLeadPartPhiEta = new TH2F("fHistLeadPartPhiEta","fHistLeadPartPhiEta", 100, -1, 1, 201, 0, TMath::Pi() * 2.01);
    fHistLeadPartPhiEta->GetXaxis()->SetTitle("#eta");
    fHistLeadPartPhiEta->GetYaxis()->SetTitle("#phi");
    fOutput->Add(fHistLeadPartPhiEta);
  }

  TString histname;

  const Int_t nbinsZ = 12;
  Float_t binsZ[nbinsZ+1] = {0,1,2,3,4,5,6,7,8,9,10,20,1000};

  Float_t *binsPt       = GenerateFixedBinArray(fNbins, fMinBinPt, fMaxBinPt);
  Float_t *binsCorrPt   = GenerateFixedBinArray(fNbins*2, -fMaxBinPt, fMaxBinPt);
  Float_t *binsArea     = GenerateFixedBinArray(40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3);

  for (Int_t i = 0; i < fNcentBins; i++) {
    if (!fTracksName.IsNull() || !fCaloName.IsNull()) {
      histname = "fHistRCPt_";
      histname += i;
      fHistRCPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
      fHistRCPt[i]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fHistRCPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistRCPt[i]);

      histname = "fHistRhoVSRCPt_";
      histname += i;
      fHistRhoVSRCPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
      fHistRhoVSRCPt[i]->GetXaxis()->SetTitle("A#rho (GeV/#it{c})");
      fHistRhoVSRCPt[i]->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fOutput->Add(fHistRhoVSRCPt[i]);

      histname = "fHistDeltaPtRC_";
      histname += i;
      fHistDeltaPtRC[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
      fHistDeltaPtRC[i]->GetXaxis()->SetTitle("#delta#it{p}_{T}^{RC} (GeV/#it{c})");
      fHistDeltaPtRC[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistDeltaPtRC[i]);
      
      if (!fJetsName.IsNull()) {
	histname = "fHistRCPtExLJ_";
	histname += i;
	fHistRCPtExLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
	fHistRCPtExLJ[i]->GetXaxis()->SetTitle("#it{p}_{T}^{RC} (GeV/#it{c})");
	fHistRCPtExLJ[i]->GetYaxis()->SetTitle("counts");
	fOutput->Add(fHistRCPtExLJ[i]);

	histname = "fHistDeltaPtRCExLJ_";
	histname += i;
	fHistDeltaPtRCExLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
	fHistDeltaPtRCExLJ[i]->GetXaxis()->SetTitle("#delta#it{p}_{T}^{RC} (GeV/#it{c})");
	fHistDeltaPtRCExLJ[i]->GetYaxis()->SetTitle("counts");
	fOutput->Add(fHistDeltaPtRCExLJ[i]);
      }
    }

    if (!fRandTracksName.IsNull() || !fRandCaloName.IsNull()) {
      histname = "fHistRCPtRand_";
      histname += i;
      fHistRCPtRand[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
      fHistRCPtRand[i]->GetXaxis()->SetTitle("#it{p}_{T}^{RC} (GeV/#it{c})");
      fHistRCPtRand[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistRCPtRand[i]);

      histname = "fHistDeltaPtRCRand_";
      histname += i;
      fHistDeltaPtRCRand[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
      fHistDeltaPtRCRand[i]->GetXaxis()->SetTitle("#delta#it{p}_{T}^{RC} (GeV/#it{c})");
      fHistDeltaPtRCRand[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistDeltaPtRCRand[i]);
    }

    if (!fEmbJetsName.IsNull()) {
      histname = "fHistEmbJetsPtArea_";
      histname += i;
      fHistEmbJetsPtArea[i] = new TH3F(histname.Data(), histname.Data(), 40, binsArea, fNbins, binsPt, nbinsZ, binsZ);
      fHistEmbJetsPtArea[i]->GetXaxis()->SetTitle("area");
      fHistEmbJetsPtArea[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb,raw} (GeV/#it{c})");
      fOutput->Add(fHistEmbJetsPtArea[i]);

      histname = "fHistEmbJetsCorrPtArea_";
      histname += i;
      fHistEmbJetsCorrPtArea[i] = new TH3F(histname.Data(), histname.Data(), 40, binsArea, fNbins * 2, binsCorrPt, nbinsZ, binsZ);
      fHistEmbJetsCorrPtArea[i]->GetXaxis()->SetTitle("area");
      fHistEmbJetsCorrPtArea[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb,corr} (GeV/#it{c})");
      fOutput->Add(fHistEmbJetsCorrPtArea[i]);

      histname = "fHistEmbPartPtvsJetPt_";
      histname += i;
      fHistEmbPartPtvsJetPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
      fHistEmbPartPtvsJetPt[i]->GetXaxis()->SetTitle("#sum#it{p}_{T,const}^{emb} (GeV/#it{c})");
      fHistEmbPartPtvsJetPt[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb} (GeV/#it{c})");
      fHistEmbPartPtvsJetPt[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistEmbPartPtvsJetPt[i]);

      histname = "fHistEmbPartPtvsJetCorrPt_";
      histname += i;
      fHistEmbPartPtvsJetCorrPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins*2, -fMaxBinPt, fMaxBinPt);
      fHistEmbPartPtvsJetCorrPt[i]->GetXaxis()->SetTitle("#sum#it{p}_{T,const}^{emb} (GeV/#it{c})");
      fHistEmbPartPtvsJetCorrPt[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb} - A#rho (GeV/#it{c})");
      fHistEmbPartPtvsJetCorrPt[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistEmbPartPtvsJetCorrPt[i]);

      histname = "fHistJetPtvsJetCorrPt_";
      histname += i;
      fHistJetPtvsJetCorrPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins*2, -fMaxBinPt, fMaxBinPt);
      fHistJetPtvsJetCorrPt[i]->GetXaxis()->SetTitle("#it{p}_{T,jet}^{emb} (GeV/#it{c})");
      fHistJetPtvsJetCorrPt[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb} - A#rho (GeV/#it{c})");
      fHistJetPtvsJetCorrPt[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistJetPtvsJetCorrPt[i]);

      histname = "fHistDistLeadPart2JetAxis_";
      histname += i;
      fHistDistLeadPart2JetAxis[i] = new TH1F(histname.Data(), histname.Data(), 50, 0, 0.5);
      fHistDistLeadPart2JetAxis[i]->GetXaxis()->SetTitle("distance");
      fHistDistLeadPart2JetAxis[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistDistLeadPart2JetAxis[i]);

      histname = "fHistEmbNotFoundPhiEta_";
      histname += i;
      fHistEmbNotFoundPhiEta[i] = new TH2F(histname.Data(), histname.Data(), 100, -1, 1, 201, 0, TMath::Pi() * 2.01);
      fHistEmbNotFoundPhiEta[i]->GetXaxis()->SetTitle("#eta");
      fHistEmbNotFoundPhiEta[i]->GetYaxis()->SetTitle("#phi");
      fOutput->Add(fHistEmbNotFoundPhiEta[i]);

      histname = "fHistEmbNotFoundPt_";
      histname += i;
      fHistEmbNotFoundPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
      fHistEmbNotFoundPt[i]->GetXaxis()->SetTitle("#it{p}_{T,const}^{emb} (GeV/#it{c})");
      fHistEmbNotFoundPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistEmbNotFoundPt[i]);

      histname = "fHistEmbRejectedJetsPhiEta_";
      histname += i;
      fHistEmbRejectedJetsPhiEta[i] = new TH2F(histname.Data(), histname.Data(), 100, -1, 1, 201, 0, TMath::Pi() * 2.01);
      fHistEmbRejectedJetsPhiEta[i]->GetXaxis()->SetTitle("#eta");
      fHistEmbRejectedJetsPhiEta[i]->GetYaxis()->SetTitle("#phi");
      fOutput->Add(fHistEmbRejectedJetsPhiEta[i]);

      histname = "fHistEmbRejectedJetsPtArea_";
      histname += i;
      fHistEmbRejectedJetsPtArea[i] = new TH2F(histname.Data(), histname.Data(), 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
      fHistEmbRejectedJetsPtArea[i]->GetXaxis()->SetTitle("area");
      fHistEmbRejectedJetsPtArea[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb,raw} (GeV/#it{c})");
      fOutput->Add(fHistEmbRejectedJetsPtArea[i]);

      histname = "fHistEmbBkgArea_";
      histname += i;
      fHistEmbBkgArea[i] = new TH2F(histname.Data(), histname.Data(), 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
      fHistEmbBkgArea[i]->GetXaxis()->SetTitle("area");
      fHistEmbBkgArea[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb} - #sum#it{p}_{T,const}^{emb} (GeV/#it{c})");
      fOutput->Add(fHistEmbBkgArea[i]);

      histname = "fHistRhoVSEmbBkg_";
      histname += i;
      fHistRhoVSEmbBkg[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
      fHistRhoVSEmbBkg[i]->GetXaxis()->SetTitle("A#rho (GeV/#it{c})");
      fHistRhoVSEmbBkg[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb} - #sum#it{p}_{T,const}^{emb} (GeV/#it{c})");
      fOutput->Add(fHistRhoVSEmbBkg[i]);
      
      histname = "fHistDeltaPtEmbArea_";
      histname += i;
      fHistDeltaPtEmbArea[i] = new TH2F(histname.Data(), histname.Data(), 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins * 2, -fMaxBinPt, fMaxBinPt);
      fHistDeltaPtEmbArea[i]->GetXaxis()->SetTitle("area");
      fHistDeltaPtEmbArea[i]->GetYaxis()->SetTitle("#delta#it{p}_{T}^{emb} (GeV/#it{c})");
      fOutput->Add(fHistDeltaPtEmbArea[i]);
    }
  }

  delete[] binsPt;
  delete[] binsCorrPt;
  delete[] binsArea;

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
Bool_t AliAnalysisTaskDeltaPt::FillHistograms()
{
  // Fill histograms.

  // ************
  // Random cones
  // _________________________________
  
  const Float_t rcArea = fJetRadius * fJetRadius * TMath::Pi();
  Float_t RCpt = 0;
  Float_t RCeta = 0;
  Float_t RCphi = 0;
  
  if (fTracks || fCaloClusters) {
    
    for (Int_t i = 0; i < fRCperEvent; i++) {
      // Simple random cones
      RCpt = 0;
      RCeta = 0;
      RCphi = 0;
      GetRandomCone(RCpt, RCeta, RCphi, 0);
      if (RCpt > 0) {
	fHistRCPhiEta->Fill(RCeta, RCphi);
	fHistRhoVSRCPt[fCentBin]->Fill(fRhoVal * rcArea, RCpt);
	
	fHistRCPt[fCentBin]->Fill(RCpt);
	fHistDeltaPtRC[fCentBin]->Fill(RCpt - rcArea * fRhoVal);
      }
      
      if (fJets) {

	// Random cones far from leading jet
	static Int_t sortedJets[9999] = {-1};
	GetSortedArray(sortedJets, fJets);
	
	AliEmcalJet* jet = 0;
	
	if (sortedJets[0] >= 0) 
	  jet = static_cast<AliEmcalJet*>(fJets->At(sortedJets[0]));
	
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
      }
    }
  }
  
  // Random cones with randomized particles
  if (fRandTracks || fRandCaloClusters) {
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

  if (fEmbJets) {
    
    AliEmcalJet *embJet = NextEmbeddedJet(0);
    
    Int_t countEmbJets = 0;
    
    while (embJet != 0) {
      AliDebug(2,Form("Elaborating embedded jet n. %d", countEmbJets));
      countEmbJets++;

      if (!AcceptJet(embJet)) {
	AliDebug(2,"Embedded jet not accepted, skipping...");
	fHistEmbRejectedJetsPhiEta[fCentBin]->Fill(embJet->Eta(), embJet->Phi());
	fHistEmbRejectedJetsPtArea[fCentBin]->Fill(embJet->Area(), embJet->Pt());
	
	embJet = NextEmbeddedJet();
	continue;
      }
      
      Double_t maxClusterPt = 0;
      Double_t maxClusterEta = 0;
      Double_t maxClusterPhi = 0;

      Double_t maxTrackPt = 0;
      Double_t maxTrackEta = 0;
      Double_t maxTrackPhi = 0;
      
      Double_t maxPartPt = 0;
      Double_t maxPartEta = 0;
      Double_t maxPartPhi = 0;
      
      if (fLeadingHadronType == 1 || fLeadingHadronType == 2) {
	AliVCluster *cluster = embJet->GetLeadingCluster(fEmbCaloClusters);
	if (cluster) {
	  TLorentzVector nPart;
	  cluster->GetMomentum(nPart, fVertex);
	  
	  maxClusterEta = nPart.Eta();
	  maxClusterPhi = nPart.Phi();
	  maxClusterPt = nPart.Pt();
	}
      }
      
      if (fLeadingHadronType == 0 || fLeadingHadronType == 2) {
	AliVParticle *track = embJet->GetLeadingTrack(fEmbTracks);
	if (track) {
	  maxTrackEta = track->Eta();
	  maxTrackPhi = track->Phi();
	  maxTrackPt = track->Pt();
	}
      }
      
      if (maxTrackPt > maxClusterPt) {
	maxPartPt = maxTrackPt;
	maxPartEta = maxTrackEta;
	maxPartPhi = maxTrackPhi;
      }
      else {
	maxPartPt = maxClusterPt;
	maxPartEta = maxClusterEta;
	maxPartPhi = maxClusterPhi;
      }
      
      Double_t distLeading2Jet = TMath::Sqrt((embJet->Eta() - maxPartEta) * (embJet->Eta() - maxPartEta) + (embJet->Phi() - maxPartPhi) * (embJet->Phi() - maxPartPhi));
      
      fHistEmbPartPtvsJetPt[fCentBin]->Fill(embJet->MCPt(), embJet->Pt());
      fHistEmbPartPtvsJetCorrPt[fCentBin]->Fill(embJet->MCPt(), embJet->Pt() - embJet->Area() * fRhoVal);
      fHistLeadPartPhiEta->Fill(maxPartEta, maxPartPhi);
      fHistDistLeadPart2JetAxis[fCentBin]->Fill(distLeading2Jet);
      
      fHistEmbJetsPtArea[fCentBin]->Fill(embJet->Area(), embJet->Pt(), maxPartPt);
      fHistEmbJetsCorrPtArea[fCentBin]->Fill(embJet->Area(), embJet->Pt() - fRhoVal * embJet->Area(), maxPartPt);
      fHistEmbJetsPhiEta->Fill(embJet->Eta(), embJet->Phi());
      fHistJetPtvsJetCorrPt[fCentBin]->Fill(embJet->Pt(), embJet->Pt() - fRhoVal * embJet->Area());
      
      fHistEmbBkgArea[fCentBin]->Fill(embJet->Area(), embJet->Pt() - embJet->MCPt());
      fHistRhoVSEmbBkg[fCentBin]->Fill(fRhoVal * embJet->Area(), embJet->Pt() - embJet->MCPt());
      fHistDeltaPtEmbArea[fCentBin]->Fill(embJet->Area(), embJet->Pt() - embJet->Area() * fRhoVal - embJet->MCPt());

      embJet = NextEmbeddedJet();
    }

    if (countEmbJets==0) {
      AliDebug(1,"No embedded jets found!");
      if (fEmbTracks) {
	DoEmbTrackLoop();
	for (Int_t i = 0; i < fEmbeddedTrackNIds; i++) {
	  AliDebug(2,Form("Embedded track %d found!",i));
	  AliVParticle *track2 = static_cast<AliVParticle*>(fEmbTracks->At(fEmbeddedTrackIds[i]));
	  if (!track2) continue;
	  fHistEmbNotFoundPhiEta[fCentBin]->Fill(track2->Eta(), track2->Phi());
	  fHistEmbNotFoundPt[fCentBin]->Fill(track2->Pt());
	}
      }
      
      if (fEmbCaloClusters) {
	DoEmbClusterLoop();
	for (Int_t i = 0; i < fEmbeddedClusterNIds; i++) {
	  AliDebug(2,Form("Embedded cluster %d found!",i));
	  AliVCluster *cluster2 = static_cast<AliVCluster*>(fEmbCaloClusters->At(fEmbeddedClusterIds[i]));
	  TLorentzVector nPart;
	  cluster2->GetMomentum(nPart, fVertex);
	  fHistEmbNotFoundPhiEta[fCentBin]->Fill(nPart.Eta(), nPart.Phi());
	  fHistEmbNotFoundPt[fCentBin]->Fill(nPart.Pt());
	}
      }
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

    if (track->GetLabel() > 0) {
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

    if (cluster->GetLabel() > 0) {
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
  else
    iJet++;

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

    if (jet->MCPt() < fMCJetPtThreshold)
      continue;
     
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
