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
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliAnalysisTaskDeltaPt.h"

ClassImp(AliAnalysisTaskDeltaPt)

//________________________________________________________________________
AliAnalysisTaskDeltaPt::AliAnalysisTaskDeltaPt() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskDeltaPt", kTRUE),
  fMCJetPtThreshold(1),
  fMinRC2LJ(-1),
  fRCperEvent(-1),
  fConeRadius(0.2),
  fConeMinEta(-0.9),
  fConeMaxEta(0.9),
  fConeMinPhi(0),
  fConeMaxPhi(TMath::Pi()*2),
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0),
  fEmbJetsCont(0),
  fEmbTracksCont(0),
  fEmbCaloClustersCont(0),
  fRandTracksCont(0),
  fRandCaloClustersCont(0),
  fHistRhovsCent(0),
  fHistRCPhiEta(0), 
  fHistRCPt(0),
  fHistRCPtExLJ(0),
  fHistRCPtExPartialLJ(0), 
  fHistRCPtRand(0),
  fHistRhoVSRCPt(0),
  fHistDeltaPtRCvsEP(0),
  fHistDeltaPtRCExLJ(0),
  fHistDeltaPtRCExPartialLJ(0),
  fHistDeltaPtRCRand(0),
  fHistEmbJetsPtArea(0),
  fHistEmbJetsCorrPtArea(0),
  fHistEmbPartPtvsJetPt(0),
  fHistEmbPartPtvsJetCorrPt(0),
  fHistJetPtvsJetCorrPt(0),
  fHistDistLeadPart2JetAxis(0),
  fHistEmbBkgArea(0),
  fHistRhoVSEmbBkg(0),
  fHistDeltaPtEmbArea(0),
  fHistDeltaPtEmbvsEP(0),
  fHistRCPtExLJVSDPhiLJ(0),
  fHistRCPtExPartialLJVSDPhiLJ(0),
  fHistEmbJetsPhiEta(0),
  fHistLeadPartPhiEta(0)
{
  // Default constructor.

  fHistRCPt = 0;
  fHistRCPtExLJ = 0;
  fHistRCPtExPartialLJ = 0;
  fHistRCPtRand = 0;
  fHistRhoVSRCPt = 0;
  fHistDeltaPtRCvsEP = 0;
  fHistDeltaPtRCExLJ = 0;
  fHistDeltaPtRCExPartialLJ = 0;
  fHistDeltaPtRCRand = 0;
  fHistEmbJetsPtArea = 0;
  fHistEmbJetsCorrPtArea = 0;
  fHistEmbPartPtvsJetPt = 0;
  fHistEmbPartPtvsJetCorrPt = 0;
  fHistJetPtvsJetCorrPt = 0;
  fHistDistLeadPart2JetAxis = 0;
  fHistEmbBkgArea = 0;
  fHistRhoVSEmbBkg = 0;
  fHistDeltaPtEmbArea = 0;
  fHistDeltaPtEmbvsEP = 0;

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskDeltaPt::AliAnalysisTaskDeltaPt(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fMCJetPtThreshold(1),
  fMinRC2LJ(-1),
  fRCperEvent(-1),
  fConeRadius(0.2),
  fConeMinEta(-0.9),
  fConeMaxEta(0.9),
  fConeMinPhi(0),
  fConeMaxPhi(TMath::Pi()*2),
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0),
  fEmbJetsCont(0),
  fEmbTracksCont(0),
  fEmbCaloClustersCont(0),
  fRandTracksCont(0),
  fRandCaloClustersCont(0),
  fHistRhovsCent(0),
  fHistRCPhiEta(0), 
  fHistRCPt(0),
  fHistRCPtExLJ(0),
  fHistRCPtExPartialLJ(0), 
  fHistRCPtRand(0),
  fHistRhoVSRCPt(0),
  fHistDeltaPtRCvsEP(0),
  fHistDeltaPtRCExLJ(0),
  fHistDeltaPtRCExPartialLJ(0),
  fHistDeltaPtRCRand(0),
  fHistEmbJetsPtArea(0),
  fHistEmbJetsCorrPtArea(0),
  fHistEmbPartPtvsJetPt(0),
  fHistEmbPartPtvsJetCorrPt(0),
  fHistJetPtvsJetCorrPt(0),
  fHistDistLeadPart2JetAxis(0),
  fHistEmbBkgArea(0),
  fHistRhoVSEmbBkg(0),
  fHistDeltaPtEmbArea(0),
  fHistDeltaPtEmbvsEP(0),
  fHistRCPtExLJVSDPhiLJ(0),
  fHistRCPtExPartialLJVSDPhiLJ(0),
  fHistEmbJetsPhiEta(0),
  fHistLeadPartPhiEta(0)
{
  // Standard constructor.

  fHistRCPt = 0;
  fHistRCPtExLJ = 0;
  fHistRCPtExPartialLJ = 0;
  fHistRCPtRand = 0;
  fHistRhoVSRCPt = 0;
  fHistDeltaPtRCvsEP = 0;
  fHistDeltaPtRCExLJ = 0;
  fHistDeltaPtRCExPartialLJ = 0;
  fHistDeltaPtRCRand = 0;
  fHistEmbJetsPtArea = 0;
  fHistEmbJetsCorrPtArea = 0;
  fHistEmbPartPtvsJetPt = 0;
  fHistEmbPartPtvsJetCorrPt = 0;
  fHistJetPtvsJetCorrPt = 0;
  fHistDistLeadPart2JetAxis = 0;
  fHistEmbBkgArea = 0;
  fHistRhoVSEmbBkg = 0;
  fHistDeltaPtEmbArea = 0;
  fHistDeltaPtEmbvsEP = 0;

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
void AliAnalysisTaskDeltaPt::AllocateHistogramArrays()
{
  fHistRCPt = new TH1*[fNcentBins];
  fHistRCPtExLJ = new TH1*[fNcentBins];
  fHistRCPtExPartialLJ = new TH1*[fNcentBins];
  fHistRCPtRand = new TH1*[fNcentBins];
  fHistRhoVSRCPt = new TH2*[fNcentBins];
  fHistDeltaPtRCvsEP = new TH2*[fNcentBins];
  fHistDeltaPtRCExLJ = new TH1*[fNcentBins];
  fHistDeltaPtRCExPartialLJ = new TH1*[fNcentBins];
  fHistDeltaPtRCRand = new TH1*[fNcentBins];
  fHistEmbJetsPtArea = new TH3*[fNcentBins];
  fHistEmbJetsCorrPtArea = new TH3*[fNcentBins];
  fHistEmbPartPtvsJetPt = new TH2*[fNcentBins];
  fHistEmbPartPtvsJetCorrPt = new TH2*[fNcentBins];
  fHistJetPtvsJetCorrPt = new TH2*[fNcentBins];
  fHistDistLeadPart2JetAxis = new TH1*[fNcentBins];
  fHistEmbBkgArea = new TH2*[fNcentBins];
  fHistRhoVSEmbBkg = new TH2*[fNcentBins];
  fHistDeltaPtEmbArea = new TH2*[fNcentBins];
  fHistDeltaPtEmbvsEP = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fHistRCPt[i] = 0;
    fHistRCPtExLJ[i] = 0;
    fHistRCPtExPartialLJ[i] = 0;
    fHistRCPtRand[i] = 0;
    fHistRhoVSRCPt[i] = 0;
    fHistDeltaPtRCvsEP[i] = 0;
    fHistDeltaPtRCExLJ[i] = 0;
    fHistDeltaPtRCExPartialLJ[i] = 0;
    fHistDeltaPtRCRand[i] = 0;
    fHistEmbJetsPtArea[i] = 0;
    fHistEmbJetsCorrPtArea[i] = 0;
    fHistEmbPartPtvsJetPt[i] = 0;
    fHistEmbPartPtvsJetCorrPt[i] = 0;
    fHistJetPtvsJetCorrPt[i] = 0;
    fHistDistLeadPart2JetAxis[i] = 0;
    fHistEmbBkgArea[i] = 0;
    fHistRhoVSEmbBkg[i] = 0;
    fHistDeltaPtEmbArea[i] = 0;
    fHistDeltaPtEmbvsEP[i] = 0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskDeltaPt::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  AllocateHistogramArrays();

  fHistRhovsCent = new TH2F("fHistRhovsCent", "fHistRhovsCent", 101, -1,  100, fNbins, 0, fMaxBinPt*2);
  fHistRhovsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistRhovsCent->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
  fOutput->Add(fHistRhovsCent);

  fJetsCont = GetJetContainer("Jets");
  fTracksCont = GetParticleContainer("Tracks");
  fCaloClustersCont = GetClusterContainer("CaloClusters");
  fEmbJetsCont = GetJetContainer("EmbJets");
  fEmbTracksCont = GetParticleContainer("EmbTracks");
  fEmbCaloClustersCont = GetClusterContainer("EmbCaloClusters");
  fRandTracksCont = GetParticleContainer("RandTracks");
  fRandCaloClustersCont = GetClusterContainer("RandCaloClusters");

  if (fTracksCont || fCaloClustersCont) {
    fHistRCPhiEta = new TH2F("fHistRCPhiEta","fHistRCPhiEta", 100, -1, 1, 201, 0, TMath::Pi() * 2.01);
    fHistRCPhiEta->GetXaxis()->SetTitle("#eta");
    fHistRCPhiEta->GetYaxis()->SetTitle("#phi");
    fOutput->Add(fHistRCPhiEta);

    if (fJetsCont) {
      fHistRCPtExLJVSDPhiLJ = new TH2F("fHistRCPtExLJVSDPhiLJ","fHistRCPtExLJVSDPhiLJ", fNbins, fMinBinPt, fMaxBinPt, 128, -1.6, 4.8);
      fHistRCPtExLJVSDPhiLJ->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fHistRCPtExLJVSDPhiLJ->GetYaxis()->SetTitle("#Delta#phi");
      fOutput->Add(fHistRCPtExLJVSDPhiLJ);

      fHistRCPtExPartialLJVSDPhiLJ = new TH2F("fHistRCPtExPartialLJVSDPhiLJ","fHistRCPtExPartialLJVSDPhiLJ", fNbins, fMinBinPt, fMaxBinPt, 128, -1.6, 4.8);
      fHistRCPtExPartialLJVSDPhiLJ->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
      fHistRCPtExPartialLJVSDPhiLJ->GetYaxis()->SetTitle("#Delta#phi");
      fOutput->Add(fHistRCPtExPartialLJVSDPhiLJ);
    }
  }

  if (fEmbJetsCont) {
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
  Double_t binsZ[nbinsZ+1] = {0,1,2,3,4,5,6,7,8,9,10,20,1000};

  Double_t *binsPt       = GenerateFixedBinArray(fNbins, fMinBinPt, fMaxBinPt);
  Double_t *binsCorrPt   = GenerateFixedBinArray(fNbins*2, -fMaxBinPt, fMaxBinPt);
  Double_t *binsArea     = GenerateFixedBinArray(50, 0, 2);

  for (Int_t i = 0; i < fNcentBins; i++) {
    if (fTracksCont || fCaloClustersCont) {
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

      histname = "fHistDeltaPtRCvsEP_";
      histname += i;
      fHistDeltaPtRCvsEP[i] = new TH2F(histname.Data(), histname.Data(), 101, 0, TMath::Pi()*1.01, fNbins * 2, -fMaxBinPt, fMaxBinPt);
      fHistDeltaPtRCvsEP[i]->GetXaxis()->SetTitle("#phi_{RC} - #psi_{RP}");
      fHistDeltaPtRCvsEP[i]->GetYaxis()->SetTitle("#delta#it{p}_{T}^{RC} (GeV/#it{c})");
      fHistDeltaPtRCvsEP[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistDeltaPtRCvsEP[i]);

      if (fJetsCont) {
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

        histname = "fHistRCPtExPartialLJ_";
        histname += i;
        fHistRCPtExPartialLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
        fHistRCPtExPartialLJ[i]->GetXaxis()->SetTitle("#it{p}_{T}^{RC} (GeV/#it{c})");
        fHistRCPtExPartialLJ[i]->GetYaxis()->SetTitle("counts");
        fOutput->Add(fHistRCPtExPartialLJ[i]);

        histname = "fHistDeltaPtRCExPartialLJ_";
        histname += i;
        fHistDeltaPtRCExPartialLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
        fHistDeltaPtRCExPartialLJ[i]->GetXaxis()->SetTitle("#delta#it{p}_{T}^{RC} (GeV/#it{c})");
        fHistDeltaPtRCExPartialLJ[i]->GetYaxis()->SetTitle("counts");
        fOutput->Add(fHistDeltaPtRCExPartialLJ[i]);
      }
    }

    if (fRandTracksCont || fRandCaloClustersCont) {
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

    if (fEmbJetsCont) {
      histname = "fHistEmbJetsPtArea_";
      histname += i;
      fHistEmbJetsPtArea[i] = new TH3F(histname.Data(), histname.Data(), 50, binsArea, fNbins, binsPt, nbinsZ, binsZ);
      fHistEmbJetsPtArea[i]->GetXaxis()->SetTitle("area");
      fHistEmbJetsPtArea[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb,raw} (GeV/#it{c})");
      fOutput->Add(fHistEmbJetsPtArea[i]);

      histname = "fHistEmbJetsCorrPtArea_";
      histname += i;
      fHistEmbJetsCorrPtArea[i] = new TH3F(histname.Data(), histname.Data(), 50, binsArea, fNbins * 2, binsCorrPt, nbinsZ, binsZ);
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
      fHistEmbPartPtvsJetCorrPt[i] = new TH2F(histname.Data(), histname.Data(), 
          fNbins, fMinBinPt, fMaxBinPt, fNbins*2, -fMaxBinPt, fMaxBinPt);
      fHistEmbPartPtvsJetCorrPt[i]->GetXaxis()->SetTitle("#sum#it{p}_{T,const}^{emb} (GeV/#it{c})");
      fHistEmbPartPtvsJetCorrPt[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb} - A#rho (GeV/#it{c})");
      fHistEmbPartPtvsJetCorrPt[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistEmbPartPtvsJetCorrPt[i]);

      histname = "fHistJetPtvsJetCorrPt_";
      histname += i;
      fHistJetPtvsJetCorrPt[i] = new TH2F(histname.Data(), histname.Data(), 
          fNbins, fMinBinPt, fMaxBinPt, fNbins*2, -fMaxBinPt, fMaxBinPt);
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

      histname = "fHistEmbBkgArea_";
      histname += i;
      fHistEmbBkgArea[i] = new TH2F(histname.Data(), histname.Data(), 50, 0, 2, fNbins, fMinBinPt, fMaxBinPt);
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
      fHistDeltaPtEmbArea[i] = new TH2F(histname.Data(), histname.Data(), 
          50, 0, 2, fNbins * 2, -fMaxBinPt, fMaxBinPt);
      fHistDeltaPtEmbArea[i]->GetXaxis()->SetTitle("area");
      fHistDeltaPtEmbArea[i]->GetYaxis()->SetTitle("#delta#it{p}_{T}^{emb} (GeV/#it{c})");
      fHistDeltaPtEmbArea[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistDeltaPtEmbArea[i]);

      histname = "fHistDeltaPtEmbvsEP_";
      histname += i;
      fHistDeltaPtEmbvsEP[i] = new TH2F(histname.Data(), histname.Data(), 101, 0, TMath::Pi()*1.01, fNbins * 2, -fMaxBinPt, fMaxBinPt);
      fHistDeltaPtEmbvsEP[i]->GetXaxis()->SetTitle("#phi_{jet} - #Psi_{EP}");
      fHistDeltaPtEmbvsEP[i]->GetYaxis()->SetTitle("#delta#it{p}_{T}^{emb} (GeV/#it{c})");
      fHistDeltaPtEmbvsEP[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistDeltaPtEmbvsEP[i]);
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

  fHistRhovsCent->Fill(fCent, fJetsCont->GetRhoVal());

  // ************
  // Random cones
  // _________________________________

  const Float_t rcArea = fConeRadius * fConeRadius * TMath::Pi();
  Float_t RCpt = 0;
  Float_t RCeta = 0;
  Float_t RCphi = 0;

  if (fTracksCont || fCaloClustersCont) {

    for (Int_t i = 0; i < fRCperEvent; i++) {
      // Simple random cones
      RCpt = 0;
      RCeta = 0;
      RCphi = 0;
      GetRandomCone(RCpt, RCeta, RCphi, fTracksCont, fCaloClustersCont, 0);
      if (RCpt > 0) {
        fHistRCPhiEta->Fill(RCeta, RCphi);
        fHistRhoVSRCPt[fCentBin]->Fill(fJetsCont->GetRhoVal() * rcArea, RCpt);

        fHistRCPt[fCentBin]->Fill(RCpt);

        Double_t ep = RCphi - fEPV0;
        while (ep < 0) ep += TMath::Pi();
        while (ep >= TMath::Pi()) ep -= TMath::Pi();

        fHistDeltaPtRCvsEP[fCentBin]->Fill(ep, RCpt - rcArea * fJetsCont->GetRhoVal());
      }

      if (fJetsCont) {

        // Random cones far from leading jet
        AliEmcalJet* jet = fJetsCont->GetLeadingJet("rho");

        RCpt = 0;
        RCeta = 0;
        RCphi = 0;
        GetRandomCone(RCpt, RCeta, RCphi, fTracksCont, fCaloClustersCont, jet);
        if (RCpt > 0) {
          if (jet) {
            Float_t dphi = RCphi - jet->Phi();
            if (dphi > 4.8) dphi -= TMath::Pi() * 2;
            if (dphi < -1.6) dphi += TMath::Pi() * 2;
            fHistRCPtExLJVSDPhiLJ->Fill(RCpt, dphi);
          }
          fHistRCPtExLJ[fCentBin]->Fill(RCpt);
          fHistDeltaPtRCExLJ[fCentBin]->Fill(RCpt - rcArea * fJetsCont->GetRhoVal());
        }

        //partial exclusion
        if(fBeamType == kpA) {

          RCpt = 0;
          RCeta = 0;
          RCphi = 0;
          GetRandomCone(RCpt, RCeta, RCphi, fTracksCont, fCaloClustersCont, jet, kTRUE);

          if (RCpt > 0) {
            if (jet) {
              Float_t dphi = RCphi - jet->Phi();
              if (dphi > 4.8) dphi -= TMath::Pi() * 2;
              if (dphi < -1.6) dphi += TMath::Pi() * 2;
              fHistRCPtExPartialLJVSDPhiLJ->Fill(RCpt, dphi);
            }
            fHistRCPtExPartialLJ[fCentBin]->Fill(RCpt);
            fHistDeltaPtRCExPartialLJ[fCentBin]->Fill(RCpt - rcArea * fJetsCont->GetRhoVal());
          }
        }
      }
    }
  }

  // Random cones with randomized particles
  if (fRandTracksCont || fRandCaloClustersCont) {
    RCpt = 0;
    RCeta = 0;
    RCphi = 0;
    GetRandomCone(RCpt, RCeta, RCphi, fRandTracksCont, fRandCaloClustersCont, 0);
    if (RCpt > 0) {
      fHistRCPtRand[fCentBin]->Fill(RCpt);
      fHistDeltaPtRCRand[fCentBin]->Fill(RCpt - rcArea * fJetsCont->GetRhoVal());
    }  
  }

  // ************
  // Embedding
  // _________________________________

  if (fEmbJetsCont) {

    AliEmcalJet *embJet = NextEmbeddedJet(kTRUE);

    while (embJet != 0) {
      TLorentzVector mom;
      fEmbJetsCont->GetLeadingHadronMomentum(mom,embJet);

      Double_t distLeading2Jet = TMath::Sqrt((embJet->Eta() - mom.Eta()) * (embJet->Eta() - mom.Eta()) + (embJet->Phi() - mom.Phi()) * (embJet->Phi() - mom.Phi()));

      fHistEmbPartPtvsJetPt[fCentBin]->Fill(embJet->MCPt(), embJet->Pt());
      fHistEmbPartPtvsJetCorrPt[fCentBin]->Fill(embJet->MCPt(), embJet->Pt() - embJet->Area() * fJetsCont->GetRhoVal());
      fHistLeadPartPhiEta->Fill(mom.Eta(), mom.Phi());
      fHistDistLeadPart2JetAxis[fCentBin]->Fill(distLeading2Jet);

      fHistEmbJetsPtArea[fCentBin]->Fill(embJet->Area(), embJet->Pt(), mom.Pt());
      fHistEmbJetsCorrPtArea[fCentBin]->Fill(embJet->Area(), embJet->Pt() - fJetsCont->GetRhoVal() * embJet->Area(), mom.Pt());
      fHistEmbJetsPhiEta->Fill(embJet->Eta(), embJet->Phi());
      fHistJetPtvsJetCorrPt[fCentBin]->Fill(embJet->Pt(), embJet->Pt() - fJetsCont->GetRhoVal() * embJet->Area());

      fHistEmbBkgArea[fCentBin]->Fill(embJet->Area(), embJet->Pt() - embJet->MCPt());
      fHistRhoVSEmbBkg[fCentBin]->Fill(fJetsCont->GetRhoVal() * embJet->Area(), embJet->Pt() - embJet->MCPt());
      fHistDeltaPtEmbArea[fCentBin]->Fill(embJet->Area(), embJet->Pt() - embJet->Area() * fJetsCont->GetRhoVal() - embJet->MCPt());

      Double_t ep = embJet->Phi() - fEPV0;
      while (ep < 0) ep += TMath::Pi();
      while (ep >= TMath::Pi()) ep -= TMath::Pi();

      fHistDeltaPtEmbvsEP[fCentBin]->Fill(ep, embJet->Pt() - embJet->Area() * fJetsCont->GetRhoVal() - embJet->MCPt());

      embJet = NextEmbeddedJet();
    }
  }

  return kTRUE;
}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskDeltaPt::NextEmbeddedJet(Bool_t reset)
{
  // Get the next accepted embedded jet.

  if(reset)
    fEmbJetsCont->ResetCurrentID();

  AliEmcalJet* jet = fEmbJetsCont->GetNextAcceptJet();
  while (jet && jet->MCPt() < fMCJetPtThreshold) jet = fEmbJetsCont->GetNextAcceptJet();

  return jet;
}

//________________________________________________________________________
void AliAnalysisTaskDeltaPt::GetRandomCone(Float_t &pt, Float_t &eta, Float_t &phi,
    AliParticleContainer* tracks, AliClusterContainer* clusters,
    AliEmcalJet *jet, Bool_t bPartialExclusion) const
{
  // Get rigid cone.

  eta = -999;
  phi = -999;
  pt = 0;

  if (!tracks && !clusters)
    return;

  Float_t LJeta = 999;
  Float_t LJphi = 999;

  if (jet) {
    LJeta = jet->Eta();
    LJphi = jet->Phi();
  }

  Float_t maxEta = fConeMaxEta;
  Float_t minEta = fConeMinEta;
  Float_t maxPhi = fConeMaxPhi;
  Float_t minPhi = fConeMinPhi;

  if (maxPhi > TMath::Pi() * 2) maxPhi = TMath::Pi() * 2;
  if (minPhi < 0) minPhi = 0;

  Float_t dLJ = 0;
  Int_t repeats = 0;
  Bool_t reject = kTRUE;
  do {
    eta = gRandom->Rndm() * (maxEta - minEta) + minEta;
    phi = gRandom->Rndm() * (maxPhi - minPhi) + minPhi;
    dLJ = TMath::Sqrt((LJeta - eta) * (LJeta - eta) + (LJphi - phi) * (LJphi - phi));

    if(bPartialExclusion) {
      reject = kFALSE;

      TRandom3 rnd;
      rnd.SetSeed(0);

      Double_t ncoll = GetNColl();

      Double_t prob = 0.;
      if(ncoll>0)
        prob = 1./ncoll;

      if(rnd.Rndm()<=prob) reject = kTRUE; //reject cone
    }

    repeats++;
  } while (dLJ < fMinRC2LJ && repeats < 999 && reject);

  if (repeats == 999) {
    AliWarning(Form("%s: Could not get random cone!", GetName()));
    return;
  }

  if (clusters) {
    clusters->ResetCurrentID();
    AliVCluster* cluster = clusters->GetNextAcceptCluster();
    while (cluster) {     
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, const_cast<Double_t*>(fVertex));

      Float_t cluseta = nPart.Eta();
      Float_t clusphi = nPart.Phi();

      if (TMath::Abs(clusphi - phi) > TMath::Abs(clusphi - phi + 2 * TMath::Pi()))
        clusphi += 2 * TMath::Pi();
      if (TMath::Abs(clusphi - phi) > TMath::Abs(clusphi - phi - 2 * TMath::Pi()))
        clusphi -= 2 * TMath::Pi();

      Float_t d = TMath::Sqrt((cluseta - eta) * (cluseta - eta) + (clusphi - phi) * (clusphi - phi));
      if (d <= fConeRadius) 
        pt += nPart.Pt();

      cluster = clusters->GetNextAcceptCluster();
    }
  }

  if (tracks) {
    tracks->ResetCurrentID();
    AliVParticle* track = tracks->GetNextAcceptParticle(); 
    while(track) { 
      Float_t tracketa = track->Eta();
      Float_t trackphi = track->Phi();

      if (TMath::Abs(trackphi - phi) > TMath::Abs(trackphi - phi + 2 * TMath::Pi()))
        trackphi += 2 * TMath::Pi();
      if (TMath::Abs(trackphi - phi) > TMath::Abs(trackphi - phi - 2 * TMath::Pi()))
        trackphi -= 2 * TMath::Pi();

      Float_t d = TMath::Sqrt((tracketa - eta) * (tracketa - eta) + (trackphi - phi) * (trackphi - phi));
      if (d <= fConeRadius)
        pt += track->Pt();

      track = tracks->GetNextAcceptParticle(); 
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskDeltaPt::SetConeEtaPhiEMCAL()
{
  // Set default cuts for full cones

  SetConeEtaLimits(-0.7+fConeRadius,0.7-fConeRadius);
  SetConePhiLimits(1.4+fConeRadius,TMath::Pi()-fConeRadius);
}

//________________________________________________________________________
void AliAnalysisTaskDeltaPt::SetConeEtaPhiTPC()
{
  // Set default cuts for charged cones

  SetConeEtaLimits(-0.9+fConeRadius, 0.9-fConeRadius);
  SetConePhiLimits(-10, 10);
}

//________________________________________________________________________
void AliAnalysisTaskDeltaPt::ExecOnce()
{
  // Initialize the analysis.

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;
  if (fEmbTracksCont && fEmbTracksCont->GetArray() == 0) fEmbTracksCont = 0;
  if (fEmbCaloClustersCont && fEmbCaloClustersCont->GetArray() == 0) fEmbCaloClustersCont = 0;
  if (fRandTracksCont && fRandTracksCont->GetArray() == 0) fRandTracksCont = 0;
  if (fRandCaloClustersCont && fRandCaloClustersCont->GetArray() == 0) fRandCaloClustersCont = 0;
  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fEmbJetsCont && fEmbJetsCont->GetArray() == 0) fEmbJetsCont = 0;

  if (fRCperEvent < 0) {
    Double_t area = (fConeMaxEta - fConeMinEta) * (fConeMaxPhi - fConeMinPhi);
    Double_t rcArea = TMath::Pi() * fConeRadius * fConeRadius;
    fRCperEvent = TMath::FloorNint(area / rcArea - 0.5);
    if (fRCperEvent == 0)
      fRCperEvent = 1;
  }

  if (fMinRC2LJ < 0)
    fMinRC2LJ = fConeRadius * 1.5;

  const Float_t maxDist = TMath::Max(fConeMaxPhi - fConeMinPhi, fConeMaxEta - fConeMinEta) / 2;
  if (fMinRC2LJ > maxDist) {
    AliWarning(Form("The parameter fMinRC2LJ = %f is too large for the considered acceptance. "
        "Will use fMinRC2LJ = %f", fMinRC2LJ, maxDist));
    fMinRC2LJ = maxDist;
  }
}

//________________________________________________________________________
Double_t AliAnalysisTaskDeltaPt::GetNColl() const {
  // Get NColl - returns value of corresponding bin
  // only works for pA
  // values taken from V0A slicing https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PACentStudies#Tables_with_centrality_bins_from

  if(fBeamType == kpA) {

    const Int_t nNCollBins = 7;

    Double_t centMin[nNCollBins] = {0.,5.,10.,20.,40.,60.,80.};
    Double_t centMax[nNCollBins] = {5.,10.,20.,40.,60.,80.,100.};

    Double_t nColl[nNCollBins] = {14.7,13.,11.7,9.38,6.49,3.96,1.52};

    for(Int_t i = 0; i<nNCollBins; i++) {
      if(fCent>=centMin[i] && fCent<centMax[i])
        return nColl[i];
    }

    return -1.;
  }
  else {
    AliWarning(Form("%s: Only works for pA analysis. Returning -1",GetName()));
    return -1.;
  }
}
