// $Id$
//
// Jet analysis task.
//
// Author: S.Aiola

#include <TObject.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TParameter.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliVEventHandler.h"
#include "AliRhoParameter.h"
#include "AliLog.h"

#include "AliAnalysisTaskSAJF.h"

ClassImp(AliAnalysisTaskSAJF)

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskSAJF", kTRUE),
  fMinRC2LJ(1.0),
  fEmbJetsName("EmbJets"),
  fEmbTracksName(""),
  fEmbCaloName(""),
  fRandTracksName("TracksRandomized"),
  fRandCaloName("CaloClustersRandomized"),
  fRhoName("Rho"),
  fEmbJets(0),
  fEmbTracks(0),
  fEmbCaloClusters(0),
  fRandTracks(0),
  fRandCaloClusters(0),
  fRho(0),
  fEmbeddedClusterId(-1),
  fEmbeddedTrackId(-1),
  fHistCentrality(0),
  fHistDeltaVectorPt(0),
  fHistRhoVSleadJetPt(0),
  fHistRCPhiEta(0),
  fHistRCPtExLJVSDPhiLJ(0),
  fHistRhoVSRCPt(0),
  fHistEmbNotFoundPhiEta(0),
  fHistEmbJetPhiEta(0),
  fHistEmbPartPhiEta(0),
  fHistRhoVSEmbBkg(0)

{
  // Default constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistEvents[i] = 0;
    fHistTracksPt[i] = 0;
    fHistClustersPt[i] = 0;
    fHistJetPhiEta[i] = 0;
    fHistJetsPt[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHist2LeadingJetPt[i] = 0;
    fHistJetsNEFvsPt[i] = 0;
    fHistJetsZvsPt[i] = 0;
    fHistMaxTrackPtvsJetPt[i] = 0;
    fHistMaxClusPtvsJetPt[i] = 0;
    fHistMaxPartPtvsJetPt[i] = 0;
    fHistMaxTrackPtvsJetCorrPt[i] = 0;
    fHistMaxClusPtvsJetCorrPt[i] = 0;
    fHistMaxPartPtvsJetCorrPt[i] = 0;
    fHistRho[i] = 0;
    fHistCorrJetsPt[i] = 0;
    fHistCorrJetsPtArea[i] = 0;
    fHistCorrLeadingJetPt[i] = 0;
    fHistRCPtRigid[i] = 0;
    fHistRCPt[i] = 0;
    fHistRCPtExLJ[i] = 0;
    fHistRCPtRand[i] = 0;
    fHistDeltaPtRCRigid[i] = 0;
    fHistDeltaPtRC[i] = 0;
    fHistDeltaPtRCExLJ[i] = 0;
    fHistDeltaPtRCRand[i] = 0;
    fHistEmbJetsPt[i] = 0;
    fHistEmbJetsCorrPt[i] = 0;
    fHistEmbPart[i] = 0;
    fHistDeltaPtEmb[i] = 0;
  }
}

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fMinRC2LJ(1.0),
  fEmbJetsName("EmbJets"),
  fEmbTracksName(""),
  fEmbCaloName(""),
  fRandTracksName("TracksRandomized"),
  fRandCaloName("CaloClustersRandomized"),
  fRhoName("Rho"),
  fEmbJets(0),
  fEmbTracks(0),
  fEmbCaloClusters(0),
  fRandTracks(0),
  fRandCaloClusters(0),
  fRho(0),
  fEmbeddedClusterId(-1),
  fEmbeddedTrackId(-1),
  fHistCentrality(0),
  fHistDeltaVectorPt(0),
  fHistRhoVSleadJetPt(0),
  fHistRCPhiEta(0),
  fHistRCPtExLJVSDPhiLJ(0),
  fHistRhoVSRCPt(0),
  fHistEmbNotFoundPhiEta(0),
  fHistEmbJetPhiEta(0),
  fHistEmbPartPhiEta(0),
  fHistRhoVSEmbBkg(0)
{
  // Standard constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistEvents[i] = 0;
    fHistTracksPt[i] = 0;
    fHistClustersPt[i] = 0;
    fHistJetPhiEta[i] = 0;
    fHistJetsPt[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHist2LeadingJetPt[i] = 0;
    fHistJetsNEFvsPt[i] = 0;
    fHistJetsZvsPt[i] = 0;
    fHistMaxTrackPtvsJetPt[i] = 0;
    fHistMaxClusPtvsJetPt[i] = 0;
    fHistMaxPartPtvsJetPt[i] = 0;
    fHistMaxTrackPtvsJetCorrPt[i] = 0;
    fHistMaxClusPtvsJetCorrPt[i] = 0;
    fHistMaxPartPtvsJetCorrPt[i] = 0;
    fHistRho[i] = 0;
    fHistCorrJetsPt[i] = 0;
    fHistCorrJetsPtArea[i] = 0;
    fHistCorrLeadingJetPt[i] = 0;
    fHistRCPtRigid[i] = 0;
    fHistRCPt[i] = 0;
    fHistRCPtExLJ[i] = 0;
    fHistRCPtRand[i] = 0;   
    fHistDeltaPtRCRigid[i] = 0;
    fHistDeltaPtRC[i] = 0;
    fHistDeltaPtRCExLJ[i] = 0;
    fHistDeltaPtRCRand[i] = 0;
    fHistEmbJetsPt[i] = 0;
    fHistEmbJetsCorrPt[i] = 0;
    fHistEmbPart[i] = 0;
    fHistDeltaPtEmb[i] = 0;
  }
}

//________________________________________________________________________
AliAnalysisTaskSAJF::~AliAnalysisTaskSAJF()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  
  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner(); 
  
  fHistCentrality = new TH1F("fHistCentrality","fHistCentrality", fNbins, 0, 100);
  fHistCentrality->GetXaxis()->SetTitle("Centrality (%)");
  fHistCentrality->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistCentrality);

  fHistDeltaVectorPt = new TH1F("fHistDeltaVectorPt", "fHistDeltaVectorPt", fNbins, -50, 50);
  fHistDeltaVectorPt->GetXaxis()->SetTitle("#deltap_{T} [GeV/c]");
  fHistDeltaVectorPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaVectorPt);

  fHistRhoVSleadJetPt = new TH2F("fHistRhoVSleadJetPt","fHistRhoVSleadJetPt", fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
  fHistRhoVSleadJetPt->GetXaxis()->SetTitle("#rho * area [GeV/c]");
  fHistRhoVSleadJetPt->GetYaxis()->SetTitle("Leading jet p_{T} [GeV/c]");
  fOutput->Add(fHistRhoVSleadJetPt);

  fHistRCPhiEta = new TH2F("fHistRCPhiEta","Phi-Eta distribution of rigid cones", 40, -2, 2, 64, 0, 6.4);
  fHistRCPhiEta->GetXaxis()->SetTitle("#eta");
  fHistRCPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistRCPhiEta);

  fHistRCPtExLJVSDPhiLJ = new TH2F("fHistRCPtExLJVSDPhiLJ","fHistRCPtExLJVSDPhiLJ", fNbins, fMinBinPt, fMaxBinPt, 128, -1.6, 4.8);
  fHistRCPtExLJVSDPhiLJ->GetXaxis()->SetTitle("rigid cone p_{T} [GeV/c]");
  fHistRCPtExLJVSDPhiLJ->GetYaxis()->SetTitle("#Delta#phi");
  fOutput->Add(fHistRCPtExLJVSDPhiLJ);

  fHistRhoVSRCPt = new TH2F("fHistRhoVSRCPt","fHistRhoVSRCPt", fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
  fHistRhoVSRCPt->GetXaxis()->SetTitle("#rho * area [GeV/c]");
  fHistRhoVSRCPt->GetYaxis()->SetTitle("rigid cone p_{T} [GeV/c]");
  fOutput->Add(fHistRhoVSRCPt);

  if (!fEmbJetsName.IsNull()) {
    fHistEmbNotFoundPhiEta = new TH2F("fHistEmbNotFoundPhiEta","Phi-Eta distribution of not found embedded particles", 40, -2, 2, 64, 0, 6.4);
    fHistEmbNotFoundPhiEta->GetXaxis()->SetTitle("#eta");
    fHistEmbNotFoundPhiEta->GetYaxis()->SetTitle("#phi");
    fOutput->Add(fHistEmbNotFoundPhiEta);

    fHistEmbJetPhiEta = new TH2F("fHistEmbJetPhiEta","Phi-Eta distribution of embedded jets", 40, -2, 2, 64, 0, 6.4);
    fHistEmbJetPhiEta->GetXaxis()->SetTitle("#eta");
    fHistEmbJetPhiEta->GetYaxis()->SetTitle("#phi");
    fOutput->Add(fHistEmbJetPhiEta);
    
    fHistEmbPartPhiEta = new TH2F("fHistEmbPartPhiEta","Phi-Eta distribution of embedded particles", 40, -2, 2, 64, 0, 6.4);
    fHistEmbPartPhiEta->GetXaxis()->SetTitle("#eta");
    fHistEmbPartPhiEta->GetYaxis()->SetTitle("#phi");
    fOutput->Add(fHistEmbPartPhiEta);
    
    fHistRhoVSEmbBkg = new TH2F("fHistRhoVSEmbBkg","fHistRhoVSEmbBkg", fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
    fHistRhoVSEmbBkg->GetXaxis()->SetTitle("rho * area [GeV/c]");
    fHistRhoVSEmbBkg->GetYaxis()->SetTitle("background of embedded track [GeV/c]");
    fOutput->Add(fHistRhoVSEmbBkg);
  }

  TString histname;

  for (Int_t i = 0; i < 4; i++) {
    histname = "fHistEvents_";
    histname += i;
    fHistEvents[i] = new TH1F(histname,histname, 6, 0, 6);
    fHistEvents[i]->GetXaxis()->SetTitle("Event state");
    fHistEvents[i]->GetYaxis()->SetTitle("counts");
    fHistEvents[i]->GetXaxis()->SetBinLabel(1, "Rho <= 0");
    fHistEvents[i]->GetXaxis()->SetBinLabel(2, "No jets");
    fHistEvents[i]->GetXaxis()->SetBinLabel(3, "Max Jet not found");
    fHistEvents[i]->GetXaxis()->SetBinLabel(4, "Max Jet <= 0");
    fHistEvents[i]->GetXaxis()->SetBinLabel(5, "Accepted");
    fOutput->Add(fHistEvents[i]);

    histname = "fHistTracksPt_";
    histname += i;
    fHistTracksPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins / 2.5, fMinBinPt, fMaxBinPt / 2.5);
    fHistTracksPt[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fHistTracksPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistTracksPt[i]);

    if (fAnaType == kEMCAL) {    
      histname = "fHistClustersPt_";
      histname += i;
      fHistClustersPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins / 2.5, fMinBinPt, fMaxBinPt / 2.5);
      fHistClustersPt[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      fHistClustersPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistClustersPt[i]);
    }

    histname = "fHistJetPhiEta_";
    histname += i;
    fHistJetPhiEta[i] = new TH2F(histname.Data(), histname.Data(), 40, -2, 2, 64, 0, 6.4);
    fHistJetPhiEta[i]->GetXaxis()->SetTitle("#eta");
    fHistJetPhiEta[i]->GetYaxis()->SetTitle("#phi");
    fOutput->Add(fHistJetPhiEta[i]);

    histname = "fHistJetsPt_";
    histname += i;
    fHistJetsPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
    fHistJetsPt[i]->GetXaxis()->SetTitle("p_{T}^{raw} [GeV/c]");
    fHistJetsPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistJetsPt[i]);

    histname = "fHistJetsPtArea_";
    histname += i;
    fHistJetsPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 20, 0, fJetRadius * fJetRadius * TMath::Pi() * 1.5);
    fHistJetsPtArea[i]->GetXaxis()->SetTitle("p_{T}^{raw} [GeV/c]");
    fHistJetsPtArea[i]->GetYaxis()->SetTitle("area");
    fOutput->Add(fHistJetsPtArea[i]);

    histname = "fHistLeadingJetPt_";
    histname += i;
    fHistLeadingJetPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
    fHistLeadingJetPt[i]->GetXaxis()->SetTitle("p_{T}^{raw} [GeV]");
    fHistLeadingJetPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistLeadingJetPt[i]);

    histname = "fHist2LeadingJetPt_";
    histname += i;
    fHist2LeadingJetPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
    fHist2LeadingJetPt[i]->GetXaxis()->SetTitle("p_{T}^{raw} [GeV]");
    fHist2LeadingJetPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHist2LeadingJetPt[i]);

    histname = "fHistJetsZvsPt_";
    histname += i;
    fHistJetsZvsPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
    fHistJetsZvsPt[i]->GetXaxis()->SetTitle("Z");
    fHistJetsZvsPt[i]->GetYaxis()->SetTitle("p_{T}^{raw} [GeV/c]");
    fOutput->Add(fHistJetsZvsPt[i]);

    if (fAnaType == kEMCAL) {
      histname = "fHistJetsNEFvsPt_";
      histname += i;
      fHistJetsNEFvsPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
      fHistJetsNEFvsPt[i]->GetXaxis()->SetTitle("NEF");
      fHistJetsNEFvsPt[i]->GetYaxis()->SetTitle("p_{T}^{raw} [GeV/c]");
      fOutput->Add(fHistJetsNEFvsPt[i]);
    }

    histname = "fHistMaxTrackPtvsJetPt_";
    histname += i;
    fHistMaxTrackPtvsJetPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins / 2.5, fMinBinPt, fMaxBinPt / 2.5);
    fHistMaxTrackPtvsJetPt[i]->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
    fHistMaxTrackPtvsJetPt[i]->GetYaxis()->SetTitle("p_{T}^{track} [GeV/c]");
    fOutput->Add(fHistMaxTrackPtvsJetPt[i]);

    if (fAnaType == kEMCAL) {
      histname = "fHistMaxClusPtvsJetPt_";
      histname += i;
      fHistMaxClusPtvsJetPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins / 2.5, fMinBinPt, fMaxBinPt / 2.5);
      fHistMaxClusPtvsJetPt[i]->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
      fHistMaxClusPtvsJetPt[i]->GetYaxis()->SetTitle("p_{T}^{clus} [GeV/c]");
      fOutput->Add(fHistMaxClusPtvsJetPt[i]);
    }

    histname = "fHistMaxPartPtvsJetPt_";
    histname += i;
    fHistMaxPartPtvsJetPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins / 2.5, fMinBinPt, fMaxBinPt / 2.5);
    fHistMaxPartPtvsJetPt[i]->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
    fHistMaxPartPtvsJetPt[i]->GetYaxis()->SetTitle("p_{T}^{part} [GeV/c]");
    fOutput->Add(fHistMaxPartPtvsJetPt[i]);

    histname = "fHistMaxTrackPtvsJetCorrPt_";
    histname += i;
    fHistMaxTrackPtvsJetCorrPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt, fNbins / 2.5, fMinBinPt, fMaxBinPt / 2.5);
    fHistMaxTrackPtvsJetCorrPt[i]->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
    fHistMaxTrackPtvsJetCorrPt[i]->GetYaxis()->SetTitle("p_{T}^{track} [GeV/c]");
    fOutput->Add(fHistMaxTrackPtvsJetCorrPt[i]);

    if (fAnaType == kEMCAL) {
      histname = "fHistMaxClusPtvsJetCorrPt_";
      histname += i;
      fHistMaxClusPtvsJetCorrPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt, fNbins / 2.5, fMinBinPt, fMaxBinPt / 2.5);
      fHistMaxClusPtvsJetCorrPt[i]->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
      fHistMaxClusPtvsJetCorrPt[i]->GetYaxis()->SetTitle("p_{T}^{clus} [GeV/c]");
      fOutput->Add(fHistMaxClusPtvsJetCorrPt[i]);
    }

    histname = "fHistMaxPartPtvsJetCorrPt_";
    histname += i;
    fHistMaxPartPtvsJetCorrPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt, fNbins / 2.5, fMinBinPt, fMaxBinPt / 2.5);
    fHistMaxPartPtvsJetCorrPt[i]->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
    fHistMaxPartPtvsJetCorrPt[i]->GetYaxis()->SetTitle("p_{T}^{part} [GeV/c]");
    fOutput->Add(fHistMaxPartPtvsJetCorrPt[i]);

    histname = "fHistRho_";
    histname += i;
    fHistRho[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
    fHistRho[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fHistRho[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistRho[i]);

    histname = "fHistCorrJetsPt_";
    histname += i;
    fHistCorrJetsPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
    fHistCorrJetsPt[i]->GetXaxis()->SetTitle("p_{T}^{corr} [GeV/c]");
    fHistCorrJetsPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistCorrJetsPt[i]);

    histname = "fHistCorrJetsPtArea_";
    histname += i;
    fHistCorrJetsPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt, 20, 0, fJetRadius * fJetRadius * TMath::Pi() * 1.5);
    fHistCorrJetsPtArea[i]->GetXaxis()->SetTitle("p_{T}^{corr} [GeV/c]");
    fHistCorrJetsPtArea[i]->GetYaxis()->SetTitle("area");
    fOutput->Add(fHistCorrJetsPtArea[i]);

    histname = "fHistCorrLeadingJetPt_";
    histname += i;
    fHistCorrLeadingJetPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
    fHistCorrLeadingJetPt[i]->GetXaxis()->SetTitle("p_{T}^{RC} [GeV/c]");
    fHistCorrLeadingJetPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistCorrLeadingJetPt[i]);

    histname = "fHistRCPtRigid_";
    histname += i;
    fHistRCPtRigid[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
    fHistRCPtRigid[i]->GetXaxis()->SetTitle("rigid cone p_{T} [GeV/c]");
    fHistRCPtRigid[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistRCPtRigid[i]);

    histname = "fHistRCPt_";
    histname += i;
    fHistRCPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
    fHistRCPt[i]->GetXaxis()->SetTitle("rigid cone p_{T} [GeV/c]");
    fHistRCPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistRCPt[i]);

    histname = "fHistRCPtExLJ_";
    histname += i;
    fHistRCPtExLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
    fHistRCPtExLJ[i]->GetXaxis()->SetTitle("rigid cone p_{T}^{RC} [GeV/c]");
    fHistRCPtExLJ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistRCPtExLJ[i]);

    histname = "fHistRCPtRand_";
    histname += i;
    fHistRCPtRand[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
    fHistRCPtRand[i]->GetXaxis()->SetTitle("rigid cone p_{T}^{RC} [GeV/c]");
    fHistRCPtRand[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistRCPtRand[i]);

    histname = "fHistDeltaPtRCRigid_";
    histname += i;
    fHistDeltaPtRCRigid[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
    fHistDeltaPtRCRigid[i]->GetXaxis()->SetTitle("#deltap_{T}^{RC} [GeV/c]");
    fHistDeltaPtRCRigid[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtRCRigid[i]);

    histname = "fHistDeltaPtRC_";
    histname += i;
    fHistDeltaPtRC[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
    fHistDeltaPtRC[i]->GetXaxis()->SetTitle("#deltap_{T}^{RC} [GeV/c]");
    fHistDeltaPtRC[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtRC[i]);

    histname = "fHistDeltaPtRCExLJ_";
    histname += i;
    fHistDeltaPtRCExLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
    fHistDeltaPtRCExLJ[i]->GetXaxis()->SetTitle("#deltap_{T}^{RC} [GeV/c]");
    fHistDeltaPtRCExLJ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtRCExLJ[i]);

    histname = "fHistDeltaPtRCRand_";
    histname += i;
    fHistDeltaPtRCRand[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
    fHistDeltaPtRCRand[i]->GetXaxis()->SetTitle("#deltap_{T}^{RC} [GeV/c]");
    fHistDeltaPtRCRand[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtRCRand[i]);

    if (!fEmbJetsName.IsNull()) {
      histname = "fHistEmbJetsPt_";
      histname += i;
      fHistEmbJetsPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
      fHistEmbJetsPt[i]->GetXaxis()->SetTitle("embedded jet p_{T}^{raw} [GeV/c]");
      fHistEmbJetsPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistEmbJetsPt[i]);

      histname = "fHistEmbJetsCorrPt_";
      histname += i;
      fHistEmbJetsCorrPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
      fHistEmbJetsCorrPt[i]->GetXaxis()->SetTitle("embedded jet p_{T}^{corr} [GeV/c]");
      fHistEmbJetsCorrPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistEmbJetsCorrPt[i]);
      
      histname = "fHistEmbPart_";
      histname += i;
      fHistEmbPart[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
      fHistEmbPart[i]->GetXaxis()->SetTitle("embedded particle p_{T}^{emb} [GeV/c]");
      fHistEmbPart[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistEmbPart[i]);
      
      histname = "fHistDeltaPtEmb_";
      histname += i;
      fHistDeltaPtEmb[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
      fHistDeltaPtEmb[i]->GetXaxis()->SetTitle("#deltap_{T}^{emb} [GeV/c]");
      fHistDeltaPtEmb[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistDeltaPtEmb[i]);
    }
  }

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::RetrieveEventObjects()
{
  // Retrieve event objects.

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;
  
  if (!fRhoName.IsNull() && !fRho) {
    fRho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRhoName));
    if (!fRho) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRhoName.Data()));
      return kFALSE;
    }
  }
  
  if (!fEmbJetsName.IsNull() && !fEmbJets) {
    fEmbJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fEmbJetsName));
    if (!fEmbJets) {
      AliError(Form("%s: Could not retrieve emb jets %s!", GetName(), fEmbJetsName.Data()));
      return kFALSE;
    }
    else if (!fEmbJets->GetClass()->GetBaseClass("AliEmcalJet")) {
      AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fEmbJetsName.Data())); 
      fEmbJets = 0;
      return kFALSE;
    }
  }

  if (!fEmbCaloName.IsNull() && fAnaType == kEMCAL && !fEmbCaloClusters) {
    fEmbCaloClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fEmbCaloName));
    if (!fEmbCaloClusters) {
      AliError(Form("%s: Could not retrieve embedded clusters %s!", GetName(), fEmbCaloName.Data()));
      return kFALSE;
    }
    else if (!fEmbCaloClusters->GetClass()->GetBaseClass("AliVCluster") && !fEmbCaloClusters->GetClass()->GetBaseClass("AliEmcalParticle")) {
      AliError(Form("%s: Collection %s does not contain AliVCluster nor AliEmcalParticle objects!", GetName(), fEmbCaloName.Data())); 
      fEmbCaloClusters = 0;
      return kFALSE;
    }
  }

  if (!fEmbTracksName.IsNull() && !fEmbTracks) {
    fEmbTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fEmbTracksName));
    if (!fEmbTracks) {
      AliError(Form("%s: Could not retrieve embedded tracks %s!", GetName(), fEmbTracksName.Data()));
      return kFALSE;
    }
    else if (!fEmbTracks->GetClass()->GetBaseClass("AliVParticle") && !fEmbTracks->GetClass()->GetBaseClass("AliEmcalParticle")) {
      AliError(Form("%s: Collection %s does not contain AliVParticle nor AliEmcalParticle objects!", GetName(), fEmbTracksName.Data())); 
      fEmbTracks = 0;
      return kFALSE;
    }
  }

  if (!fRandCaloName.IsNull() && fAnaType == kEMCAL && !fRandCaloClusters) {
    fRandCaloClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fRandCaloName));
    if (!fRandCaloClusters) {
      AliError(Form("%s: Could not retrieve randomized clusters %s!", GetName(), fRandCaloName.Data()));
      return kFALSE;
    }
    else if (!fRandCaloClusters->GetClass()->GetBaseClass("AliVCluster") && !fRandCaloClusters->GetClass()->GetBaseClass("AliEmcalParticle")) {
      AliError(Form("%s: Collection %s does not contain AliVCluster nor AliEmcalParticle objects!", GetName(), fRandCaloName.Data())); 
      fRandCaloClusters = 0;
      return kFALSE;
    }
  }

  if (!fRandTracksName.IsNull() && !fRandTracks) {
    fRandTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fRandTracksName));
    if (!fRandTracks) {
      AliError(Form("%s: Could not retrieve randomized tracks %s!", GetName(), fRandTracksName.Data()));
      return kFALSE;
    }
    else if (!fRandTracks->GetClass()->GetBaseClass("AliVParticle") && !fRandTracks->GetClass()->GetBaseClass("AliEmcalParticle")) {
      AliError(Form("%s: Collection %s does not contain AliVParticle nor AliEmcalParticle objects!", GetName(), fRandTracksName.Data())); 
      fRandTracks = 0;
      return kFALSE;
    }
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::FillHistograms()
{
  // Fill histograms.

  // Before filling any histogram, check if the event is interesting
  Double_t rho = fRho->GetVal();

  if (rho < 0) {  // rho < 0, probably a diffractive event, skipping
    fHistEvents[fCentBin]->Fill("Rho <= 0", 1);
    return kTRUE;
  }

  Int_t maxJetIndex  = -1;
  Int_t max2JetIndex = -1;

  GetLeadingJets(maxJetIndex, max2JetIndex);
  
  if (maxJetIndex < 0) { // no accept jet, skipping
    fHistEvents[fCentBin]->Fill("No jets", 1);
    return kTRUE;
  }

  AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(maxJetIndex));
  if (!jet) {  // error, I cannot get the lead jet from collection (should never happen), skipping
    fHistEvents[fCentBin]->Fill("Max Jet not found", 1);
    return kTRUE;
  }

  // OK, event accepted

  Float_t maxJetCorrPt = jet->Pt() - rho * jet->Area();

  if (maxJetCorrPt <= 0)
    fHistEvents[fCentBin]->Fill("Max Jet <= 0", 1);

  fHistEvents[fCentBin]->Fill("Accepted", 1);

  // ************
  // General histograms
  // _________________________________

  DoJetLoop();
  DoTrackLoop();
  DoClusterLoop();

  fHistCentrality->Fill(fCent);
  fHistRho[fCentBin]->Fill(rho);

  if (jet) {
    fHistLeadingJetPt[fCentBin]->Fill(jet->Pt());
    fHistRhoVSleadJetPt->Fill(rho * jet->Area(), jet->Pt());
    fHistCorrLeadingJetPt[fCentBin]->Fill(maxJetCorrPt);
  }
  
  AliEmcalJet* jet2 = 0;
  if (max2JetIndex >= 0)
    jet2 = static_cast<AliEmcalJet*>(fJets->At(max2JetIndex));

  if (jet2)
    fHist2LeadingJetPt[fCentBin]->Fill(jet2->Pt());

  // ************
  // Random cones
  // _________________________________
  
  const Float_t rcArea = fJetRadius * fJetRadius * TMath::Pi();
  const Bool_t IsMCEvent = (Bool_t)(fTracksName.Contains("Randomized") || fTracksName.Contains("Embedded"));
 
  // Simple random cones
  Float_t RCpt = 0;
  Float_t RCptRigid = 0;
  Float_t RCeta = 0;
  Float_t RCphi = 0;
  GetRigidCone(RCpt, RCptRigid, RCeta, RCphi, IsMCEvent, 0);
  if (RCpt > 0) {
    fHistRCPt[fCentBin]->Fill(RCpt / rcArea);
    fHistDeltaPtRC[fCentBin]->Fill(RCpt - rcArea * rho);
  }
  if (RCptRigid > 0) {
    fHistRCPtRigid[fCentBin]->Fill(RCptRigid / rcArea);
    fHistDeltaPtRCRigid[fCentBin]->Fill(RCptRigid - rcArea * rho);
  }
  
  // Random cones far from leading jet
  RCpt = 0;
  RCptRigid = 0;
  RCeta = 0;
  RCphi = 0;
  GetRigidCone(RCpt, RCptRigid, RCeta, RCphi, IsMCEvent, jet);
  if (RCpt > 0) {
    fHistRCPhiEta->Fill(RCeta, RCphi);
    fHistRhoVSRCPt->Fill(rho, RCpt / rcArea);

    Float_t dphi = RCphi - jet->Phi();
    if (dphi > 4.8) dphi -= TMath::Pi() * 2;
    if (dphi < -1.6) dphi += TMath::Pi() * 2; 
    fHistRCPtExLJVSDPhiLJ->Fill(RCpt, dphi);
    
    fHistRCPtExLJ[fCentBin]->Fill(RCpt / rcArea);
    fHistDeltaPtRCExLJ[fCentBin]->Fill(RCpt - rcArea * rho);
  }

  // Random cones with randomized particles
  RCpt = 0;
  RCptRigid = 0;
  RCeta = 0;
  RCphi = 0;
  GetRigidCone(RCpt, RCptRigid, RCeta, RCphi, kTRUE, 0, fRandTracks, fRandCaloClusters);
  if (RCpt > 0) {
    fHistRCPtRand[fCentBin]->Fill(RCpt / rcArea);
    fHistDeltaPtRCRand[fCentBin]->Fill(RCpt - rcArea * rho);
  }

  // ************
  // Embedding
  // _________________________________

  if (!fEmbJets)
    return kTRUE;

  AliEmcalJet *embJet  = 0;
  TObject     *embPart = 0;

  DoEmbJetLoop(embJet, embPart);

  if (embJet) {
    Double_t probePt = 0;

    AliVCluster *cluster = dynamic_cast<AliVCluster*>(embPart);
    if (cluster) {
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);
      
      fHistEmbPartPhiEta->Fill(nPart.Eta(), nPart.Phi());
      probePt = nPart.Pt();
    }
    else {
      AliVTrack *track = dynamic_cast<AliVTrack*>(embPart);
      if (track) {
	fHistEmbPartPhiEta->Fill(track->Eta(), track->Phi());
	probePt = track->Pt();
      }
      else {
	  AliWarning(Form("%s - Embedded jet found but embedded particle not found (wrong type?)!", GetName()));
	return kTRUE;
      }
    }

    fHistEmbJetsPt[fCentBin]->Fill(embJet->Pt());
    fHistEmbJetsCorrPt[fCentBin]->Fill(embJet->Pt() - rho * embJet->Area());
    fHistEmbPart[fCentBin]->Fill(probePt);
    fHistEmbJetPhiEta->Fill(embJet->Eta(), embJet->Phi());
    fHistDeltaPtEmb[fCentBin]->Fill(embJet->Pt() - embJet->Area() * rho - probePt);
    fHistRhoVSEmbBkg->Fill(embJet->Area() * rho, embJet->Pt() - probePt);

  }
  else {
    if (fEmbeddedTrackId >= 0) {
      AliVTrack *track2 = static_cast<AliVTrack*>(fTracks->At(fEmbeddedTrackId));
      fHistEmbNotFoundPhiEta->Fill(track2->Eta(), track2->Phi());
    }
    else if (fEmbeddedClusterId >= 0) {
      AliVCluster *cluster2 = static_cast<AliVCluster*>(fCaloClusters->At(fEmbeddedClusterId));
      TLorentzVector nPart;
      cluster2->GetMomentum(nPart, fVertex);
      fHistEmbNotFoundPhiEta->Fill(nPart.Eta(), nPart.Phi());
    }
    else {
      AliWarning(Form("%s - Embedded particle not found!", GetName()));
    }
  }

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::GetLeadingJets(Int_t &maxJetIndex, Int_t &max2JetIndex)
{
  // Get the leading jets.

  if (!fJets)
    return;

  const Double_t rho = fRho->GetVal();

  const Int_t njets = fJets->GetEntriesFast();

  Float_t maxJetPt = -999;
  Float_t max2JetPt = -999;
  for (Int_t ij = 0; ij < njets; ij++) {

    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(ij));

    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    }  

    if (!AcceptJet(jet))
      continue;

    Float_t corrPt = jet->Pt() - rho * jet->Area();

    if (maxJetIndex == -1 || maxJetPt < corrPt) {
      max2JetPt = maxJetPt;
      max2JetIndex = maxJetIndex;
      maxJetPt = corrPt;
      maxJetIndex = ij;
    }
    else if (max2JetIndex == -1 || max2JetPt < corrPt) {
      max2JetPt = corrPt;
      max2JetIndex = ij;
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::DoClusterLoop()
{
  // Do cluster loop.

  if (!fCaloClusters)
    return;

  fEmbeddedClusterId = -1;

  Int_t nclusters =  fCaloClusters->GetEntriesFast();

  for (Int_t iClusters = 0; iClusters < nclusters; iClusters++) {
    AliVCluster* cluster = static_cast<AliVCluster*>(fCaloClusters->At(iClusters));
    if (!cluster) {
      AliError(Form("Could not receive cluster %d", iClusters));
      continue;
    }  

    if (!AcceptCluster(cluster, kTRUE)) 
      continue;

    if (cluster->Chi2() == 100)
      fEmbeddedClusterId = iClusters;

    fHistClustersPt[fCentBin]->Fill(cluster->E());
  }
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::DoTrackLoop()
{
  // Do track loop.

  if (!fTracks)
    return;

  fEmbeddedTrackId = -1;

  Int_t ntracks = fTracks->GetEntriesFast();

  for (Int_t i = 0; i < ntracks; i++) {

    AliVParticle* track = static_cast<AliVParticle*>(fTracks->At(i)); // pointer to reconstructed to track  

    if (!track) {
      AliError(Form("Could not retrieve track %d",i)); 
      continue; 
    }

    AliVTrack* vtrack = dynamic_cast<AliVTrack*>(track); 
    
    if (vtrack && !AcceptTrack(vtrack, kTRUE)) 
      continue;

    if (track->GetLabel() == 100)
      fEmbeddedTrackId = i;
    
    fHistTracksPt[fCentBin]->Fill(track->Pt());
  }
}


//________________________________________________________________________
void AliAnalysisTaskSAJF::DoJetLoop()
{
  // Do the jet loop.

  if (!fJets)
    return;

  const Double_t rho = fRho->GetVal();

  const Int_t njets = fJets->GetEntriesFast();

  for (Int_t ij = 0; ij < njets; ij++) {

    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(ij));

    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    }  

    if (!AcceptJet(jet))
      continue;

    Float_t corrPt = jet->Pt() - rho * jet->Area();

    fHistJetsPt[fCentBin]->Fill(jet->Pt());
    fHistJetsPtArea[fCentBin]->Fill(jet->Pt(), jet->Area());
    fHistCorrJetsPt[fCentBin]->Fill(corrPt);
    fHistCorrJetsPtArea[fCentBin]->Fill(corrPt, jet->Area());

    fHistMaxTrackPtvsJetPt[fCentBin]->Fill(jet->Pt(), jet->MaxTrackPt());
    fHistMaxPartPtvsJetPt[fCentBin]->Fill(jet->Pt(), jet->MaxPartPt());

    fHistMaxTrackPtvsJetCorrPt[fCentBin]->Fill(corrPt, jet->MaxTrackPt());
    fHistMaxPartPtvsJetCorrPt[fCentBin]->Fill(corrPt, jet->MaxPartPt());

    fHistJetPhiEta[fCentBin]->Fill(jet->Eta(), jet->Phi());

    if (fAnaType == kEMCAL) {
      fHistMaxClusPtvsJetPt[fCentBin]->Fill(jet->Pt(), jet->MaxClusterPt());
      fHistMaxClusPtvsJetCorrPt[fCentBin]->Fill(corrPt, jet->MaxClusterPt());
      fHistJetsNEFvsPt[fCentBin]->Fill(jet->NEF(), jet->Pt());
    }

    Float_t scalarpt = 0;

    if (fTracks) {
      for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
	AliVParticle *track = jet->TrackAt(it, fTracks);
	if (track) {
	  fHistJetsZvsPt[fCentBin]->Fill(track->Pt() / jet->Pt(), jet->Pt());
	  scalarpt += track->Pt();
	}
      }
    }

    if (fCaloClusters) {
      for (Int_t ic = 0; ic < jet->GetNumberOfClusters(); ic++) {
	AliVCluster *cluster = jet->ClusterAt(ic, fCaloClusters);
	
	if (cluster) {
	  TLorentzVector nPart;
	  cluster->GetMomentum(nPart, fVertex);
	  fHistJetsZvsPt[fCentBin]->Fill(nPart.Et() / jet->Pt(), jet->Pt());
	  scalarpt += nPart.Pt();
	}
      }
    }

    fHistDeltaVectorPt->Fill(scalarpt - jet->Pt());
  } //jet loop 
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::DoEmbJetLoop(AliEmcalJet* &embJet, TObject* &embPart)
{
  // Do the embedded jet loop.

  if (!fEmbJets)
    return;

  TLorentzVector maxClusVect;

  Int_t nembjets = fEmbJets->GetEntriesFast();

  for (Int_t ij = 0; ij < nembjets; ij++) {
      
    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fEmbJets->At(ij));
      
    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    } 
      
    if (!AcceptJet(jet))
      continue;

    if (!jet->IsMC())
      continue;

    AliVParticle *maxTrack = 0;

    for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
      AliVParticle *track = jet->TrackAt(it, fEmbTracks);
      
      if (!track) 
	continue;

      if (track->GetLabel() != 100)
	continue;
      
      if (!maxTrack || track->Pt() > maxTrack->Pt())
	maxTrack = track;
    }
    
    AliVCluster *maxClus = 0;

    for (Int_t ic = 0; ic < jet->GetNumberOfClusters(); ic++) {
      AliVCluster *cluster = jet->ClusterAt(ic, fEmbCaloClusters);
      
      if (!cluster) 
	continue;

      if (cluster->Chi2() != 100)
	continue;
      
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);
      
      if (!maxClus || nPart.Et() > maxClusVect.Et()) {
	new (&maxClusVect) TLorentzVector(nPart);
	maxClus = cluster;
      } 
    }

    if ((maxClus && maxTrack && maxClusVect.Et() > maxTrack->Pt()) || (maxClus && !maxTrack)) {
      embPart = maxClus;
      embJet = jet;
    }
    else if (maxTrack) {
      embPart = maxTrack;
      embJet = jet;
    }

    return;  // MC jets found, exit
  }
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::GetRigidCone(Float_t &pt, Float_t &ptrigid, Float_t &eta, Float_t &phi, Bool_t acceptMC,
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
  ptrigid = 0;

  if (!tracks && !clusters)
    return;

  Float_t LJeta = 999;
  Float_t LJphi = 999;

  if (jet) {
    LJeta = jet->Eta();
    LJphi = jet->Phi();
  }

  Float_t maxEta = fMaxEta;
  Float_t minEta = fMinEta;
  Float_t maxPhi = fMaxPhi;
  Float_t minPhi = fMinPhi;

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
    AliWarning("Could not get random cone!");
    return;
  }

  TVector3 rigidAxis;
  rigidAxis.SetPtEtaPhi(1, eta, phi);
  rigidAxis = rigidAxis.Unit();

  if (fAnaType == kEMCAL && clusters) {
    Int_t nclusters =  clusters->GetEntriesFast();
    for (Int_t iClusters = 0; iClusters < nclusters; iClusters++) {
      AliVCluster* cluster = static_cast<AliVCluster*>(clusters->At(iClusters));
      if (!cluster) {
	AliError(Form("Could not receive cluster %d", iClusters));
	continue;
      }  
      
      if (!AcceptCluster(cluster, acceptMC))
	continue;
      
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, const_cast<Double_t*>(fVertex));
     
      Float_t d = TMath::Sqrt((nPart.Eta() - eta) * (nPart.Eta() - eta) + (nPart.Phi() - phi) * (nPart.Phi() - phi));

      if (d <= fJetRadius) {
	TVector3 vect = nPart.Vect();
	vect *= vect * rigidAxis / vect.Mag();
	ptrigid += vect.Pt();
	
	pt += nPart.Pt();
      }
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
      
      if (!AcceptTrack(track, acceptMC)) 
	continue;
      
      Float_t tracketa = track->Eta();
      Float_t trackphi = track->Phi();
      
      if (TMath::Abs(trackphi - phi) > TMath::Abs(trackphi - phi + 2 * TMath::Pi()))
	trackphi += 2 * TMath::Pi();
      if (TMath::Abs(trackphi - phi) > TMath::Abs(trackphi - phi - 2 * TMath::Pi()))
	trackphi -= 2 * TMath::Pi();
      
      Float_t d = TMath::Sqrt((tracketa - eta) * (tracketa - eta) + (trackphi - phi) * (trackphi - phi));
      if (d <= fJetRadius){
	TVector3 vect(track->Px(), track->Py(), track->Pz());
	vect *= vect * rigidAxis / vect.Mag();
	ptrigid += vect.Pt();
	
	pt += track->Pt();
      }
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskSAJF::ExecOnce()
{
  // Initialize the analysis.

  if (!fEmbJetsName.IsNull()) {
    if (fEmbTracksName.IsNull()) {
      fEmbTracksName = fTracksName;
      fEmbTracksName += "Embedded";
    }
    if (fEmbCaloName.IsNull() && fAnaType == kEMCAL) {
      fEmbCaloName = fCaloName;
      fEmbCaloName += "Embedded";
    }
  }

  AliAnalysisTaskEmcalJet::ExecOnce();

  const Float_t semiDiag = TMath::Sqrt((fMaxPhi - fMinPhi) * (fMaxPhi - fMinPhi) + 
                                       (fMaxEta - fMinEta) * (fMaxEta - fMinEta)) / 2;
  if (fMinRC2LJ > semiDiag * 0.5) {
    AliWarning(Form("The parameter fMinRC2LJ = %f is too large for the considered acceptance. "
                    "Will use fMinRC2LJ = %f", fMinRC2LJ, semiDiag * 0.5));
    fMinRC2LJ = semiDiag * 0.5;
  }
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
