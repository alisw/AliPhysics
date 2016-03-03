/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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
//
// Analysis task for Dijet-hadron correlations
//
// Author: T.Kobayashi
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFormula.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile2D.h>
#include <THnSparse.h>
#include <TROOT.h>
#include <TTree.h>
#include <TArrayI.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TList.h>

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliVParticle.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliNamedArrayI.h"
#include "AliNamedString.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliEmcalJet.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

const Double_t pi = TMath::Pi();
//const Double_t areaCut[4] = {0.1, 0.23, 0.4, 0.63};

#include "AliAnalysisTaskDijetHadron.h"

ClassImp(AliAnalysisTaskDijetHadron)

//________________________________________________________________________
AliAnalysisTaskDijetHadron::AliAnalysisTaskDijetHadron() : AliAnalysisTaskEmcalJet("AliAnalysisTaskDijetHadron", kTRUE),
  fMCJetPtThreshold(1), fleadingHadronPtcut1(0.0), fleadingHadronPtcut2(3.0), fleadingHadronPtcut3(5.0), fJet1Ptcut1(10.0), fJet1Ptcut2(20.0), fJet1Ptcut3(30.0), fJet2Ptcut1(10.0), fJet2Ptcut2(20.0), fJet2Ptcut3(30.0), fConeRadius(0.2), fConeMinEta(-0.9), fConeMaxEta(0.9), fConeMinPhi(0), fConeMaxPhi(TMath::Pi()*2), fJetsCont(0), fTracksCont(0), fCaloClustersCont(0), fMCJetsCont(0), fMCTracksCont(0), fMCCaloClustersCont(0), fEmbJetsCont(0), fEmbTracksCont(0), fEmbCaloClustersCont(0), fCent_V0(0), fVertex_z_cut(0), fEP2(0), fJetBG_rho(0), fJetBG_rho_Cent(0), fTrackPt_PbPb(0), fTrackPhi_PbPb(0), fTrackEta_PbPb(0), fTrack_Phi_Eta_PbPb(0), fTrackPt_MC(0), fTrackPhi_MC(0), fTrackEta_MC(0), fTrack_Phi_Eta_MC(0), fTrackPt_EMB(0), fTrackPhi_EMB(0), fTrackEta_EMB(0), fTrack_Phi_Eta_EMB(0), fJetPt_PbPb(), fJetPhi_PbPb(), fJetEta_PbPb(), fJet_Phi_Eta_PbPb(), fJetPt_BG_PbPb(), fJetDeltaEP_PbPb(), fJet1Pt_PbPb(), fJet2Pt_PbPb(), fJet1Pt_BG_PbPb(), fJet2Pt_BG_PbPb(), fJetDeltaPhi_PbPb(), fJetDeltaEta_PbPb(), fJet1SelectPt_BG_PbPb(), fJet2SelectPt_BG_PbPb(), fJet1EP_PbPb(), fAj_PbPb(), fJetPt_MC(), fJetPhi_MC(), fJetEta_MC(), fJet_Phi_Eta_MC(), fJetDeltaEP_MC(), fJet1Pt_MC(), fJet2Pt_MC(), fJetDeltaPhi_MC(), fJetDeltaEta_MC(), fJet1EP_MC(), fAj_MC(), fJetPt_EMB(), fJetPhi_EMB(), fJetEta_EMB(), fJet_Phi_Eta_EMB(), fJetPt_BG_EMB(), fJetDeltaPt(), fJetDeltaEP_EMB(), fJet1Pt_EMB(), fJet2Pt_EMB(), fJet1Pt_BG_EMB(), fJet2Pt_BG_EMB(), fJet1DeltaPt(), fJet2DeltaPt(), fJetDeltaPhi_EMB(), fJetDeltaEta_EMB(), fJet1SelectPt_BG_EMB(), fJet2SelectPt_BG_EMB(), fJet1SelectDeltaPt(), fJet2SelectDeltaPt(), fJet1EP_EMB(), fAj_EMB(), fHJetDeltaPhi_Aj0_PbPb(), fHJetDeltaPhi_Aj1_PbPb(), fHJetDeltaPhi_Aj2_PbPb(), fHJetDeltaPhi_Aj3_PbPb(), fHJetDeltaPhi_Aj4_PbPb(), fHJet_EP_Aj0_PbPb(), fHJet_EP_Aj1_PbPb(), fHJet_EP_Aj2_PbPb(), fHJet_EP_Aj3_PbPb(), fHJet_EP_Aj4_PbPb(), fHJetDeltaPhi_Aj0_MC(), fHJetDeltaPhi_Aj1_MC(), fHJetDeltaPhi_Aj2_MC(), fHJetDeltaPhi_Aj3_MC(), fHJetDeltaPhi_Aj4_MC(), fHJet_EP_Aj0_MC(), fHJet_EP_Aj1_MC(), fHJet_EP_Aj2_MC(), fHJet_EP_Aj3_MC(), fHJet_EP_Aj4_MC(), fHJetDeltaPhi_Aj0_EMB(), fHJetDeltaPhi_Aj1_EMB(), fHJetDeltaPhi_Aj2_EMB(), fHJetDeltaPhi_Aj3_EMB(), fHJetDeltaPhi_Aj4_EMB(), fHJet_EP_Aj0_EMB(), fHJet_EP_Aj1_EMB(), fHJet_EP_Aj2_EMB(), fHJet_EP_Aj3_EMB(), fHJet_EP_Aj4_EMB(), fEvent(0), fCentrality(0) {

  // Default constructor.
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskDijetHadron::AliAnalysisTaskDijetHadron(const char *name) : AliAnalysisTaskEmcalJet(name, kTRUE),
  fMCJetPtThreshold(1), fleadingHadronPtcut1(0.0), fleadingHadronPtcut2(3.0), fleadingHadronPtcut3(5.0), fJet1Ptcut1(10.0), fJet1Ptcut2(20.0), fJet1Ptcut3(30.0), fJet2Ptcut1(10.0), fJet2Ptcut2(20.0), fJet2Ptcut3(30.0), fConeRadius(0.2), fConeMinEta(-0.9), fConeMaxEta(0.9), fConeMinPhi(0), fConeMaxPhi(TMath::Pi()*2), fJetsCont(0), fTracksCont(0), fCaloClustersCont(0), fMCJetsCont(0), fMCTracksCont(0), fMCCaloClustersCont(0), fEmbJetsCont(0), fEmbTracksCont(0), fEmbCaloClustersCont(0), fCent_V0(0), fVertex_z_cut(0), fEP2(0), fJetBG_rho(0), fJetBG_rho_Cent(0), fTrackPt_PbPb(0), fTrackPhi_PbPb(0), fTrackEta_PbPb(0), fTrack_Phi_Eta_PbPb(0), fTrackPt_MC(0), fTrackPhi_MC(0), fTrackEta_MC(0), fTrack_Phi_Eta_MC(0), fTrackPt_EMB(0), fTrackPhi_EMB(0), fTrackEta_EMB(0), fTrack_Phi_Eta_EMB(0), fJetPt_PbPb(), fJetPhi_PbPb(), fJetEta_PbPb(), fJet_Phi_Eta_PbPb(), fJetPt_BG_PbPb(), fJetDeltaEP_PbPb(), fJet1Pt_PbPb(), fJet2Pt_PbPb(), fJet1Pt_BG_PbPb(), fJet2Pt_BG_PbPb(), fJetDeltaPhi_PbPb(), fJetDeltaEta_PbPb(), fJet1SelectPt_BG_PbPb(), fJet2SelectPt_BG_PbPb(), fJet1EP_PbPb(), fAj_PbPb(), fJetPt_MC(), fJetPhi_MC(), fJetEta_MC(), fJet_Phi_Eta_MC(), fJetDeltaEP_MC(), fJet1Pt_MC(), fJet2Pt_MC(), fJetDeltaPhi_MC(), fJetDeltaEta_MC(), fJet1EP_MC(), fAj_MC(), fJetPt_EMB(), fJetPhi_EMB(), fJetEta_EMB(), fJet_Phi_Eta_EMB(), fJetPt_BG_EMB(), fJetDeltaPt(), fJetDeltaEP_EMB(), fJet1Pt_EMB(), fJet2Pt_EMB(), fJet1Pt_BG_EMB(), fJet2Pt_BG_EMB(), fJet1DeltaPt(), fJet2DeltaPt(), fJetDeltaPhi_EMB(), fJetDeltaEta_EMB(), fJet1SelectPt_BG_EMB(), fJet2SelectPt_BG_EMB(), fJet1SelectDeltaPt(), fJet2SelectDeltaPt(), fJet1EP_EMB(), fAj_EMB(), fHJetDeltaPhi_Aj0_PbPb(), fHJetDeltaPhi_Aj1_PbPb(), fHJetDeltaPhi_Aj2_PbPb(), fHJetDeltaPhi_Aj3_PbPb(), fHJetDeltaPhi_Aj4_PbPb(), fHJet_EP_Aj0_PbPb(), fHJet_EP_Aj1_PbPb(), fHJet_EP_Aj2_PbPb(), fHJet_EP_Aj3_PbPb(), fHJet_EP_Aj4_PbPb(), fHJetDeltaPhi_Aj0_MC(), fHJetDeltaPhi_Aj1_MC(), fHJetDeltaPhi_Aj2_MC(), fHJetDeltaPhi_Aj3_MC(), fHJetDeltaPhi_Aj4_MC(), fHJet_EP_Aj0_MC(), fHJet_EP_Aj1_MC(), fHJet_EP_Aj2_MC(), fHJet_EP_Aj3_MC(), fHJet_EP_Aj4_MC(), fHJetDeltaPhi_Aj0_EMB(), fHJetDeltaPhi_Aj1_EMB(), fHJetDeltaPhi_Aj2_EMB(), fHJetDeltaPhi_Aj3_EMB(), fHJetDeltaPhi_Aj4_EMB(), fHJet_EP_Aj0_EMB(), fHJet_EP_Aj1_EMB(), fHJet_EP_Aj2_EMB(), fHJet_EP_Aj3_EMB(), fHJet_EP_Aj4_EMB(), fEvent(0), fCentrality(0) {

  // Standard constructor.
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
void AliAnalysisTaskDijetHadron::AllocateHistogramArrays()
{
  fTrackPt_PbPb = new TH1*[fNcentBins]; fTrackPhi_PbPb = new TH1*[fNcentBins]; fTrackEta_PbPb = new TH1*[fNcentBins]; fTrack_Phi_Eta_PbPb = new TH2*[fNcentBins];
  fTrackPt_MC = new TH1*[fNcentBins]; fTrackPhi_MC = new TH1*[fNcentBins]; fTrackEta_MC = new TH1*[fNcentBins]; fTrack_Phi_Eta_MC = new TH2*[fNcentBins];
  fTrackPt_EMB = new TH1*[fNcentBins]; fTrackPhi_EMB = new TH1*[fNcentBins]; fTrackEta_EMB = new TH1*[fNcentBins]; fTrack_Phi_Eta_EMB = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
  fTrackPt_PbPb[i] = 0; fTrackPhi_PbPb[i] = 0; fTrackEta_PbPb[i] = 0; fTrack_Phi_Eta_PbPb[i] = 0;

  fTrackPt_MC[i] = 0; fTrackPhi_MC[i] = 0; fTrackEta_MC[i] = 0; fTrack_Phi_Eta_MC[i] = 0;

  fTrackPt_EMB[i] = 0; fTrackPhi_EMB[i] = 0; fTrackEta_EMB[i] = 0; fTrack_Phi_Eta_EMB[i] = 0;
  }

  for (Int_t i = 0; i < fNcentBins; i++) {
  for (Int_t j = 0; j < 3; j++) {
  fJetPt_PbPb[i][j] = 0; fJetPhi_PbPb[i][j] = 0; fJetEta_PbPb[i][j] = 0; fJet_Phi_Eta_PbPb[i][j] = 0; fJetPt_BG_PbPb[i][j] = 0; fJetDeltaEP_PbPb[i][j] = 0;

  fJetPt_MC[i][j] = 0; fJetPhi_MC[i][j] = 0; fJetEta_MC[i][j] = 0; fJet_Phi_Eta_MC[i][j] = 0; fJetDeltaEP_MC[i][j] = 0;

  fJetPt_EMB[i][j] = 0; fJetPhi_EMB[i][j] = 0; fJetEta_EMB[i][j] = 0; fJet_Phi_Eta_EMB[i][j] = 0; fJetPt_BG_EMB[i][j] = 0; fJetDeltaEP_EMB[i][j] = 0;
  }
  }

  for (Int_t i = 0; i < fNcentBins; i++) {
  for (Int_t j = 0; j < 3; j++) {
  for (Int_t k = 0; k < 4; k++) {
  for (Int_t l = 0; l < k+1; l++) {
  fJet1Pt_PbPb[i][j][k][l] = 0; fJet2Pt_PbPb[i][j][k][l] = 0; fJet1Pt_BG_PbPb[i][j][k][l] = 0; fJet2Pt_BG_PbPb[i][j][k][l] = 0; fJetDeltaPhi_PbPb[i][j][k][l] = 0; fJetDeltaEta_PbPb[i][j][k][l] = 0; fJet1SelectPt_BG_PbPb[i][j][k][l] = 0; fJet2SelectPt_BG_PbPb[i][j][k][l] = 0; fJet1EP_PbPb[i][j][k][l] = 0; fAj_PbPb[i][j][k][l] = 0;

  fJet1Pt_MC[i][j][k][l] = 0; fJet2Pt_MC[i][j][k][l] = 0; fJetDeltaPhi_MC[i][j][k][l] = 0; fJetDeltaEta_MC[i][j][k][l] = 0; fJet1EP_MC[i][j][k][l] = 0; fAj_MC[i][j][k][l] = 0;

  fJet1Pt_EMB[i][j][k][l] = 0; fJet2Pt_EMB[i][j][k][l] = 0; fJet1Pt_BG_EMB[i][j][k][l] = 0; fJet2Pt_BG_EMB[i][j][k][l] = 0; fJet1DeltaPt[i][j][k][l] = 0; fJet2DeltaPt[i][j][k][l] = 0; fJetDeltaPhi_EMB[i][j][k][l] = 0; fJetDeltaEta_EMB[i][j][k][l] = 0; fJet1SelectPt_BG_EMB[i][j][k][l] = 0; fJet2SelectPt_BG_EMB[i][j][k][l] = 0; fJet1SelectDeltaPt[i][j][k][l] = 0; fJet2SelectDeltaPt[i][j][k][l] = 0; fJet1EP_EMB[i][j][k][l] = 0; fAj_EMB[i][j][k][l] = 0;
  }
  }
  }
  }

  for (Int_t i = 0; i < fNcentBins; i++) {
  for (Int_t j = 0; j < 3; j++) {
  for (Int_t k = 0; k < 4; k++) {
  for (Int_t l = 0; l < 4; l++) {
  for (Int_t m = 0; m < l+1; m++) {
  fHJetDeltaPhi_Aj0_PbPb[i][j][k][l][m] = 0; fHJetDeltaPhi_Aj1_PbPb[i][j][k][l][m] = 0; fHJetDeltaPhi_Aj2_PbPb[i][j][k][l][m] = 0; fHJetDeltaPhi_Aj3_PbPb[i][j][k][l][m] = 0; fHJetDeltaPhi_Aj4_PbPb[i][j][k][l][m] = 0; fHJet_EP_Aj0_PbPb[i][j][k][l][m] = 0; fHJet_EP_Aj1_PbPb[i][j][k][l][m] = 0; fHJet_EP_Aj2_PbPb[i][j][k][l][m] = 0; fHJet_EP_Aj3_PbPb[i][j][k][l][m] = 0; fHJet_EP_Aj4_PbPb[i][j][k][l][m] = 0;

  fHJetDeltaPhi_Aj0_MC[i][j][k][l][m] = 0; fHJetDeltaPhi_Aj1_MC[i][j][k][l][m] = 0; fHJetDeltaPhi_Aj2_MC[i][j][k][l][m] = 0; fHJetDeltaPhi_Aj3_MC[i][j][k][l][m] = 0; fHJetDeltaPhi_Aj4_MC[i][j][k][l][m] = 0; fHJet_EP_Aj0_MC[i][j][k][l][m] = 0; fHJet_EP_Aj1_MC[i][j][k][l][m] = 0; fHJet_EP_Aj2_MC[i][j][k][l][m] = 0; fHJet_EP_Aj3_MC[i][j][k][l][m] = 0; fHJet_EP_Aj4_MC[i][j][k][l][m] = 0;

  fHJetDeltaPhi_Aj0_EMB[i][j][k][l][m] = 0; fHJetDeltaPhi_Aj1_EMB[i][j][k][l][m] = 0; fHJetDeltaPhi_Aj2_EMB[i][j][k][l][m] = 0; fHJetDeltaPhi_Aj3_EMB[i][j][k][l][m] = 0; fHJetDeltaPhi_Aj4_EMB[i][j][k][l][m] = 0; fHJet_EP_Aj0_EMB[i][j][k][l][m] = 0; fHJet_EP_Aj1_EMB[i][j][k][l][m] = 0; fHJet_EP_Aj2_EMB[i][j][k][l][m] = 0; fHJet_EP_Aj3_EMB[i][j][k][l][m] = 0; fHJet_EP_Aj4_EMB[i][j][k][l][m] = 0;
  }
  }
  }
  }
  }

}

//________________________________________________________________________
void AliAnalysisTaskDijetHadron::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  AllocateHistogramArrays();

  fJetsCont = GetJetContainer("Jets");
  fTracksCont = GetParticleContainer("Tracks");
  fCaloClustersCont = GetClusterContainer("CaloClusters");
  fMCJetsCont = GetJetContainer("MCJets");
  fMCTracksCont = GetParticleContainer("MCTracks");
  fMCCaloClustersCont = GetClusterContainer("MCCaloClusters");
  fEmbJetsCont = GetJetContainer("EmbJets");
  fEmbTracksCont = GetParticleContainer("EmbTracks");
  fEmbCaloClustersCont = GetClusterContainer("EmbCaloClusters");

   //User Task
  fCent_V0  = new TH1F("fCent_V0", "Centrality (all) by V0M", 103,-2,101);
  fOutput->Add(fCent_V0);
  fVertex_z_cut = new TH1F("fVertex_z_cut", "SPD vertex z (cut)", 120,-30,30);
  fOutput->Add(fVertex_z_cut);
  fEP2     = new TH1F("fEP2", "fEP2", 80,-pi,pi);
  fOutput->Add(fEP2);
  fJetBG_rho = new TH1F("fJetBG_rho","fJetBG_rho",300,0,300);
  fOutput->Add(fJetBG_rho);
  fJetBG_rho_Cent = new TH2F("fJetBG_rho_Cent","fJetBG_rho_Cent",100,0,100,300,0,300);
  fOutput->Add(fJetBG_rho_Cent);

  // Track histograms...
  for (Int_t i = 0; i < fNcentBins; i++) {
  //PbPb
  fTrackPt_PbPb[i]     = new TH1F(Form("fTrackPt_PbPb[%d]",i),Form("fTrackPt_PbPb[%d]",i),100,-80,120);
  fOutput->Add(fTrackPt_PbPb[i]);
  fTrackPhi_PbPb[i]     = new TH1F(Form("fTrackPhi_PbPb[%d]",i),Form("fTrackPhi_PbPb[%d]",i),40,0.0,2.*pi);
  fOutput->Add(fTrackPhi_PbPb[i]);
  fTrackEta_PbPb[i]     = new TH1F(Form("fTrackEta_PbPb[%d]",i),Form("fTrackEta_PbPb[%d]",i),40,-1.,1.);
  fOutput->Add(fTrackEta_PbPb[i]);
  fTrack_Phi_Eta_PbPb[i]     = new TH2F(Form("fTrack_Phi_Eta_PbPb[%d]",i),Form("fTrack_Phi_Eta_PbPb[%d]",i),40,0.0,2.*pi,40,-1.,1.);
  fOutput->Add(fTrack_Phi_Eta_PbPb[i]);

  //MC
  fTrackPt_MC[i]     = new TH1F(Form("fTrackPt_MC[%d]",i),Form("fTrackPt_MC[%d]",i),100,-80,120);
  fOutput->Add(fTrackPt_MC[i]);
  fTrackPhi_MC[i]     = new TH1F(Form("fTrackPhi_MC[%d]",i),Form("fTrackPhi_MC[%d]",i),40,0.0,2.*pi);
  fOutput->Add(fTrackPhi_MC[i]);
  fTrackEta_MC[i]     = new TH1F(Form("fTrackEta_MC[%d]",i),Form("fTrackEta_MC[%d]",i),40,-1.,1.);
  fOutput->Add(fTrackEta_MC[i]);
  fTrack_Phi_Eta_MC[i]     = new TH2F(Form("fTrack_Phi_Eta_MC[%d]",i),Form("fTrack_Phi_Eta_MC[%d]",i),40,0.0,2.*pi,40,-1.,1.);
  fOutput->Add(fTrack_Phi_Eta_MC[i]);

  //EMB
  fTrackPt_EMB[i]     = new TH1F(Form("fTrackPt_EMB[%d]",i),Form("fTrackPt_EMB[%d]",i),100,-80,120);
  fOutput->Add(fTrackPt_EMB[i]);
  fTrackPhi_EMB[i]     = new TH1F(Form("fTrackPhi_EMB[%d]",i),Form("fTrackPhi_EMB[%d]",i),40,0.0,2.*pi);
  fOutput->Add(fTrackPhi_EMB[i]);
  fTrackEta_EMB[i]     = new TH1F(Form("fTrackEta_EMB[%d]",i),Form("fTrackEta_EMB[%d]",i),40,-1.,1.);
  fOutput->Add(fTrackEta_EMB[i]);
  fTrack_Phi_Eta_EMB[i]     = new TH2F(Form("fTrack_Phi_Eta_EMB[%d]",i),Form("fTrack_Phi_Eta_EMB[%d]",i),40,0.0,2.*pi,40,-1.,1.);
  fOutput->Add(fTrack_Phi_Eta_EMB[i]);
  }

  for (Int_t i = 0; i < fNcentBins; i++) {
  for (Int_t j = 0; j < 3; j++) {
//Jet Histgrams...
  //PbPb
  fJetPt_PbPb[i][j] = new TH1F(Form("fJetPt_PbPb[%d][%d]",i,j),Form("fJetPt_PbPb[%d][%d]",i,j),100,-80,120);
  fOutput->Add(fJetPt_PbPb[i][j]);
  fJetPhi_PbPb[i][j] = new TH1F(Form("fJetPhi_PbPb[%d][%d]",i,j),Form("fJetPhi_PbPb[%d][%d]",i,j),40,0.0, 2*pi);
  fOutput->Add(fJetPhi_PbPb[i][j]);
  fJetEta_PbPb[i][j] = new TH1F(Form("fJetEta_PbPb[%d][%d]",i,j),Form("fJetEta_PbPb[%d][%d]",i,j),40,-1.,1.);
  fOutput->Add(fJetEta_PbPb[i][j]);
  fJet_Phi_Eta_PbPb[i][j] = new TH2F(Form("fJet_Phi_Eta_PbPb[%d][%d]",i,j),Form("fJet_Phi_Eta_PbPb[%d][%d]",i,j),40,0.0, 2*pi,40,-1.,1.);
  fOutput->Add(fJet_Phi_Eta_PbPb[i][j]);
  fJetPt_BG_PbPb[i][j] = new TH1F(Form("fJetPt_BG_PbPb[%d][%d]",i,j),Form("fJetPt_BG_PbPb[%d][%d]",i,j),100,-80,120);
  fOutput->Add(fJetPt_BG_PbPb[i][j]);
  fJetDeltaEP_PbPb[i][j] = new TH2F(Form("fJetDeltaEP_PbPb[%d][%d]",i,j),Form("fJetDeltaEP_PbPb[%d][%d]",i,j),40,-1./2.*pi,3./2.*pi,100,-80,120);
  fOutput->Add(fJetDeltaEP_PbPb[i][j]);

  //MC
  fJetPt_MC[i][j] = new TH1F(Form("fJetPt_MC[%d][%d]",i,j),Form("fJetPt_MC[%d][%d]",i,j),100,-80,120);
  fOutput->Add(fJetPt_MC[i][j]);
  fJetPhi_MC[i][j] = new TH1F(Form("fJetPhi_MC[%d][%d]",i,j),Form("fJetPhi_MC[%d][%d]",i,j),40,0.0, 2*pi);
  fOutput->Add(fJetPhi_MC[i][j]);
  fJetEta_MC[i][j] = new TH1F(Form("fJetEta_MC[%d][%d]",i,j),Form("fJetEta_MC[%d][%d]",i,j),40,-1.,1.);
  fOutput->Add(fJetEta_MC[i][j]);
  fJet_Phi_Eta_MC[i][j] = new TH2F(Form("fJet_Phi_Eta_MC[%d][%d]",i,j),Form("fJet_Phi_Eta_MC[%d][%d]",i,j),40,0.0, 2*pi,40,-1.,1.);
  fOutput->Add(fJet_Phi_Eta_MC[i][j]);
  fJetDeltaEP_MC[i][j] = new TH2F(Form("fJetDeltaEP_MC[%d][%d]",i,j),Form("fJetDeltaEP_MC[%d][%d]",i,j),40,-1./2.*pi,3./2.*pi,100,-80,120);
  fOutput->Add(fJetDeltaEP_MC[i][j]);

  //EMB
  fJetPt_EMB[i][j] = new TH1F(Form("fJetPt_EMB[%d][%d]",i,j),Form("fJetPt_EMB[%d][%d]",i,j),100,-80,120);
  fOutput->Add(fJetPt_EMB[i][j]);
  fJetPhi_EMB[i][j] = new TH1F(Form("fJetPhi_EMB[%d][%d]",i,j),Form("fJetPhi_EMB[%d][%d]",i,j),40,0.0, 2*pi);
  fOutput->Add(fJetPhi_EMB[i][j]);
  fJetEta_EMB[i][j] = new TH1F(Form("fJetEta_EMB[%d][%d]",i,j),Form("fJetEta_EMB[%d][%d]",i,j),40,-1.,1.);
  fOutput->Add(fJetEta_EMB[i][j]);
  fJet_Phi_Eta_EMB[i][j] = new TH2F(Form("fJet_Phi_Eta_EMB[%d][%d]",i,j),Form("fJet_Phi_Eta_EMB[%d][%d]",i,j),40,0.0, 2*pi,40,-1.,1.);
  fOutput->Add(fJet_Phi_Eta_EMB[i][j]);
  fJetPt_BG_EMB[i][j] = new TH1F(Form("fJetPt_BG_EMB[%d][%d]",i,j),Form("fJetPt_BG_EMB[%d][%d]",i,j),100,-80,120);
  fOutput->Add(fJetPt_BG_EMB[i][j]);
  fJetDeltaPt[i][j] = new TH1F(Form("fJetDeltaPt[%d][%d]",i,j),Form("fJetDeltaPt[%d][%d]",i,j),100,-80,120);
  fOutput->Add(fJetDeltaPt[i][j]);
  fJetDeltaEP_EMB[i][j] = new TH2F(Form("fJetDeltaEP_EMB[%d][%d]",i,j),Form("fJetDeltaEP_EMB[%d][%d]",i,j),40,-1./2.*pi,3./2.*pi,100,-80,120);
  fOutput->Add(fJetDeltaEP_EMB[i][j]);

 }
 }

  for (Int_t i = 0; i < fNcentBins; i++) {
  for (Int_t j = 0; j < 3; j++) {
  for (Int_t k = 0; k < 4; k++) {
  for (Int_t l = 0; l < k+1; l++) {
// Jet Histgrams...
  //PbPb
  fJet1Pt_PbPb[i][j][k][l] = new TH1F(Form("fJet1Pt_PbPb[%d][%d][%d][%d]",i,j,k,l),Form("fJet1Pt_PbPb[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet1Pt_PbPb[i][j][k][l]);
  fJet2Pt_PbPb[i][j][k][l] = new TH1F(Form("fJet2Pt_PbPb[%d][%d][%d][%d]",i,j,k,l),Form("fJet2Pt_PbPb[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet2Pt_PbPb[i][j][k][l]);
  fJet1Pt_BG_PbPb[i][j][k][l] = new TH1F(Form("fJet1Pt_BG_PbPb[%d][%d][%d][%d]",i,j,k,l),Form("fJet1Pt_BG_PbPb[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet1Pt_BG_PbPb[i][j][k][l]);
  fJet2Pt_BG_PbPb[i][j][k][l] = new TH1F(Form("fJet2Pt_BG_PbPb[%d][%d][%d][%d]",i,j,k,l),Form("fJet2Pt_BG_PbPb[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet2Pt_BG_PbPb[i][j][k][l]);
  fJetDeltaPhi_PbPb[i][j][k][l] = new TH1F(Form("fJetDeltaPhi_PbPb[%d][%d][%d][%d]",i,j,k,l),Form("fJetDeltaPhi_PbPb[%d][%d][%d][%d]",i,j,k,l),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fJetDeltaPhi_PbPb[i][j][k][l]);
  fJetDeltaEta_PbPb[i][j][k][l] = new TH1F(Form("fJetDeltaEta_PbPb[%d][%d][%d][%d]",i,j,k,l),Form("fJetDeltaEta_PbPb[%d][%d][%d][%d]",i,j,k,l),40,-1.0,1.0);
  fOutput->Add(fJetDeltaEta_PbPb[i][j][k][l]);
  fJet1SelectPt_BG_PbPb[i][j][k][l] = new TH1F(Form("fJet1SelectPt_BG_PbPb[%d][%d][%d][%d]",i,j,k,l),Form("fJet1SelectPt_BG_PbPb[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet1SelectPt_BG_PbPb[i][j][k][l]);
  fJet2SelectPt_BG_PbPb[i][j][k][l] = new TH1F(Form("fJet2SelectPt_BG_PbPb[%d][%d][%d][%d]",i,j,k,l),Form("fJet2SelectPt_BG_PbPb[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet2SelectPt_BG_PbPb[i][j][k][l]);
  fJet1EP_PbPb[i][j][k][l] = new TH2F(Form("fJet1EP_PbPb[%d][%d][%d][%d]",i,j,k,l),Form("fJet1EP_PbPb[%d][%d][%d][%d]",i,j,k,l),40,-1./2.*pi,3./2.*pi,100,-80,120);
  fOutput->Add(fJet1EP_PbPb[i][j][k][l]);
  fAj_PbPb[i][j][k][l] = new TH1F(Form("fAj_PbPb[%d][%d][%d][%d]",i,j,k,l),Form("fAj_PbPb[%d][%d][%d][%d]",i,j,k,l),20,0.,1.);
  fOutput->Add(fAj_PbPb[i][j][k][l]);

  //MC
  fJet1Pt_MC[i][j][k][l] = new TH1F(Form("fJet1Pt_MC[%d][%d][%d][%d]",i,j,k,l),Form("fJet1Pt_MC[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet1Pt_MC[i][j][k][l]);
  fJet2Pt_MC[i][j][k][l] = new TH1F(Form("fJet2Pt_MC[%d][%d][%d][%d]",i,j,k,l),Form("fJet2Pt_MC[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet2Pt_MC[i][j][k][l]);
  fJetDeltaPhi_MC[i][j][k][l] = new TH1F(Form("fJetDeltaPhi_MC[%d][%d][%d][%d]",i,j,k,l),Form("fJetDeltaPhi_MC[%d][%d][%d][%d]",i,j,k,l),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fJetDeltaPhi_MC[i][j][k][l]);
  fJetDeltaEta_MC[i][j][k][l] = new TH1F(Form("fJetDeltaEta_MC[%d][%d][%d][%d]",i,j,k,l),Form("fJetDeltaEta_MC[%d][%d][%d][%d]",i,j,k,l),40,-1.0,1.0);
  fOutput->Add(fJetDeltaEta_MC[i][j][k][l]);
  fJet1EP_MC[i][j][k][l] = new TH2F(Form("fJet1EP_MC[%d][%d][%d][%d]",i,j,k,l),Form("fJet1EP_MC[%d][%d][%d][%d]",i,j,k,l),40,-1./2.*pi,3./2.*pi,100,-80,120);
  fOutput->Add(fJet1EP_MC[i][j][k][l]);
  fAj_MC[i][j][k][l] = new TH1F(Form("fAj_MC[%d][%d][%d][%d]",i,j,k,l),Form("fAj_MC[%d][%d][%d][%d]",i,j,k,l),20,0.,1.);
  fOutput->Add(fAj_MC[i][j][k][l]);

  //EMB
  fJet1Pt_EMB[i][j][k][l] = new TH1F(Form("fJet1Pt_EMB[%d][%d][%d][%d]",i,j,k,l),Form("fJet1Pt_EMB[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet1Pt_EMB[i][j][k][l]);
  fJet2Pt_EMB[i][j][k][l] = new TH1F(Form("fJet2Pt_EMB[%d][%d][%d][%d]",i,j,k,l),Form("fJet2Pt_EMB[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet2Pt_EMB[i][j][k][l]);
  fJet1Pt_BG_EMB[i][j][k][l] = new TH1F(Form("fJet1Pt_BG_EMB[%d][%d][%d][%d]",i,j,k,l),Form("fJet1Pt_BG_EMB[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet1Pt_BG_EMB[i][j][k][l]);
  fJet2Pt_BG_EMB[i][j][k][l] = new TH1F(Form("fJet2Pt_BG_EMB[%d][%d][%d][%d]",i,j,k,l),Form("fJet2Pt_BG_EMB[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet2Pt_BG_EMB[i][j][k][l]);
  fJet1DeltaPt[i][j][k][l] = new TH1F(Form("fJet1DeltaPt[%d][%d][%d][%d]",i,j,k,l),Form("fJet1DeltaPt[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet1DeltaPt[i][j][k][l]);
  fJet2DeltaPt[i][j][k][l] = new TH1F(Form("fJet2DeltaPt[%d][%d][%d][%d]",i,j,k,l),Form("fJet2DeltaPt[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet2DeltaPt[i][j][k][l]);
  fJetDeltaPhi_EMB[i][j][k][l] = new TH1F(Form("fJetDeltaPhi_EMB[%d][%d][%d][%d]",i,j,k,l),Form("fJetDeltaPhi_EMB[%d][%d][%d][%d]",i,j,k,l),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fJetDeltaPhi_EMB[i][j][k][l]);
  fJetDeltaEta_EMB[i][j][k][l] = new TH1F(Form("fJetDeltaEta_EMB[%d][%d][%d][%d]",i,j,k,l),Form("fJetDeltaEta_EMB[%d][%d][%d][%d]",i,j,k,l),40,-1.0,1.0);
  fOutput->Add(fJetDeltaEta_EMB[i][j][k][l]);
  fJet1SelectPt_BG_EMB[i][j][k][l] = new TH1F(Form("fJet1SelectPt_BG_EMB[%d][%d][%d][%d]",i,j,k,l),Form("fJet1SelectPt_BG_EMB[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet1SelectPt_BG_EMB[i][j][k][l]);
  fJet2SelectPt_BG_EMB[i][j][k][l] = new TH1F(Form("fJet2SelectPt_BG_EMB[%d][%d][%d][%d]",i,j,k,l),Form("fJet2SelectPt_BG_EMB[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet2SelectPt_BG_EMB[i][j][k][l]);
  fJet1SelectDeltaPt[i][j][k][l] = new TH1F(Form("fJet1SelectDeltaPt[%d][%d][%d][%d]",i,j,k,l),Form("fJet1SelectDeltaPt[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet1SelectDeltaPt[i][j][k][l]);
  fJet2SelectDeltaPt[i][j][k][l] = new TH1F(Form("fJet2SelectDeltaPt[%d][%d][%d][%d]",i,j,k,l),Form("fJet2SelectDeltaPt[%d][%d][%d][%d]",i,j,k,l),100,-80,120);
  fOutput->Add(fJet2SelectDeltaPt[i][j][k][l]);
  fJet1EP_EMB[i][j][k][l] = new TH2F(Form("fJet1EP_EMB[%d][%d][%d][%d]",i,j,k,l),Form("fJet1EP_EMB[%d][%d][%d][%d]",i,j,k,l),40,-1./2.*pi,3./2.*pi,100,-80,120);
  fOutput->Add(fJet1EP_EMB[i][j][k][l]);
  fAj_EMB[i][j][k][l] = new TH1F(Form("fAj_EMB[%d][%d][%d][%d]",i,j,k,l),Form("fAj_EMB[%d][%d][%d][%d]",i,j,k,l),20,0.,1.);
  fOutput->Add(fAj_EMB[i][j][k][l]);


 }
 }
 }
 }

  for (Int_t i = 0; i < fNcentBins; i++) {
  for (Int_t j = 0; j < 3; j++) {
  for (Int_t k = 0; k < 4; k++) {
  for (Int_t l = 0; l < 4; l++) {
  for (Int_t m = 0; m < l+1; m++) {
// Jet-Hadron Histgrams...
  //PbPb
  fHJetDeltaPhi_Aj0_PbPb[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj0_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj0_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj0_PbPb[i][j][k][l][m]);
  fHJetDeltaPhi_Aj1_PbPb[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj1_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj1_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj1_PbPb[i][j][k][l][m]);
  fHJetDeltaPhi_Aj2_PbPb[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj2_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj2_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj2_PbPb[i][j][k][l][m]);
  fHJetDeltaPhi_Aj3_PbPb[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj3_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj3_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj3_PbPb[i][j][k][l][m]);
  fHJetDeltaPhi_Aj4_PbPb[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj4_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj4_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj4_PbPb[i][j][k][l][m]);

  fHJet_EP_Aj0_PbPb[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj0_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj0_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj0_PbPb[i][j][k][l][m]);
  fHJet_EP_Aj1_PbPb[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj1_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj1_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj1_PbPb[i][j][k][l][m]);
  fHJet_EP_Aj2_PbPb[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj2_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj2_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj2_PbPb[i][j][k][l][m]);
  fHJet_EP_Aj3_PbPb[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj3_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj3_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj3_PbPb[i][j][k][l][m]);
  fHJet_EP_Aj4_PbPb[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj4_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj4_PbPb[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj4_PbPb[i][j][k][l][m]);

  //MC
  fHJetDeltaPhi_Aj0_MC[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj0_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj0_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj0_MC[i][j][k][l][m]);
  fHJetDeltaPhi_Aj1_MC[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj1_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj1_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj1_MC[i][j][k][l][m]);
  fHJetDeltaPhi_Aj2_MC[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj2_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj2_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj2_MC[i][j][k][l][m]);
  fHJetDeltaPhi_Aj3_MC[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj3_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj3_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj3_MC[i][j][k][l][m]);
  fHJetDeltaPhi_Aj4_MC[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj4_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj4_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj4_MC[i][j][k][l][m]);

  fHJet_EP_Aj0_MC[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj0_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj0_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj0_MC[i][j][k][l][m]);
  fHJet_EP_Aj1_MC[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj1_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj1_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj1_MC[i][j][k][l][m]);
  fHJet_EP_Aj2_MC[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj2_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj2_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj2_MC[i][j][k][l][m]);
  fHJet_EP_Aj3_MC[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj3_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj3_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj3_MC[i][j][k][l][m]);
  fHJet_EP_Aj4_MC[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj4_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj4_MC[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj4_MC[i][j][k][l][m]);

  //EMB
  fHJetDeltaPhi_Aj0_EMB[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj0_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj0_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj0_EMB[i][j][k][l][m]);
  fHJetDeltaPhi_Aj1_EMB[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj1_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj1_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj1_EMB[i][j][k][l][m]);
  fHJetDeltaPhi_Aj2_EMB[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj2_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj2_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj2_EMB[i][j][k][l][m]);
  fHJetDeltaPhi_Aj3_EMB[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj3_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj3_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj3_EMB[i][j][k][l][m]);
  fHJetDeltaPhi_Aj4_EMB[i][j][k][l][m] = new TH1F(Form("fHJetDeltaPhi_Aj4_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJetDeltaPhi_Aj4_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJetDeltaPhi_Aj4_EMB[i][j][k][l][m]);

  fHJet_EP_Aj0_EMB[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj0_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj0_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj0_EMB[i][j][k][l][m]);
  fHJet_EP_Aj1_EMB[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj1_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj1_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj1_EMB[i][j][k][l][m]);
  fHJet_EP_Aj2_EMB[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj2_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj2_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj2_EMB[i][j][k][l][m]);
  fHJet_EP_Aj3_EMB[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj3_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj3_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj3_EMB[i][j][k][l][m]);
  fHJet_EP_Aj4_EMB[i][j][k][l][m] = new TH1F(Form("fHJet_EP_Aj4_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),Form("fHJet_EP_Aj4_EMB[%d][%d][%d][%d][%d]",i,j,k,l,m),40,-1./2.*pi,3./2.*pi);
  fOutput->Add(fHJet_EP_Aj4_EMB[i][j][k][l][m]);
  }
  }
  }
  }
  }

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
Bool_t AliAnalysisTaskDijetHadron::FillHistograms()
{

  // Fill histograms.
  // Avoid TPCHole(lhc11h)
  Int_t runNumber = InputEvent()->GetRunNumber();
  Int_t fAvoidTpcHole = 0;
      Int_t runs_iroc[28] = {169975, 169981, 170038, 170040, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170163, 170193, 170195, 170203, 170204, 170205, 170228, 170230, 170264, 170268, 170269, 170270, 170306, 170308, 170309};
      for(Int_t i=0; i<28; i++)
	{
	  if(runNumber==runs_iroc[i])
	    {
	      fAvoidTpcHole = 1;
	      break;
	    }
	}

  //if(fAvoidTpcHole==1 && !(trigPhi>3.89 && trigPhi<5.53)) trigIndex = -1;

  //Get Event
  fEvent = InputEvent();
  if(fEvent && fAvoidTpcHole == 0){

  //trigger
  Int_t fTriggerType =-1;
  UInt_t trigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

      if (trigger & AliVEvent::kAnyINT)      { fTriggerType=0; }
      else if (trigger & AliVEvent::kCentral)     { fTriggerType=0; }
      else if (trigger & AliVEvent::kSemiCentral) { fTriggerType=0; }
      else if (trigger & AliVEvent::kEMCEGA)      { fTriggerType=1; }
      else if (trigger & AliVEvent::kEMCEJE)      { fTriggerType=2; }

  // Vertex cut 
  const AliVVertex* vtx = fEvent->GetPrimaryVertex();

  if(fTriggerType==0){
  if (vtx && vtx->GetNContributors()>1.){
  if (TMath::Abs(vtx->GetZ())<10.){
  fVertex_z_cut->Fill(vtx->GetZ());

  // GetCentrality
      AliCentrality *centrality = fEvent->GetCentrality();
      if (centrality)
      	fCentrality = centrality->GetCentralityPercentile("V0M");
      else 
	fCentrality = 99;
  fCent_V0->Fill(fCentrality);

  // Get EP(2) and Rhovalue
  Double_t EP_2 = fEPV0;
  fEP2->Fill(EP_2);
  fJetBG_rho->Fill(fRhoVal);
  fJetBG_rho_Cent->Fill(fCentrality,fRhoVal);

  // Prepare Jet value
//PbPb
  // jet counting value
  Int_t c_jet1_PbPb[3] ={0,0,0}; Int_t c_jet2_PbPb[3] ={0,0,0};
  Int_t leading_jet_count0[3]={0,0,0}; Int_t subleading_jet_count0[3]={0,0,0};
  Double_t dEPJet0[4]  ={-999.,-999.,-999.,-999.};
  // Jet1: Leading jet
  Double_t jet1_pt0[4]  ={-999.,-999.,-999.,-999.}; Double_t jet1_pt_BG0[4]  ={-999.,-999.,-999.,-999.};
  Double_t jet1_phi0[4]  ={-999.,-999.,-999.,-999.}; Double_t jet1_eta0[4]  ={-999.,-999.,-999.,-999.};
  // Jet2: Sub Leading jet
  Double_t jet2_pt0[4]  ={-999.,-999.,-999.,-999.}; Double_t jet2_pt_BG0[4]  ={-999.,-999.,-999.,-999.};
  Double_t jet2_phi0[4]  ={-999.,-999.,-999.,-999.}; Double_t jet2_eta0[4]  ={-999.,-999.,-999.,-999.};
  // Correlation of jet1 and jet2
  Double_t Delta_phi0[4]  ={-999.,-999.,-999.,-999.}; Double_t Delta_epjet1_0[4] ={-999.,-999.,-999.,-999.}; Double_t Delta_eta0[4]  ={-999.,-999.,-999.,-999.}; Double_t Aj0[4]  ={-999.,-999.,-999.,-999.};
//MC
  // jet counting value
  Int_t c_jet1_MC[3] ={0,0,0}; Int_t c_jet2_MC[3] ={0,0,0};
  Int_t leading_jet_count1[3]={0,0,0}; Int_t subleading_jet_count1[3]={0,0,0};
  Double_t dEPJet1[4]  ={-999.,-999.,-999.,-999.};
  // Jet1: Leading jet
  Double_t jet1_pt1[4]  ={-999.,-999.,-999.,-999.};
  Double_t jet1_phi1[4]  ={-999.,-999.,-999.,-999.}; Double_t jet1_eta1[4]  ={-999.,-999.,-999.,-999.};
  // Jet2: Sub Leading jet
  Double_t jet2_pt1[4]  ={-999.,-999.,-999.,-999.};
  Double_t jet2_phi1[4]  ={-999.,-999.,-999.,-999.}; Double_t jet2_eta1[4]  ={-999.,-999.,-999.,-999.};
  // Correlation of jet1 and jet2
  Double_t Delta_phi1[4]  ={-999.,-999.,-999.,-999.}; Double_t Delta_epjet1_1[4] ={-999.,-999.,-999.,-999.}; Double_t Delta_eta1[4]  ={-999.,-999.,-999.,-999.}; Double_t Aj1[4]  ={-999.,-999.,-999.,-999.};
//EMB
  // jet counting value
  Int_t c_jet1_EMB[3] ={0,0,0}; Int_t c_jet2_EMB[3] ={0,0,0};
  Int_t leading_jet_count2[3]={0,0,0}; Int_t subleading_jet_count2[3]={0,0,0};
  Double_t dEPJet2[4]  ={-999.,-999.,-999.,-999.};
  // Jet1: Leading jet
  Double_t jet1_pt2[4]  ={-999.,-999.,-999.,-999.}; Double_t jet1_pt_BG2[4]  ={-999.,-999.,-999.,-999.}; Double_t jet1_Deltapt[4]  ={-999.,-999.,-999.,-999.};
  Double_t jet1_phi2[4]  ={-999.,-999.,-999.,-999.}; Double_t jet1_eta2[4]  ={-999.,-999.,-999.,-999.};
  // Jet2: Sub Leading jet
  Double_t jet2_pt2[4]  ={-999.,-999.,-999.,-999.}; Double_t jet2_pt_BG2[4]  ={-999.,-999.,-999.,-999.}; Double_t jet2_Deltapt[4]  ={-999.,-999.,-999.,-999.};
  Double_t jet2_phi2[4]  ={-999.,-999.,-999.,-999.}; Double_t jet2_eta2[4]  ={-999.,-999.,-999.,-999.};
  // Correlation of jet1 and jet2
  Double_t Delta_phi2[4]  ={-999.,-999.,-999.,-999.}; Double_t Delta_epjet1_2[4] ={-999.,-999.,-999.,-999.}; Double_t Delta_eta2[4]  ={-999.,-999.,-999.,-999.}; Double_t Aj2[4]  ={-999.,-999.,-999.,-999.};


  //threshold
  double Jet1_threshold[4]; double Jet2_threshold[4];
  Jet1_threshold[0]=0.0; Jet2_threshold[0]=0.0;
  Jet1_threshold[1]=fJet1Ptcut1; Jet2_threshold[1]=fJet2Ptcut1;
  Jet1_threshold[2]=fJet1Ptcut2; Jet2_threshold[2]=fJet2Ptcut2;
  Jet1_threshold[3]=fJet1Ptcut3; Jet2_threshold[3]=fJet2Ptcut3;

  // ************
  // PbPb
  // _________________________________

  //Track histogram
  if (fTracksCont) {
    fTracksCont->ResetCurrentID();
    AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()); 
    while(track) {
      if(track->GetLabel()==0){
      fTrackPt_PbPb[fCentBin]->Fill(track->Pt()); 
      fTrackPhi_PbPb[fCentBin]->Fill(track->Phi()); 
      fTrackEta_PbPb[fCentBin]->Fill(track->Eta()); 
      fTrack_Phi_Eta_PbPb[fCentBin]->Fill(track->Phi(),track->Eta());
      }
      track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    }
  }

  // Jet and Jet-HadronHistgram
  if (fJetsCont) {
    fJetsCont->ResetCurrentID();
    AliEmcalJet *Jet = fJetsCont->GetNextAcceptJet();
    while(Jet) {

      //leading track cut
      Int_t leading_track_count[3]={0,0,0};
      if(fJetsCont->GetLeadingHadronPt(Jet) > fleadingHadronPtcut1) leading_track_count[0] += 1;
      if(fJetsCont->GetLeadingHadronPt(Jet) > fleadingHadronPtcut2) leading_track_count[1] += 1;
      if(fJetsCont->GetLeadingHadronPt(Jet) > fleadingHadronPtcut3) leading_track_count[2] += 1;

      for(int m=0;m<3;m++){
      c_jet1_PbPb[m]++; // jet count in acceptance.
      if(leading_track_count[m] > 0){
      dEPJet0[m] = GetDPhi(Jet->Phi(),EP_2);
      fJetPt_PbPb[fCentBin][m]->Fill(Jet->Pt());
      fJetPhi_PbPb[fCentBin][m]->Fill(Jet->Phi());
      fJetEta_PbPb[fCentBin][m]->Fill(Jet->Eta());
      fJet_Phi_Eta_PbPb[fCentBin][m]->Fill(Jet->Phi(),Jet->Eta());
      fJetPt_BG_PbPb[fCentBin][m]->Fill(Jet->Pt() - Jet->Area() * fJetsCont->GetRhoVal());
      fJetDeltaEP_PbPb[fCentBin][m]->Fill(dEPJet0[m],Jet->Pt());
      }
      }

for(int m=0;m<3;m++){
      if(c_jet1_PbPb[m] == 1)
	{
        if(leading_track_count[m] > 0){
          jet1_pt0[m] = Jet->Pt(); jet1_pt_BG0[m] = Jet->Pt() - Jet->Area() * fJetsCont->GetRhoVal(); jet1_phi0[m] = Jet->Phi(); jet1_eta0[m] = Jet->Eta(); //Get Leading Jet(Jet1) value
        leading_jet_count0[m] += 1;
        }
	}

      else if(c_jet1_PbPb[m] > 1 && c_jet2_PbPb[m] == 0  && leading_jet_count0[m] > 0 && leading_track_count[m] > 0)// sub leading
	{
	  jet2_pt0[m] = Jet->Pt(); jet2_pt_BG0[m] = Jet->Pt() - Jet->Area() * fJetsCont->GetRhoVal(); jet2_phi0[m] = Jet->Phi(); jet2_eta0[m] = Jet->Eta(); //Get Sub Leading Jet(Jet2) value
          Delta_phi0[m] = GetDPhi(jet1_phi0[m],jet2_phi0[m]); Aj0[m] = (jet1_pt_BG0[m] - jet2_pt_BG0[m]) / (jet1_pt_BG0[m] + jet2_pt_BG0[m]);
	  Delta_epjet1_0[m] = GetDPhi(jet1_phi0[m],EP_2);
	  Delta_eta0[m] = jet1_eta0[m] - jet2_eta0[m]; //Get Correlations of jet1, jet2

	 // Correlations of jet1, jet2
	for(int count1=0;count1<4;count1++){
		for(int count2=0;count2<count1+1;count2++){
		if(jet1_pt_BG0[m] > Jet1_threshold[count1] && jet2_pt_BG0[m] > Jet2_threshold[count2]){
          fJet1Pt_PbPb[fCentBin][m][count1][count2]->Fill(jet1_pt0[m]);
          fJet2Pt_PbPb[fCentBin][m][count1][count2]->Fill(jet2_pt0[m]);
          fJet1Pt_BG_PbPb[fCentBin][m][count1][count2]->Fill(jet1_pt_BG0[m]);
          fJet2Pt_BG_PbPb[fCentBin][m][count1][count2]->Fill(jet2_pt_BG0[m]);
          fJetDeltaPhi_PbPb[fCentBin][m][count1][count2]->Fill(Delta_phi0[m]);
          fJetDeltaEta_PbPb[fCentBin][m][count1][count2]->Fill(Delta_eta0[m]);
		}
	}//count2
	}//count1


          // Find delta_phi
          if(Delta_phi0[m] > (2./3.)* TMath::Pi() && Delta_phi0[m] < (4./3.)* TMath::Pi()){
	        for(int count1=0;count1<4;count1++){
		for(int count2=0;count2<count1+1;count2++){
		if(jet1_pt_BG0[m] > Jet1_threshold[count1] && jet2_pt_BG0[m] > Jet2_threshold[count2]){
                fJet1SelectPt_BG_PbPb[fCentBin][m][count1][count2]->Fill(jet1_pt_BG0[m]);
                fJet2SelectPt_BG_PbPb[fCentBin][m][count1][count2]->Fill(jet2_pt_BG0[m]);
                fJet1EP_PbPb[fCentBin][m][count1][count2]->Fill(Delta_epjet1_0[m],jet1_pt_BG0[m]);
                fAj_PbPb[fCentBin][m][count1][count2]->Fill((jet1_pt_BG0[m] - jet2_pt_BG0[m]) / (jet1_pt_BG0[m] + jet2_pt_BG0[m]));
		        }
	        }//count2
	        }//count1
          }// Find delta_phi

	c_jet2_PbPb[m]++;
	subleading_jet_count0[m] += 1;
         }// sub leading

  }

  Jet = fJetsCont->GetNextAcceptJet();
  }// jet while

//jet-hadron
for(int m=0;m<3;m++){//leading track cut
  if (fTracksCont) {//track cont
	int c_subleading_jet = 0;
	c_subleading_jet = subleading_jet_count0[m];
	if(c_subleading_jet > 0){//if find sub leading

    fTracksCont->ResetCurrentID();
    AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    while(track) {
      if(track->GetLabel()==0){

      Double_t pt,phi,dphi,dep;
      pt = -999.0; phi = -999.0; dphi = -999.0; dep = -999.0;
      pt = track->Pt(); phi = track->Phi();
      dphi = GetDPhi(jet1_phi0[m],phi);
      dep = GetDPhi(EP_2,phi);

      //dphi cut
      Bool_t dphi_cut[4];
      //devide(jet side)
      dphi_cut[0]= ((-1./3.*TMath::Pi()) <= dphi && dphi <= (1./3.*TMath::Pi()));
      dphi_cut[1]= ((2./3.*TMath::Pi()) <= dphi && dphi <= (4./3.*TMath::Pi()));
      //devide(out jet side)
      dphi_cut[2]= (((1./3.*TMath::Pi()) <= dphi && dphi <= (1./2.*TMath::Pi())) || ((-1./2.*TMath::Pi()) < dphi && dphi < (-1./3.*TMath::Pi())));//leadingjet
      dphi_cut[3]= (((1./2.*TMath::Pi()) <= dphi && dphi <= (2./3.*TMath::Pi())) || ((4./3.*TMath::Pi()) < dphi && dphi < (3./2.*TMath::Pi())));//subleadingjet

      //pt switch
      Bool_t pt_switch[4];
      pt_switch[0]= (pt > 0.15);
      pt_switch[1]= (pt > 0.15 && pt <= 2.0);
      pt_switch[2]= (pt > 2.0 && pt <= 4.0);
      pt_switch[3]= (pt > 4.0);

      //jetdphi switch
      Double_t jet_dphi = -999.0;
      jet_dphi = Delta_phi0[m];
      Bool_t jet_dphi_switch[3];
      jet_dphi_switch[0]= (jet_dphi > (2./3.)* TMath::Pi() && jet_dphi < (4./3.)* TMath::Pi());
      jet_dphi_switch[1]= (jet_dphi > (5./6.)* TMath::Pi() && jet_dphi < (7./12.)* TMath::Pi());
      jet_dphi_switch[2]= (jet_dphi > (11./12.)* TMath::Pi() && jet_dphi < (13./12.)* TMath::Pi());

//hadron-dphi
for(int pt_cut=0;pt_cut<4;pt_cut++){
if(pt_switch[pt_cut]){
		for(int count1=0;count1<4;count1++){
		if(jet_dphi_switch[0]){
			for(int count2=0;count2<count1+1;count2++){
        		if(jet1_pt_BG0[m] > Jet1_threshold[count1] && jet2_pt_BG0[m] > Jet2_threshold[count2]){
				fHJetDeltaPhi_Aj0_PbPb[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj0_PbPb[fCentBin][m][pt_cut][count1][count2]->Fill(dep);

                		if(Aj0[m] >= 0.0 && Aj0[m] < 0.2){
				fHJetDeltaPhi_Aj1_PbPb[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj1_PbPb[fCentBin][m][pt_cut][count1][count2]->Fill(dep);
                        	}

                		if(Aj0[m] >= 0.2 && Aj0[m] < 0.4){
				fHJetDeltaPhi_Aj2_PbPb[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj2_PbPb[fCentBin][m][pt_cut][count1][count2]->Fill(dep);
                        	}

                		if(Aj0[m] >= 0.4 && Aj0[m] < 0.6){
				fHJetDeltaPhi_Aj3_PbPb[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj3_PbPb[fCentBin][m][pt_cut][count1][count2]->Fill(dep);
                        	}

                		if(Aj0[m] >= 0.6 && Aj0[m] <= 0.8){
				fHJetDeltaPhi_Aj4_PbPb[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj4_PbPb[fCentBin][m][pt_cut][count1][count2]->Fill(dep);
                        	}
                	}


		}//count2
		}
        	}//count1
}//pt cut
}//pt for

      }// if label
      track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    }// track while

}// if sub leading jet
}// tracks Cont
}// jet leading track cut

       }// jetCont

  // ************
  // MC
  // _________________________________

  //Track histogram
  if (fMCTracksCont) {
    fMCTracksCont->ResetCurrentID();
    AliVTrack *MCtrack = static_cast<AliVTrack*>(fMCTracksCont->GetNextAcceptParticle()); 
    while(MCtrack) {
      if(MCtrack->GetLabel()!=0){
      fTrackPt_MC[fCentBin]->Fill(MCtrack->Pt()); 
      fTrackPhi_MC[fCentBin]->Fill(MCtrack->Phi()); 
      fTrackEta_MC[fCentBin]->Fill(MCtrack->Eta()); 
      fTrack_Phi_Eta_MC[fCentBin]->Fill(MCtrack->Phi(),MCtrack->Eta());
      }
      MCtrack = static_cast<AliVTrack*>(fMCTracksCont->GetNextAcceptParticle());
    }
  }

  // Jet and Jet-HadronHistgram
  if (fMCJetsCont) {
    fMCJetsCont->ResetCurrentID();
    AliEmcalJet *MCJet = fMCJetsCont->GetNextAcceptJet();
    while(MCJet) {

      //leading track cut
      Int_t leading_track_count[3]={0,0,0};
      if(fMCJetsCont->GetLeadingHadronPt(MCJet) > fleadingHadronPtcut1) leading_track_count[0] += 1;
      if(fMCJetsCont->GetLeadingHadronPt(MCJet) > fleadingHadronPtcut2) leading_track_count[1] += 1;
      if(fMCJetsCont->GetLeadingHadronPt(MCJet) > fleadingHadronPtcut3) leading_track_count[2] += 1;

      for(int m=0;m<3;m++){
      c_jet1_MC[m]++; // jet count in acceptance.
      if(leading_track_count[m] > 0){
      dEPJet1[m] = GetDPhi(MCJet->Phi(),EP_2);
      fJetPt_MC[fCentBin][m]->Fill(MCJet->Pt());
      fJetPhi_MC[fCentBin][m]->Fill(MCJet->Phi());
      fJetEta_MC[fCentBin][m]->Fill(MCJet->Eta());
      fJet_Phi_Eta_MC[fCentBin][m]->Fill(MCJet->Phi(),MCJet->Eta());
      fJetDeltaEP_MC[fCentBin][m]->Fill(dEPJet1[m],MCJet->Pt());
      }
      }

for(int m=0;m<3;m++){
      if(c_jet1_MC[m] == 1)
	{
        if(leading_track_count[m] > 0){
          jet1_pt1[m] = MCJet->Pt(); jet1_phi1[m] = MCJet->Phi(); jet1_eta1[m] = MCJet->Eta(); //Get Leading Jet(Jet1) value

        leading_jet_count1[m] += 1;
        }
	}

      else if(c_jet1_MC[m] > 1 && c_jet2_MC[m] == 0  && leading_jet_count1[m] > 0 && leading_track_count[m] > 0)// sub leading
	{
	  jet2_pt1[m] = MCJet->Pt(); jet2_phi1[m] = MCJet->Phi(); jet2_eta1[m] = MCJet->Eta(); //Get Sub Leading Jet(Jet2) value
          Delta_phi1[m] = GetDPhi(jet1_phi1[m],jet2_phi1[m]); Aj1[m] = (jet1_pt1[m] - jet2_pt1[m]) / (jet1_pt1[m] + jet2_pt1[m]);
	  Delta_epjet1_1[m] = GetDPhi(jet1_phi1[m],EP_2);
	  Delta_eta1[m] = jet1_eta1[m] - jet2_eta1[m]; //Get Correlations of jet1, jet2

	 // Correlations of jet1, jet2
	for(int count1=0;count1<4;count1++){
		for(int count2=0;count2<count1+1;count2++){
		if(jet1_pt1[m] > Jet1_threshold[count1] && jet2_pt1[m] > Jet2_threshold[count2]){
          fJet1Pt_MC[fCentBin][m][count1][count2]->Fill(jet1_pt1[m]);
          fJet2Pt_MC[fCentBin][m][count1][count2]->Fill(jet2_pt1[m]);
          fJetDeltaPhi_MC[fCentBin][m][count1][count2]->Fill(Delta_phi1[m]);
          fJetDeltaEta_MC[fCentBin][m][count1][count2]->Fill(Delta_eta1[m]);
		}
	}//count2
	}//count1


          // Find delta_phi
          if(Delta_phi1[m] > (2./3.)* TMath::Pi() && Delta_phi1[m] < (4./3.)* TMath::Pi()){
	        for(int count1=0;count1<4;count1++){
		for(int count2=0;count2<count1+1;count2++){
		if(jet1_pt1[m] > Jet1_threshold[count1] && jet2_pt1[m] > Jet2_threshold[count2]){
                fJet1EP_MC[fCentBin][m][count1][count2]->Fill(Delta_epjet1_1[m],jet1_pt1[m]);
                fAj_MC[fCentBin][m][count1][count2]->Fill((jet1_pt1[m] - jet2_pt1[m]) / (jet1_pt1[m] + jet2_pt1[m]));
		        }
	        }//count2
	        }//count1
          }// Find delta_phi

	c_jet2_MC[m]++;
	subleading_jet_count1[m] += 1;
         }// sub leading

  }

MCJet = fMCJetsCont->GetNextAcceptJet(); 
}//while jet

//jet-hadron
for(int m=0;m<3;m++){
  if (fMCTracksCont) {//track cont
	int c_subleading_jet = 0;
	c_subleading_jet = subleading_jet_count1[m];
	if(c_subleading_jet > 0){//if find sub leading

    fMCTracksCont->ResetCurrentID();
    AliVTrack *MCtrack = static_cast<AliVTrack*>(fMCTracksCont->GetNextAcceptParticle());
    while(MCtrack) {
      if(MCtrack->GetLabel()!=0){

      Double_t pt,phi,dphi,dep;
      pt = -999.0; phi = -999.0; dphi = -999.0; dep = -999.0;
      pt = MCtrack->Pt(); phi = MCtrack->Phi();
      dphi = GetDPhi(jet1_phi1[m],phi);
      dep = GetDPhi(EP_2,phi);

      //dphi cut
      Bool_t dphi_cut[4];
      //devide(jet side)
      dphi_cut[0]= ((-1./3.*TMath::Pi()) <= dphi && dphi <= (1./3.*TMath::Pi()));
      dphi_cut[1]= ((2./3.*TMath::Pi()) <= dphi && dphi <= (4./3.*TMath::Pi()));
      //devide(out jet side)
      dphi_cut[2]= (((1./3.*TMath::Pi()) <= dphi && dphi <= (1./2.*TMath::Pi())) || ((-1./2.*TMath::Pi()) < dphi && dphi < (-1./3.*TMath::Pi())));//leadingjet
      dphi_cut[3]= (((1./2.*TMath::Pi()) <= dphi && dphi <= (2./3.*TMath::Pi())) || ((4./3.*TMath::Pi()) < dphi && dphi < (3./2.*TMath::Pi())));//subleadingjet

      //pt switch
      Bool_t pt_switch[4];
      pt_switch[0]= (pt > 0.15);
      pt_switch[1]= (pt > 0.15 && pt <= 2.0);
      pt_switch[2]= (pt > 2.0 && pt <= 4.0);
      pt_switch[3]= (pt > 4.0);

      //jetdphi switch
      Double_t jet_dphi = -999.0;
      jet_dphi = Delta_phi1[m];
      Bool_t jet_dphi_switch[3];
      jet_dphi_switch[0]= (jet_dphi > (2./3.)* TMath::Pi() && jet_dphi < (4./3.)* TMath::Pi());
      jet_dphi_switch[1]= (jet_dphi > (5./6.)* TMath::Pi() && jet_dphi < (7./12.)* TMath::Pi());
      jet_dphi_switch[2]= (jet_dphi > (11./12.)* TMath::Pi() && jet_dphi < (13./12.)* TMath::Pi());

//hadron-dphi
for(int pt_cut=0;pt_cut<4;pt_cut++){
if(pt_switch[pt_cut]){
		for(int count1=0;count1<4;count1++){
		if(jet_dphi_switch[0]){
			for(int count2=0;count2<count1+1;count2++){
        		if(jet1_pt1[m] > Jet1_threshold[count1] && jet2_pt1[m] > Jet2_threshold[count2]){
				fHJetDeltaPhi_Aj0_MC[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj0_MC[fCentBin][m][pt_cut][count1][count2]->Fill(dep);

                		if(Aj1[m] >= 0.0 && Aj1[m] < 0.2){
				fHJetDeltaPhi_Aj1_MC[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj1_MC[fCentBin][m][pt_cut][count1][count2]->Fill(dep);
                        	}

                		if(Aj1[m] >= 0.2 && Aj1[m] < 0.4){
				fHJetDeltaPhi_Aj2_MC[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj2_MC[fCentBin][m][pt_cut][count1][count2]->Fill(dep);
                        	}

                		if(Aj1[m] >= 0.4 && Aj1[m] < 0.6){
				fHJetDeltaPhi_Aj3_MC[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj3_MC[fCentBin][m][pt_cut][count1][count2]->Fill(dep);
                        	}

                		if(Aj1[m] >= 0.6 && Aj1[m] <= 0.8){
				fHJetDeltaPhi_Aj4_MC[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj4_MC[fCentBin][m][pt_cut][count1][count2]->Fill(dep);
                        	}
                	}


		}//count2
		}
        	}//count1
}//pt cut
}//pt for

      }// if label
      MCtrack = static_cast<AliVTrack*>(fMCTracksCont->GetNextAcceptParticle());
    }// track while

}// if sub leading jet
}// tracks Cont
}// jet leading track cut

}//jet Cont

  // ************
  // Embedding
  // _________________________________

  //Track histogram
  if (fEmbTracksCont) {
    fEmbTracksCont->ResetCurrentID();
    AliVTrack *EMBtrack = static_cast<AliVTrack*>(fEmbTracksCont->GetNextAcceptParticle()); 
    while(EMBtrack) {
      fTrackPt_EMB[fCentBin]->Fill(EMBtrack->Pt()); 
      fTrackPhi_EMB[fCentBin]->Fill(EMBtrack->Phi()); 
      fTrackEta_EMB[fCentBin]->Fill(EMBtrack->Eta()); 
      fTrack_Phi_Eta_EMB[fCentBin]->Fill(EMBtrack->Phi(),EMBtrack->Eta());
      EMBtrack = static_cast<AliVTrack*>(fEmbTracksCont->GetNextAcceptParticle());
    }
  }

  if (fEmbJetsCont) {
    
    AliEmcalJet *embJet = NextEmbeddedJet(kTRUE);
    
    while (embJet != 0) {

      //leading track cut
      Int_t leading_track_count[3]={0,0,0};
      if(fEmbJetsCont->GetLeadingHadronPt(embJet) > fleadingHadronPtcut1) leading_track_count[0] += 1;
      if(fEmbJetsCont->GetLeadingHadronPt(embJet) > fleadingHadronPtcut2) leading_track_count[1] += 1;
      if(fEmbJetsCont->GetLeadingHadronPt(embJet) > fleadingHadronPtcut3) leading_track_count[2] += 1;

      for(int m=0;m<3;m++){
      c_jet1_EMB[m]++; // jet count in acceptance.
      if(leading_track_count[m] > 0){
      dEPJet2[m] = GetDPhi(embJet->Phi(),EP_2);
      fJetPt_EMB[fCentBin][m]->Fill(embJet->Pt());
      fJetPhi_EMB[fCentBin][m]->Fill(embJet->Phi());
      fJetEta_EMB[fCentBin][m]->Fill(embJet->Eta());
      fJet_Phi_Eta_EMB[fCentBin][m]->Fill(embJet->Phi(),embJet->Eta());
      fJetPt_BG_EMB[fCentBin][m]->Fill(embJet->Pt() - embJet->Area() * fEmbJetsCont->GetRhoVal());
      fJetDeltaPt[fCentBin][m]->Fill(embJet->Pt() - embJet->Area() * fEmbJetsCont->GetRhoVal() - embJet->MCPt());
      fJetDeltaEP_EMB[fCentBin][m]->Fill(dEPJet2[m],embJet->Pt());
      }
      }

for(int m=0;m<3;m++){
      if(c_jet1_EMB[m] == 1)
	{
        if(leading_track_count[m] > 0){
          jet1_pt2[m] = embJet->Pt(); jet1_pt_BG2[m] = embJet->Pt() - embJet->Area() * fEmbJetsCont->GetRhoVal(); jet1_Deltapt[m] = embJet->Pt() - embJet->Area() * fEmbJetsCont->GetRhoVal() - embJet->MCPt(); jet1_phi2[m] = embJet->Phi(); jet1_eta2[m] = embJet->Eta(); //Get Leading Jet(Jet1) value

        leading_jet_count2[m] += 1;
        }
	}

      else if(c_jet1_EMB[m] > 1 && c_jet2_EMB[m] == 0  && leading_jet_count2[m] > 0 && leading_track_count[m] > 0)// sub leading
	{
	  jet2_pt2[m] = embJet->Pt(); jet2_pt_BG2[m] = embJet->Pt() - embJet->Area() * fEmbJetsCont->GetRhoVal(); jet2_Deltapt[m] = embJet->Pt() - embJet->Area() * fEmbJetsCont->GetRhoVal() - embJet->MCPt(); jet2_phi2[m] = embJet->Phi(); jet2_eta2[m] = embJet->Eta(); //Get Sub Leading Jet(Jet2) value
          Delta_phi2[m] = GetDPhi(jet1_phi2[m],jet2_phi2[m]); Aj2[m] = (jet1_pt_BG2[m] - jet2_pt_BG2[m]) / (jet1_pt_BG2[m] + jet2_pt_BG2[m]);
	  Delta_epjet1_2[m] = GetDPhi(jet1_phi2[m],EP_2);
	  Delta_eta2[m] = jet1_eta2[m] - jet2_eta2[m]; //Get Correlations of jet1, jet2

	 // Correlations of jet1, jet2
	for(int count1=0;count1<4;count1++){
		for(int count2=0;count2<count1+1;count2++){
		if(jet1_pt_BG2[m] > Jet1_threshold[count1] && jet2_pt_BG2[m] > Jet2_threshold[count2]){
          fJet1Pt_EMB[fCentBin][m][count1][count2]->Fill(jet1_pt2[m]);
          fJet2Pt_EMB[fCentBin][m][count1][count2]->Fill(jet2_pt2[m]);
          fJet1Pt_BG_EMB[fCentBin][m][count1][count2]->Fill(jet1_pt_BG2[m]);
          fJet2Pt_BG_EMB[fCentBin][m][count1][count2]->Fill(jet2_pt_BG2[m]);
          fJet1DeltaPt[fCentBin][m][count1][count2]->Fill(jet1_Deltapt[m]);
          fJet2DeltaPt[fCentBin][m][count1][count2]->Fill(jet2_Deltapt[m]);
          fJetDeltaPhi_EMB[fCentBin][m][count1][count2]->Fill(Delta_phi2[m]);
          fJetDeltaEta_EMB[fCentBin][m][count1][count2]->Fill(Delta_eta2[m]);
		}
	}//count2
	}//count1


          // Find delta_phi
          if(Delta_phi2[m] > (2./3.)* TMath::Pi() && Delta_phi2[m] < (4./3.)* TMath::Pi()){
	        for(int count1=0;count1<4;count1++){
		for(int count2=0;count2<count1+1;count2++){
		if(jet1_pt_BG2[m] > Jet1_threshold[count1] && jet2_pt_BG2[m] > Jet2_threshold[count2]){
                fJet1SelectPt_BG_EMB[fCentBin][m][count1][count2]->Fill(jet1_pt_BG2[m]);
                fJet2SelectPt_BG_EMB[fCentBin][m][count1][count2]->Fill(jet2_pt_BG2[m]);
                fJet1SelectDeltaPt[fCentBin][m][count1][count2]->Fill(jet1_Deltapt[m]);
                fJet2SelectDeltaPt[fCentBin][m][count1][count2]->Fill(jet2_Deltapt[m]);
                fJet1EP_EMB[fCentBin][m][count1][count2]->Fill(Delta_epjet1_2[m],jet1_pt_BG2[m]);
                fAj_EMB[fCentBin][m][count1][count2]->Fill((jet1_pt_BG2[m] - jet2_pt_BG2[m]) / (jet1_pt_BG2[m] + jet2_pt_BG2[m]));
		        }
	        }//count2
	        }//count1
          }// Find delta_phi

	c_jet2_EMB[m]++;
	subleading_jet_count2[m] += 1;
         }// sub leading

  }

      embJet = NextEmbeddedJet();
    }// jet while


//jet-hadron
for(int m=0;m<3;m++){
  if (fEmbTracksCont) {//track cont
	int c_subleading_jet = 0;
	c_subleading_jet = subleading_jet_count2[m];
	if(c_subleading_jet > 0){//if find sub leading

    fEmbTracksCont->ResetCurrentID();
    AliVTrack *EMBtrack = static_cast<AliVTrack*>(fEmbTracksCont->GetNextAcceptParticle());
    while(EMBtrack) {

      Double_t pt,phi,dphi,dep;
      pt = -999.0; phi = -999.0; dphi = -999.0; dep = -999.0;
      pt = EMBtrack->Pt(); phi = EMBtrack->Phi();
      dphi = GetDPhi(jet1_phi2[m],phi);
      dep = GetDPhi(EP_2,phi);

      //dphi cut
      Bool_t dphi_cut[4];
      //devide(jet side)
      dphi_cut[0]= ((-1./3.*TMath::Pi()) <= dphi && dphi <= (1./3.*TMath::Pi()));
      dphi_cut[1]= ((2./3.*TMath::Pi()) <= dphi && dphi <= (4./3.*TMath::Pi()));
      //devide(out jet side)
      dphi_cut[2]= (((1./3.*TMath::Pi()) <= dphi && dphi <= (1./2.*TMath::Pi())) || ((-1./2.*TMath::Pi()) < dphi && dphi < (-1./3.*TMath::Pi())));//leadingjet
      dphi_cut[3]= (((1./2.*TMath::Pi()) <= dphi && dphi <= (2./3.*TMath::Pi())) || ((4./3.*TMath::Pi()) < dphi && dphi < (3./2.*TMath::Pi())));//subleadingjet

      //pt switch
      Bool_t pt_switch[4];
      pt_switch[0]= (pt > 0.15);
      pt_switch[1]= (pt > 0.15 && pt <= 2.0);
      pt_switch[2]= (pt > 2.0 && pt <= 4.0);
      pt_switch[3]= (pt > 4.0);

      //jetdphi switch
      Double_t jet_dphi = -999.0;
      jet_dphi = Delta_phi2[m];
      Bool_t jet_dphi_switch[3];
      jet_dphi_switch[0]= (jet_dphi > (2./3.)* TMath::Pi() && jet_dphi < (4./3.)* TMath::Pi());
      jet_dphi_switch[1]= (jet_dphi > (5./6.)* TMath::Pi() && jet_dphi < (7./12.)* TMath::Pi());
      jet_dphi_switch[2]= (jet_dphi > (11./12.)* TMath::Pi() && jet_dphi < (13./12.)* TMath::Pi());

//hadron-dphi
for(int pt_cut=0;pt_cut<4;pt_cut++){
if(pt_switch[pt_cut]){
		for(int count1=0;count1<4;count1++){
		if(jet_dphi_switch[0]){
			for(int count2=0;count2<count1+1;count2++){
        		if(jet1_pt_BG2[m] > Jet1_threshold[count1] && jet2_pt_BG2[m] > Jet2_threshold[count2]){
				fHJetDeltaPhi_Aj0_EMB[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj0_EMB[fCentBin][m][pt_cut][count1][count2]->Fill(dep);

                		if(Aj2[m] >= 0.0 && Aj2[m] < 0.2){
				fHJetDeltaPhi_Aj1_EMB[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj1_EMB[fCentBin][m][pt_cut][count1][count2]->Fill(dep);
                        	}

                		if(Aj2[m] >= 0.2 && Aj2[m] < 0.4){
				fHJetDeltaPhi_Aj2_EMB[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj2_EMB[fCentBin][m][pt_cut][count1][count2]->Fill(dep);
                        	}

                		if(Aj2[m] >= 0.4 && Aj2[m] < 0.6){
				fHJetDeltaPhi_Aj3_EMB[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj3_EMB[fCentBin][m][pt_cut][count1][count2]->Fill(dep);
                        	}

                		if(Aj2[m] >= 0.6 && Aj2[m] <= 0.8){
				fHJetDeltaPhi_Aj4_EMB[fCentBin][m][pt_cut][count1][count2]->Fill(dphi);
				fHJet_EP_Aj4_EMB[fCentBin][m][pt_cut][count1][count2]->Fill(dep);
                        	}
                	}


		}//count2
		}
        	}//count1
}//pt cut
}//pt for

      EMBtrack = static_cast<AliVTrack*>(fEmbTracksCont->GetNextAcceptParticle());
    }// track while

}// if sub leading jet
}// tracks Cont
}// jet leading track cut

  }// jet Cont


}//vertex
}//vertex
}//trigger

}//event
  return kTRUE;
}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskDijetHadron::NextEmbeddedJet(Bool_t reset)
{
  // Get the next accepted embedded jet.

  if(reset)
    fEmbJetsCont->ResetCurrentID();
      
  AliEmcalJet* jet = fEmbJetsCont->GetNextAcceptJet();
  while (jet && jet->MCPt() < fMCJetPtThreshold) jet = fEmbJetsCont->GetNextAcceptJet();

  return jet;
}

//________________________________________________________________________
void AliAnalysisTaskDijetHadron::SetConeEtaPhiEMCAL()
{
  // Set default cuts for full cones

  SetConeEtaLimits(-0.7+fConeRadius,0.7-fConeRadius);
  SetConePhiLimits(1.4+fConeRadius,TMath::Pi()-fConeRadius);
}

//________________________________________________________________________
void AliAnalysisTaskDijetHadron::SetConeEtaPhiTPC()
{
  // Set default cuts for charged cones

  SetConeEtaLimits(-0.9+fConeRadius, 0.9-fConeRadius);
  SetConePhiLimits(-10, 10);
}

//________________________________________________________________________
void AliAnalysisTaskDijetHadron::ExecOnce()
{
  // Initialize the analysis.

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;
  if (fEmbTracksCont && fEmbTracksCont->GetArray() == 0) fEmbTracksCont = 0;
  if (fEmbCaloClustersCont && fEmbCaloClustersCont->GetArray() == 0) fEmbCaloClustersCont = 0;
  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fMCJetsCont && fMCJetsCont->GetArray() == 0) fMCJetsCont = 0;
  if (fEmbJetsCont && fEmbJetsCont->GetArray() == 0) fEmbJetsCont = 0;
}

//________________________________________________________________________
Double_t AliAnalysisTaskDijetHadron::GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz, const Double_t jetPx, const Double_t jetPy, const Double_t jetPz)
{
  return (trkPx*jetPx+trkPy*jetPy+trkPz*jetPz)/(jetPx*jetPx+jetPy*jetPy+jetPz*jetPz);
}

//________________________________________________________________________
Double_t AliAnalysisTaskDijetHadron::GetDPhi(Double_t mphi,Double_t vphi)
{
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  double delta_phi = mphi-vphi;
  if (delta_phi < (-1./2*TMath::Pi())) delta_phi += (2*TMath::Pi());
  else if (delta_phi > (3./2*TMath::Pi())) delta_phi -= (2*TMath::Pi());

  return delta_phi;//delta_phi in [-1/2*Pi, 3/2*Pi]                                                                                                    
}

