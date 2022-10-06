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

// Analysis task for antineutron measurement in PHOS and mesurment of charged Sigma bar
// Authors: Pavel Gordeev, Dmitri Blau, Dmitri Peresunko
#include "AliAnalysisSigmaBarCharged.h"

#include <cstdio>
#include <iostream>

#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliCaloPhoton.h"
#include "AliCascadeVertexer.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"
#include "AliPHOSGeometry.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliTriggerAnalysis.h"
#include "AliVertexerTracks.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "THashList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TParticle.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TTree.h"

ClassImp(AliAnalysisSigmaBarCharged);

AliAnalysisSigmaBarCharged::AliAnalysisSigmaBarCharged(const char* name)
  : AliAnalysisTaskSE(name),
    fOutputContainer(nullptr),
    fPIDResponse(nullptr),
    fGamma(nullptr),
    fPi(nullptr),
    fCentBin(0),
    fNCenBin(5),
    fCurrentMixedList(nullptr),
    fCluMinECut(0.3),
    fCluTimeCut(25.e-9),
    fCluNbarMinE(1.),
    fOptDepth(10.),
    fCPACut(0.9),
    fCPVCut(4.),
    fDispCut(2.),
    fDispA(-1),
    fDispB(3.5),
    fTracksBits(21),
    fIsMC(false),
    fPHOSClusterTOFOption(0),
    fPHOSClusterEnergyOption(0),
    fStack(nullptr),
    fEvent(nullptr)
{
  for (Int_t i = 0; i < 1; i++)
    for (Int_t j = 0; j < 5; j++)
      fPHOSEvents[i][j] = 0x0; // Container for PHOS photons
  fCenBinEdges.Set(fNCenBin);
  for (int cen = 1; cen <= fNCenBin; cen++)
    fCenBinEdges.AddAt(int(100. * cen / fNCenBin), cen - 1);
  DefineOutput(1, THashList::Class());
}

AliAnalysisSigmaBarCharged::AliAnalysisSigmaBarCharged(const AliAnalysisSigmaBarCharged& rh)
  : AliAnalysisTaskSE(rh.GetName()),
    fOutputContainer(nullptr),
    fPIDResponse(nullptr),
    fGamma(nullptr),
    fPi(nullptr),
    fCentBin(0),
    fNCenBin(5),
    fCurrentMixedList(nullptr),
    fCluMinECut(0.3),
    fCluTimeCut(25.e-9),
    fCluNbarMinE(1.),
    fOptDepth(10.),
    fCPACut(0.9),
    fCPVCut(4.),
    fDispCut(2.),
    fDispA(-1),
    fDispB(3.5),
    fTracksBits(21),
    fIsMC(false),
    fPHOSClusterTOFOption(0),
    fPHOSClusterEnergyOption(0),
    fStack(nullptr),
    fEvent(nullptr)
{
  for (Int_t i = 0; i < 1; i++)
    for (Int_t j = 0; j < 5; j++)
      fPHOSEvents[i][j] = 0x0; // Container for PHOS photons
  fCenBinEdges.Set(rh.fCenBinEdges.GetSize(), rh.fCenBinEdges.GetArray());
  if (fOutputContainer)
    delete fOutputContainer;
  fOutputContainer = new THashList();
}

AliAnalysisSigmaBarCharged& AliAnalysisSigmaBarCharged::operator=(const AliAnalysisSigmaBarCharged& rh)
{
  this->~AliAnalysisSigmaBarCharged();
  new (this) AliAnalysisSigmaBarCharged(rh);
  return *this;
}

AliAnalysisSigmaBarCharged::~AliAnalysisSigmaBarCharged()
{
  if (fOutputContainer) {
    delete fOutputContainer;
    fOutputContainer = 0x0;
  }
  for (Int_t i = 0; i < 10; i++)
    for (Int_t j = 0; j < 2; j++)
      if (fPHOSEvents[i][j]) {
        delete fPHOSEvents[i][j];
        fPHOSEvents[i][j] = 0x0;
      }
}

void AliAnalysisSigmaBarCharged::UserCreateOutputObjects()
{
  if (fOutputContainer != nullptr) {
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);
  // General QC

  // Criteria
  char cPID[4][15];
  snprintf(cPID[0], 15, "%s", ""); // to avoid gcc warnings
  snprintf(cPID[1], 15, "CPV2_");
  snprintf(cPID[2], 15, "Disp2_");
  snprintf(cPID[3], 15, "CPV2_Disp2_");

  char cIM[5][15];
  snprintf(cIM[0], 15, "%s", "");
  snprintf(cIM[1], 15, "CPV2_");
  snprintf(cIM[2], 15, "Disp2_");
  snprintf(cIM[3], 15, "All3_");
  snprintf(cIM[4], 15, "All4_");

  char cMC[2][15];
  snprintf(cMC[0], 15, "%s", "");
  snprintf(cMC[1], 15, "MC_");

  char cCH[6][30];
  snprintf(cCH[0], 30, "_Charge1");
  snprintf(cCH[1], 30, "_Charge-1");
  snprintf(cCH[2], 30, "_AntiSigmaPlus");
  snprintf(cCH[3], 30, "_AntiSigmaMinus");
  snprintf(cCH[4], 30, "_AntiSigmaPlus_ParentCheck");
  snprintf(cCH[5], 30, "_AntiSigmaMinus_ParentCheck");

  // Binning
  int dispmin = 0, dispmax = 15, dispbins = 200; // dispersion: m02 and m20
  double tofmin = -250.e-9, tofmax = 250.e-9;
  int tofbins = 500; // Time Of Flight
  double dcaxymin = -2.4, dcaxymax = 2.4, dcazmin = -3.2, dcazmax = 3.2;
  int dcabins = 800; // DCA XY and DCA Z
  double invmin = 1.0, invmax = 1.5;
  int invbins = 1000;                           // invariant mass
  int recmin = -20, recmax = 20, recbins = 400; // reconstructed momentum
  int dcamin = 0, dcamax = 15;                  // dca between daughters
  int cpamin = -1, cpamax = 1, cpabins = 400;   // cosine of pointing angle
  int radmin = 0, radmax = 15, radbins = 800;   // distance between primary and secondary vertex
  int ebins = 400;
  float emin = 0, emax = 20.;

  // Dispersion and check of dispersion cuts
  fOutputContainer->Add(new TH3F("All_Disp_MinCut", "Dispersion of all particles min cut", dispbins, dispmin, dispmax,
                                 dispbins, dispmin, dispmax, ebins, emin, emax));
  fOutputContainer->Add(new TH3F("AntiNeutron_Disp_MinCut", "Dispersion of antineutron min cut", dispbins, dispmin,
                                 dispmax, dispbins, dispmin, dispmax, ebins, emin, emax));
  fOutputContainer->Add(new TH3F("All_Disp", "Dispersion of all particles", dispbins, dispmin, dispmax, dispbins,
                                 dispmin, dispmax, ebins, emin, emax));
  fOutputContainer->Add(new TH3F("AntiNeutron_Disp", "Dispersion of antineutron", dispbins, dispmin, dispmax, dispbins,
                                 dispmin, dispmax, ebins, emin, emax));
  fOutputContainer->Add(new TH3F("Disp2_Cut", "Dispersion of all after cut", dispbins, dispmin, dispmax, dispbins,
                                 dispmin, dispmax, ebins, emin, emax));
  fOutputContainer->Add(new TH3F("AntiNeutron_Disp2_Cut", "Dispersion of all after cut", dispbins, dispmin, dispmax,
                                 dispbins, dispmin, dispmax, ebins, emin, emax));

  // Relation between cluster and MC energy
  fOutputContainer->Add(new TH2F("EnergyclustervsMC", "Cluster vs MC energy", ebins, emin, emax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("AntiNeutron_EnergyclustervsMC", "Antineutron cluster vs MC energy", ebins, emin, emax,
                                 ebins, emin, emax));

  // Spectra of identified clusters with CPV и Disp cuts; antineutron spectra
  fOutputContainer->Add(new TH1F("Spectrum_Cluster_TOF", "Spectrum of clusters", ebins, emin, emax));

  for (int iter = 0; iter < 8; iter++) {
    fOutputContainer->Add(new TH1F(Form("MC_A_%d", iter), "MC spectrum of antineutron", ebins, emin, emax));
    if (iter > 0) {
      fOutputContainer->Add(new TH1F(Form("Spectrum_Cluster_%d", iter), "Spectrum of clusters", ebins, emin, emax));
    }
    if (iter > 2) {
      fOutputContainer->Add(
        new TH1F(Form("MC_APlus_%d", iter + 1), "MC spectrum of antineutron from AntiSigmaPlus", ebins, emin, emax));
      fOutputContainer->Add(new TH1F(Form("Rec_APlus_%d", iter + 1),
                                     "Reconstructed spectrum of antineutron from AntiSigmaPlus", ebins, emin, emax));
      fOutputContainer->Add(
        new TH1F(Form("MC_AMinus_%d", iter + 1), "MC spectrum of antineutron from AntiSigmaMinus", ebins, emin, emax));
      fOutputContainer->Add(new TH1F(Form("Rec_AMinus_%d", iter + 1),
                                     "Reconstructed spectrum of antineutron from AntiSigmaMinus", ebins, emin, emax));
    }
    if (iter > 3) {
      fOutputContainer->Add(
        new TH1F(Form("Rec_A_%d", iter), "Reconstructed spectrum of antineutron", ebins, emin, emax));
      fOutputContainer->Add(
        new TH1F(Form("Rec_Spectrum_Cluster_%d", iter), "Reconstructed spectrum of clusters", ebins, emin, emax));
    }
  }

  // Track spectra; pion spectra
  fOutputContainer->Add(new TH1F("MC_PionPlus_0", "MC spectrum of pions with pos charge", ebins, emin, emax));
  fOutputContainer->Add(new TH1F("MC_PionMinus_0", "MC spectrum of pions with neg charge", ebins, emin, emax));

  for (int iter = 1; iter < 6; iter++) {
    if (iter < 3) {
      fOutputContainer->Add(
        new TH1F(Form("TrackPlus_%d", iter), "Spectrum of tracks with pos charge", ebins, emin, emax));
      fOutputContainer->Add(
        new TH1F(Form("TrackMinus_%d", iter), "Spectrum of tracks with neg charge", ebins, emin, emax));
    }
    fOutputContainer->Add(new TH1F(Form("PionPlus_%d", iter), "Spectrum of pions with pos charge", ebins, emin, emax));
    fOutputContainer->Add(new TH1F(Form("PionMinus_%d", iter), "Spectrum of pions with neg charge", ebins, emin, emax));
  }

  // AntiSigma spectra
  for (int iter = 0; iter < 8; iter++) {
    fOutputContainer->Add(new TH1F(Form("MC_AntiSigmaPlus_%d", iter), "Spectrum of AntiSigmaPlus", ebins, emin, emax));
    fOutputContainer->Add(new TH1F(Form("MC_AntiSigmaMinus_%d", iter), "Spectrum of AntiSigmaPlus", ebins, emin, emax));
  }

  fOutputContainer->Add(new TH1F("MC_AntiSigmaPlus_5_CPV2", "Spectrum of AntiSigmaPlus", ebins, emin, emax));
  fOutputContainer->Add(new TH1F("MC_AntiSigmaPlus_5_Disp2", "Spectrum of AntiSigmaPlus", ebins, emin, emax));
  fOutputContainer->Add(new TH1F("MC_AntiSigmaPlus_5_CPV2_Disp2", "Spectrum of AntiSigmaPlus", ebins, emin, emax));
  fOutputContainer->Add(new TH1F("MC_AntiSigmaMinus_5_CPV2", "Spectrum of AntiSigmaMinus", ebins, emin, emax));
  fOutputContainer->Add(new TH1F("MC_AntiSigmaMinus_5_Disp2", "Spectrum of AntiSigmaMinus", ebins, emin, emax));
  fOutputContainer->Add(new TH1F("MC_AntiSigmaMinus_5_CPV2_Disp2", "Spectrum of AntiSigmaMinus", ebins, emin, emax));

  // Spectra of identified clusters with CPV и Disp cuts
  for (int iPID = 0; iPID < 4; iPID++) {
    fOutputContainer->Add(
      new TH1F(Form("%sSpectrum_Photon", cPID[iPID]), "Spectrum of clusters from photon", ebins, emin, emax));
    fOutputContainer->Add(
      new TH1F(Form("%sSpectrum_Electron", cPID[iPID]), "Spectrum of clusters from electron", ebins, emin, emax));
    fOutputContainer->Add(
      new TH1F(Form("%sSpectrum_Positron", cPID[iPID]), "Spectrum of clusters from positron", ebins, emin, emax));
    fOutputContainer->Add(
      new TH1F(Form("%sSpectrum_Proton", cPID[iPID]), "Spectrum of clusters from proton", ebins, emin, emax));
    fOutputContainer->Add(
      new TH1F(Form("%sSpectrum_AntiProton", cPID[iPID]), "Spectrum of clusters from antiproton", ebins, emin, emax));
    fOutputContainer->Add(
      new TH1F(Form("%sSpectrum_PiPlus", cPID[iPID]), "Spectrum of clusters from pi plus", ebins, emin, emax));
    fOutputContainer->Add(
      new TH1F(Form("%sSpectrum_PiMinus", cPID[iPID]), "Spectrum of clusters from pi minus", ebins, emin, emax));
    fOutputContainer->Add(
      new TH1F(Form("%sSpectrum_Neutron", cPID[iPID]), "Spectrum of clusters from neutron", ebins, emin, emax));
    fOutputContainer->Add(
      new TH1F(Form("%sSpectrum_AntiNeutron", cPID[iPID]), "Spectrum of clusters from antineutron", ebins, emin, emax));
    fOutputContainer->Add(
      new TH1F(Form("%sSpectrum_KPlus", cPID[iPID]), "Spectrum of clusters from K plus", ebins, emin, emax));
    fOutputContainer->Add(
      new TH1F(Form("%sSpectrum_KMinus", cPID[iPID]), "Spectrum of clusters from K minus", ebins, emin, emax));
    fOutputContainer->Add(
      new TH1F(Form("%sSpectrum_KLong", cPID[iPID]), "Spectrum of clusters from K long", ebins, emin, emax));
    fOutputContainer->Add(
      new TH1F(Form("%sSpectrum_Other", cPID[iPID]), "Spectrum of clusters from other", ebins, emin, emax));
    fOutputContainer->Add(new TH1F(Form("%sSpectrum_Sum", cPID[iPID]), "Spectrum of clusters sum", ebins, emin, emax));
  }

  // Time of flight vs energy
  fOutputContainer->Add(
    new TH2F("hClusterTOFvsE", "Cluster time vs energy", tofbins, tofmin, tofmax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("AntiNeutronClusterTOFvsE", "Antineutron cluster time vs energy", tofbins, tofmin,
                                 tofmax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("AntiNeutronClusterTOFvsMCE", "Antineutron cluster time vs MC energy", tofbins, tofmin,
                                 tofmax, ebins, emin, emax));

  // Reconstruction of momentum; dependences
  fOutputContainer->Add(new TH2F("AntiNeutronPolar_clustervsMC", "Cluster vs MC polar angle", 100, 1, 2., 400, 0, 4));
  fOutputContainer->Add(
    new TH2F("AntiNeutronPolar_clustervscluster-MC", "Cluster vs Cluster-MC polar angle", 100, 1, 2., 400, -2, 2.));
  fOutputContainer->Add(
    new TH2F("AntiNeutronAzimuth_clustervsMC", "Cluster vs MC azimuth angle", 100, -2.5, -0.5, 200, -3, 1.));
  fOutputContainer->Add(new TH2F("AntiNeutronAzimuth_clustervscluster-MC", "Cluster vs Cluster-MC azimuth angle", 100,
                                 -2.5, -0.5, 300, -3, 3.));
  fOutputContainer->Add(
    new TH2F("AntiNeutronMomentum_reconvsMC", "Reconstructed vs MC momentum", ebins, emin, emax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("AntiNeutronMomentum_MCvsMC-recon_Plus",
                                 "Reconstructed vs Reconstructed momentum - MC", ebins, emin, emax, recbins, recmin,
                                 recmax));
  fOutputContainer->Add(new TH2F("AntiNeutronMomentum_MCvsMC-recon_Minus",
                                 "Reconstructed vs Reconstructed momentum - MC", ebins, emin, emax, recbins, recmin,
                                 recmax));

  for (int istep = 0; istep < 7; istep++) {
    fOutputContainer->Add(new TH2F(Form("AntiNeutronMomentum_MCvsRecon%d-MC_Im", 5 * istep),
                                   "Reconstructed vs Reconstructed momentum - MC imaginary", ebins, emin, emax, recbins,
                                   recmin, recmax));
    fOutputContainer->Add(new TH2F(Form("AntiNeutronMomentum_MCvsRecon%d-MC", 5 * istep),
                                   "Reconstructed vs Reconstructed momentum - MC", ebins, emin, emax, recbins, recmin,
                                   recmax));
  }

  fOutputContainer->Add(
    new TH2F("AntiNeutronMomentum_clustervsMC", "Cluster vs MC momentum", ebins, emin, emax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("AntiNeutronMomentum_MCvscluster-MC", "MC vs Cluster-MC momentum", ebins, emin, emax,
                                 recbins, recmin, recmax));

  // Reconstructed and unrec. momentum, TOF dependences
  fOutputContainer->Add(
    new TH2F("NegativeTime+", "TOF vs Energy of cluster if real", tofbins, tofmin, tofmax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("NegativeTime-", "TOF vs Energy of cluster if complex", tofbins, tofmin, tofmax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("AntineutronNegativeTime+", "TOF vs Energy of antineutron cluster if real", tofbins,
                                 tofmin, tofmax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("AntineutronNegativeTime-", "TOF vs Energy of antineutron cluster if complex", tofbins,
                                 tofmin, tofmax, ebins, emin, emax));

  // Topological dependences
  fOutputContainer->Add(
    new TH2F("DCA_Daughters_Plus", "DCA between daughters", dcabins, dcamin, dcamax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("DCA_Daughters_Minus", "DCA between daughters", dcabins, dcamin, dcamax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("DCA_Daughters_AntiSigmaPlus", "DCA between daughters", dcabins, dcamin, dcamax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("DCA_Daughters_AntiSigmaMinus", "DCA between daughters", dcabins, dcamin, dcamax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("CPA_Daughters_Plus", "DCA between daughters", cpabins, cpamin, cpamax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("CPA_Daughters_Minus", "DCA between daughters", cpabins, cpamin, cpamax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("CPA_Daughters_AntiSigmaPlus", "DCA between daughters", cpabins, cpamin, cpamax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("CPA_Daughters_AntiSigmaMinus", "DCA between daughters", cpabins, cpamin, cpamax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("RAD_Daughters_Plus", "RAD between daughters", radbins, radmin, radmax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("RAD_Daughters_Minus", "RAD between daughters", radbins, radmin, radmax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("RAD_Daughters_AntiSigmaPlus", "RAD between daughters", radbins, radmin, radmax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("RAD_Daughters_AntiSigmaMinus", "RAD between daughters", radbins, radmin, radmax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("DCA_XY_Plus", "DCA XY vs pT of track", dcabins, dcaxymin, dcaxymax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("DCA_XY_Minus", "DCA XY vs pT of track", dcabins, dcaxymin, dcaxymax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("DCA_XY_AntiSigmaPlus", "DCA XY vs pT of pion from AntiSigma", dcabins, dcaxymin,
                                 dcaxymax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("DCA_XY_AntiSigmaMinus", "DCA XY vs pT of pion from AntiSigma", dcabins, dcaxymin,
                                 dcaxymax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("DCA_Z_Plus", "DCA Z vs pT of track", dcabins, dcazmin, dcazmax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("DCA_Z_Minus", "DCA Z vs pT of track", dcabins, dcazmin, dcazmax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("DCA_Z_AntiSigmaPlus", "DCA Z vs pT of pion from AntiSigma", dcabins, dcazmin, dcazmax,
                                 ebins, emin, emax));
  fOutputContainer->Add(new TH2F("DCA_Z_AntiSigmaMinus", "DCA Z vs pT of pion from AntiSigma", dcabins, dcazmin,
                                 dcazmax, ebins, emin, emax));

  // Invariant mass
  for (int iIM = 0; iIM < 5; iIM++) {
    for (int iMC = 0; iMC < 2; iMC++) {
      for (int iCH = 0; iCH < 6; iCH++) {
        fOutputContainer->Add(new TH2F(Form("%s%sInvMass%s", cIM[iIM], cMC[iMC], cCH[iCH]), "Invariant mass", invbins,
                                       invmin, invmax, ebins, emin, emax));
      }
    }
  }

  // Mixed invariant mass
  for (int iIM = 0; iIM < 5; iIM++) {
    for (int iCH = 0; iCH < 2; iCH++) {
      fOutputContainer->Add(new TH2F(Form("%sMixed_InvMass%s", cIM[iIM], cCH[iCH]), "Invariant mass", invbins, invmin,
                                     invmax, ebins, emin, emax));
    }
  }

  // Invariant mass
  fOutputContainer->Add(
    new TH2F("InvMass_AntiNeutron_Charge1", "Invariant mass", invbins, invmin, invmax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("InvMass_AntiNeutron_Charge-1", "Invariant mass", invbins, invmin, invmax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("InvMass_PiMinus", "Invariant mass", invbins, invmin, invmax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("InvMass_PiPlus", "Invariant mass", invbins, invmin, invmax, ebins, emin, emax));
  //МС invariant mass
  fOutputContainer->Add(
    new TH2F("MC_InvMass_AntiNeutron_Charge1", "Invariant mass", invbins, invmin, invmax, ebins, emin, emax));
  fOutputContainer->Add(
    new TH2F("MC_InvMass_AntiNeutron_Charge-1", "Invariant mass", invbins, invmin, invmax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("MC_InvMass_PiMinus", "Invariant mass", invbins, invmin, invmax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("MC_InvMass_PiPlus", "Invariant mass", invbins, invmin, invmax, ebins, emin, emax));

  //========================================//
  fOutputContainer->Add(new TH1F("hSelEvents", "Events selected", 10, 0., 10.));
  fOutputContainer->Add(new TH1F("hCellMultEvent", "PHOS cell multiplicity per event", 2000, 0, 2000));
  fOutputContainer->Add(new TH1F("hClusterMult", "CaloCluster multiplicity", 100, 0, 100));
  fOutputContainer->Add(new TH1F("hPHOSClusterMult", "PHOS cluster multiplicity", 100, 0, 100));
  fOutputContainer->Add(new TH1F("hCellEnergy", "Cell energy", 5000, 0., 50.));
  fOutputContainer->Add(new TH1F("hClusterEnergy", "Cluster energy", 5000, 0., 50.));
  fOutputContainer->Add(new TH2F("hClusterEvsN", "Cluster energy vs digit multiplicity", 5000, 0., 50., 40, 0., 40.));
  fOutputContainer->Add(new TH1F("hCellMultClu", "Cell multiplicity per cluster", 200, 0, 200));
  fOutputContainer->Add(new TH1F("hModule", "Module events", 5, 0., 5.));
  for (Int_t module = 1; module < 5; module++) {
    fOutputContainer->Add(
      new TH2F(Form("hModule%d", module), Form("Cluster occupancy in module %d", module), 64, 0., 64, 56, 0., 56.));
  }
  fOutputContainer->Add(new TH1F("hZvertex", "Z vertex", 200, -50., +50.));
  fOutputContainer->Add(new TH1F("hNvertexTracks", "N of primary tracks from the primary vertex", 150, 0., 150.));
  fOutputContainer->Add(new TH1F("hT0TOF", "T0 time (s)", 2000, -2.0e-9, 2.0e-9));
  fOutputContainer->Add(new TH1F("hTrackMult", "Charged track multiplicity", 150, 0., 150.));

  for (Int_t i = 0; i < 10; i++)
    for (Int_t j = 0; j < fNCenBin; j++)
      fPHOSEvents[i][j] = 0x0; // Container for PHOS photons

  PostData(1, fOutputContainer);
}

void AliAnalysisSigmaBarCharged::UserExec(Option_t*)
{
  // analyze one event
  FillHistogram("hSelEvents", 1);
  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fEvent) {
    Printf("ERROR: Could not retrieve event");
    PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents", 2);

  // Test vertex
  const AliAODVertex* esdVertex5 = fEvent->GetPrimaryVertex();
  fvtx5[0] = esdVertex5->GetX();
  fvtx5[1] = esdVertex5->GetY();
  fvtx5[2] = esdVertex5->GetZ();

  FillHistogram("hNvertexTracks", esdVertex5->GetNContributors());
  FillHistogram("hZvertex", fvtx5[2]);
  if (TMath::Abs(fvtx5[2]) > 10.) {
    PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents", 4);
  Int_t zvtx = (Int_t)((fvtx5[2] + 10.) / 2.);
  if (zvtx < 0)
    zvtx = 0;
  if (zvtx > 9)
    zvtx = 9;

  // Pileup selection

  FillHistogram("hSelEvents", 5);

  if (!fPHOSEvents[zvtx][fCentBin])
    fPHOSEvents[zvtx][fCentBin] = new TList();
  fCurrentMixedList = fPHOSEvents[zvtx][fCentBin];

  // particle identification
  if (!fPIDResponse) {
    AliAODInputHandler* aodH =
      dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (aodH) {
      fPIDResponse = aodH->GetPIDResponse();
    } else {
      Printf("ERROR: Could not get AODInputHandler");
      return;
    }
  }

  SelectSigma();

  // Remove old events
  fCurrentMixedList->AddFirst(fGamma);
  fGamma = 0x0;
  if (fCurrentMixedList->GetSize() > 100) {
    TClonesArray* tmp = static_cast<TClonesArray*>(fCurrentMixedList->Last());
    fCurrentMixedList->Remove(tmp);
    delete tmp;
  }

  PostData(1, fOutputContainer);
}

void AliAnalysisSigmaBarCharged::SelectSigma()
{
  Int_t inPHOS = 0;
  // List of antineutron clusters
  if (fGamma)
    fGamma->Clear();
  else
    fGamma = new TClonesArray("AliCaloPhoton", 100);

  // Number of tracks and clusters
  Int_t multClust = fEvent->GetNumberOfCaloClusters();
  Int_t multTracks = fEvent->GetNumberOfTracks();

  if (fIsMC) {
    // MC stack
    fStack = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());

    // Generated particles
    Int_t multMC = fStack->GetEntriesFast();
    for (Int_t i = 0; i < multMC; i++) {
      AliAODMCParticle* prim = (AliAODMCParticle*)fStack->At(i);
      if (prim->GetPdgCode() == -3222 && TMath::Abs(prim->Y()) < 1) {
        FillHistogram("MC_AntiSigmaMinus_0", prim->Pt());
      }
      if (prim->GetPdgCode() == -3112 && TMath::Abs(prim->Y()) < 1) {
        FillHistogram("MC_AntiSigmaPlus_0", prim->Pt());
      }
      if (prim->GetPdgCode() == -2112) {
        FillHistogram("MC_A_0", prim->Pt());
        Int_t iparent = prim->GetMother();
        if (iparent > -1) {
          AliAODMCParticle* parent = (AliAODMCParticle*)fStack->At(iparent);
          if (parent->GetPdgCode() == -3222 && TMath::Abs(prim->Y()) < 1) {
            FillHistogram("MC_AntiSigmaMinus_1", parent->Pt());
            // cout << "an- label: " << prim->GetLabel() << "parent label: " << parent->GetLabel() << endl;
          }
          if (parent->GetPdgCode() == -3112 && TMath::Abs(prim->Y()) < 1) {
            FillHistogram("MC_AntiSigmaPlus_1", parent->Pt());
          }
        }
      }
      if (prim->GetPdgCode() == 211 && TMath::Abs(prim->Eta()) < 0.8) {
        FillHistogram("MC_PionPlus_0", prim->Pt());
      }
      if (prim->GetPdgCode() == -211 && TMath::Abs(prim->Eta()) < 0.8) {
        FillHistogram("MC_PionMinus_0", prim->Pt());
      }
    }
  } // End of MC particles loop

  // Scan clusters
  for (Int_t i = 0; i < multClust; i++) {
    AliAODCaloCluster* clu = fEvent->GetCaloCluster(i);
    if (clu->GetType() != AliVCluster::kPHOSNeutral)
      continue;
    Double_t cluE = clu->E();
    if (cluE < fCluMinECut)
      continue;
    FillHistogram("All_Disp_MinCut", clu->GetM20(), clu->GetM02(), cluE);
    FillHistogram("Spectrum_Cluster_1", cluE);
    // Dispersion and spectrum of antineutron clusters
    if (fIsMC) {
      Int_t primLabel0 = clu->GetLabelAt(0);
      if (primLabel0 > -1) {
        AliAODMCParticle* prim0 = (AliAODMCParticle*)fStack->At(primLabel0);
        FillHistogram("EnergyclustervsMC", cluE, prim0->E());
        Int_t iparent = prim0->GetMother();
        if (prim0->GetPdgCode() == -2112) {
          FillHistogram("AntiNeutron_Disp_MinCut", clu->GetM20(), clu->GetM02(), cluE);
          FillHistogram("AntiNeutron_EnergyclustervsMC", cluE, prim0->E());
          FillHistogram("MC_A_1", prim0->Pt());
          if (iparent > -1) {
            AliAODMCParticle* parent = (AliAODMCParticle*)fStack->At(iparent);
            if (parent->GetPdgCode() == -3222) {
              FillHistogram("MC_AntiSigmaMinus_2", parent->Pt());
            }
            if (parent->GetPdgCode() == -3112) {
              FillHistogram("MC_AntiSigmaPlus_2", parent->Pt());
            }
          }
        }
      }
    }
    // Exotic or single cell cluster
    if (clu->GetM02() < 0.2)
      continue;
    FillHistogram("Spectrum_Cluster_2", cluE);
    if (fIsMC) {
      Int_t primLabel0 = clu->GetLabelAt(0);
      if (primLabel0 > -1) {
        AliAODMCParticle* prim0 = (AliAODMCParticle*)fStack->At(primLabel0);
        Int_t iparent = prim0->GetMother();
        if (prim0->GetPdgCode() == -2112) {
          FillHistogram("MC_A_2", prim0->Pt());
          if (iparent > -1) {
            AliAODMCParticle* parent = (AliAODMCParticle*)fStack->At(iparent);
            if (parent->GetPdgCode() == -3222) {
              FillHistogram("MC_AntiSigmaMinus_3", parent->Pt());
            }
            if (parent->GetPdgCode() == -3112) {
              FillHistogram("MC_AntiSigmaPlus_3", parent->Pt());
            }
          }
        }
      }
    }
    // Time cut in data
    if (!fIsMC) {
      if (TMath::Abs(clu->GetTOF()) > fCluTimeCut) {
        continue;
      }
      FillHistogram("Spectrum_Cluster_TOF", cluE);
    }
    // Energy cut for antineutrons
    if (cluE < fCluNbarMinE) { // append realistic resolution
      continue;
    }

    FillHistogram("All_Disp", clu->GetM20(), clu->GetM02(), cluE);
    FillHistogram("Spectrum_Cluster_3", cluE);
    if (fIsMC) {
      Int_t primLabel0 = clu->GetLabelAt(0);
      if (primLabel0 > -1) {
        AliAODMCParticle* prim0 = (AliAODMCParticle*)fStack->At(primLabel0);
        Int_t iparent = prim0->GetMother();
        if (prim0->GetPdgCode() == -2112) {
          FillHistogram("AntiNeutron_Disp", clu->GetM20(), clu->GetM02(), cluE);
          FillHistogram("MC_A_3", prim0->Pt());
          if (iparent > -1) {
            AliAODMCParticle* parent = (AliAODMCParticle*)fStack->At(iparent);
            if (parent->GetPdgCode() == -3222) {
              FillHistogram("MC_AntiSigmaMinus_4", parent->Pt());
            }
            if (parent->GetPdgCode() == -3112) {
              FillHistogram("MC_AntiSigmaPlus_4", parent->Pt());
            }
          }
        }
      }
    }
    // Reconstructed (photon) momentum
    TLorentzVector lvclu;
    clu->GetMomentum(lvclu, fvtx5);
    // Position in PHOS surface
    Float_t pos[3];
    clu->GetPosition(pos);
    const Double_t c = 29979245800.; // speed of light in cm/sec
    const Double_t mbar = 0.939485;  // neutron mass
    Double_t t = clu->GetTOF();
    Double_t r = TMath::Sqrt((fvtx5[0] - pos[0]) * (fvtx5[0] - pos[0]) + (fvtx5[1] - pos[1]) * (fvtx5[1] - pos[1]) +
                             (fvtx5[2] - pos[2]) * (fvtx5[2] - pos[2]));
    Double_t tgamma = r / c;
    if (!fIsMC) { // Real data calibrated wrt photon arrival
      t = t + tgamma;
    }
    // Time resolution
    Double_t sigt = 0.;
    if (fIsMC && fPHOSClusterTOFOption == 1) { // append realistic resolution
      sigt = RealRes(cluE);
      t = gRandom->Gaus(t, sigt);
    }

    FillHistogram("hClusterTOFvsE", t, cluE);

    // Reconstruct nbar momentum
    if (fIsMC) {
      Int_t primLabel = clu->GetLabelAt(0);
      if (primLabel > -1) {
        AliAODMCParticle* prim = (AliAODMCParticle*)fStack->At(primLabel);
        if (prim->GetPdgCode() == -2112) {
          Double_t px = prim->Px();
          Double_t py = prim->Py();
          Double_t pz = prim->Pz();
          TVector3 vecmc(px, py, pz);

          // Vary depth of shower maximum: without shift
          for (int istep = 0; istep < 7; istep++) {
            Double_t precs = mbar * mbar / (pow(t * c / (r + 5. * istep), 2) - 1);
            FillHistogram(Form("AntiNeutronMomentum_MCvsRecon%d-MC_Im", 5 * istep), vecmc.Mag(),
                          -vecmc.Mag() * vecmc.Mag() + precs);
            if (pow(t * c / (r + 5. * istep), 2) - 1 > 0) {
              Double_t prec = mbar / sqrt(pow(t * c / (r + 5. * istep), 2) - 1);
              FillHistogram(Form("AntiNeutronMomentum_MCvsRecon%d-MC", 5 * istep), vecmc.Mag(), -vecmc.Mag() + prec);
            }
          }
        }
      }
    }
    // If momentum can be reconstructed at optimal depth
    if (pow(t * c / (r + fOptDepth), 2) - 1 > 0) {
      FillHistogram("NegativeTime+", t, cluE);
      FillHistogram("Spectrum_Cluster_4", cluE);
      Double_t prec = mbar / sqrt(pow(t * c / (r + fOptDepth), 2) - 1);
      // 4-momentum of nbar
      TVector3 recon;
      recon.SetMagThetaPhi(prec, lvclu.Theta(), lvclu.Phi());
      Double_t reconE = sqrt(mbar * mbar + recon.Mag() * recon.Mag());

      if (inPHOS >= fGamma->GetSize()) {
        fGamma->Expand(inPHOS + 50);
      }
      AliCaloPhoton* ph = new ((*fGamma)[inPHOS++]) AliCaloPhoton(recon.X(), recon.Y(), recon.Z(), reconE);
      // ph->SetTOFBit((clu->GetTOF()>-50.e-9) && (clu->GetTOF() <50.e-9) );

      // CPV cut
      ph->SetCPV2Bit(clu->GetEmcCpvDistance() > fCPVCut);

      // Anti-photon dispersion cut
      ph->SetDisp2Bit((clu->Chi2() > fDispCut * fDispCut) && (clu->GetM20() >= fDispA * clu->GetM02() + fDispB));

      Int_t primLabel = clu->GetLabelAt(0);
      ph->SetPrimaryAtVertex(primLabel);

      Float_t pt = ph->Pt();
      FillHistogram("Rec_Spectrum_Cluster_4", pt);

      // CPV2 cut
      if (ph->IsCPV2OK()) {
        FillHistogram("Spectrum_Cluster_5", cluE);
        FillHistogram("Rec_Spectrum_Cluster_5", pt);
      }
      // Disp2 cut
      if (ph->IsDisp2OK()) {
        FillHistogram("Spectrum_Cluster_6", cluE);
        FillHistogram("Rec_Spectrum_Cluster_6", pt);
        FillHistogram("Disp2_Cut", clu->GetM20(), clu->GetM02(), cluE);
        if (fIsMC) {
          if (primLabel > -1) {
            AliAODMCParticle* prim0 = (AliAODMCParticle*)fStack->At(primLabel);
            if (prim0->GetPdgCode() == -2112) {
              FillHistogram("AntiNeutron_Disp2_Cut", clu->GetM20(), clu->GetM02(), cluE);
            }
          }
        }
        // CPV2 and Disp2 cut
        if (ph->IsCPV2OK()) {
          FillHistogram("Spectrum_Cluster_7", cluE);
          FillHistogram("Rec_Spectrum_Cluster_7", pt);
        }
      }

      if (fIsMC) {
        if (primLabel > -1) {
          AliAODMCParticle* prim = (AliAODMCParticle*)fStack->At(primLabel);
          Int_t iparent = prim->GetMother();
          TVector3 vecmc(prim->Px(), prim->Py(), prim->Pz());
          Double_t primE = prim->E();
          Double_t primPt = prim->Pt();
          Double_t primP = vecmc.Mag();

          // Classify particle
          TString partName;
          switch (prim->GetPdgCode()) {
            case 22:
              partName = "Photon";
              break;
            case 11:
              partName = "Electron";
              break;
            case -11:
              partName = "Positron";
              break;
            case 2212:
              partName = "Proton";
              break;
            case -2212:
              partName = "AntiProton";
              break;
            case 211:
              partName = "PiPlus";
              break;
            case -211:
              partName = "PiMinus";
              break;
            case 112:
              partName = "Neutron";
              break;
            case -2112:
              partName = "AntiNeutron";
              break;
            case 321:
              partName = "KPlus";
              break;
            case -321:
              partName = "KMinus";
              break;
            case 130:
              partName = "KLong";
              break;
            default:
              partName = "Other";
              FillHistogram("PDG_Other", prim->GetPdgCode());
          }

          // Spectrum of particles in PHOS
          FillHistogram("Spectrum_Sum", primE);
          FillHistogram(Form("Spectrum_%s", partName.Data()), primE);
          // CPV2 cut
          if (ph->IsCPV2OK()) {
            FillHistogram("CPV2_Spectrum_Sum", primE);
            FillHistogram(Form("CPV2_Spectrum_%s", partName.Data()), primE);
          }
          // Disp2 cut
          if (ph->IsDisp2OK()) {
            FillHistogram("Disp2_Spectrum_Sum", primE);
            FillHistogram(Form("Disp2_Spectrum_%s", partName.Data()), primE);
            if (ph->IsCPV2OK()) {
              FillHistogram("CPV2_Disp2_Spectrum_Sum", primE);
              FillHistogram(Form("CPV2_Disp2_Spectrum_%s", partName.Data()), primE);
            }
          }
          // For antineutrons
          if (prim->GetPdgCode() == -2112) {
            // Difference rec and true momenta
            // FillHistogram("AntiNeutronMomentum_MCvsRecon10-MC", primP, prec - primP);
            // spectra
            FillHistogram("MC_A_4", primPt);
            FillHistogram("Rec_A_4", pt);
            // Angles rec vs true
            FillHistogram("AntiNeutronPolar_clustervsMC", lvclu.Theta(), vecmc.Theta());
            FillHistogram("AntiNeutronPolar_clustervscluster-MC", lvclu.Theta(), lvclu.Theta() - vecmc.Theta());
            FillHistogram("AntiNeutronAzimuth_clustervsMC", lvclu.Phi(), vecmc.Phi());
            FillHistogram("AntiNeutronAzimuth_clustervscluster-MC", lvclu.Phi(), lvclu.Phi() - vecmc.Phi());
            // Momemta rec vs true
            FillHistogram("AntiNeutronMomentum_reconvsMC", recon.Mag(), primP);
            FillHistogram("AntiNeutronMomentum_clustervsMC", lvclu.Mag(), primP);
            FillHistogram("AntiNeutronMomentum_MCvscluster-MC", primP, lvclu.Mag() - primP);
            // Time dependences
            FillHistogram("AntiNeutronClusterTOFvsE", t, cluE);
            FillHistogram("AntiNeutronClusterTOFvsMCE", t, primE);
            //  Antineutrons
            FillHistogram("AntineutronNegativeTime+", t, cluE);
            if (ph->IsCPV2OK()) {
              FillHistogram("MC_A_5", primPt);
              FillHistogram("Rec_A_5", pt);
            }
            if (ph->IsDisp2OK()) {
              FillHistogram("MC_A_6", primPt);
              FillHistogram("Rec_A_6", pt);
              if (ph->IsCPV2OK()) {
                FillHistogram("MC_A_7", primPt);
                FillHistogram("Rec_A_7", pt);
              }
            }
            if (iparent > -1) {
              AliAODMCParticle* parent = (AliAODMCParticle*)fStack->At(iparent);
              Double_t parentPt = parent->Pt();
              TVector3 vecpar(parent->Px(), parent->Py(), parent->Pz());

              if (parent->GetPdgCode() == -3222) {
                FillHistogram("AntiNeutronMomentum_MCvsMC-recon_Minus", primP, -primP + recon.Mag());
                FillHistogram("MC_AntiSigmaMinus_5", parentPt);
                // CPV2 cut
                if (ph->IsCPV2OK()) {
                  FillHistogram("MC_AntiSigmaMinus_5_CPV2", parentPt);
                }
                // Disp2 cut
                if (ph->IsDisp2OK()) {
                  FillHistogram("MC_AntiSigmaMinus_5_Disp2", parentPt);
                  // CPV2 and Disp2 cut
                  if (ph->IsCPV2OK()) {
                    FillHistogram("MC_AntiSigmaMinus_5_CPV2_Disp2", parentPt);
                  }
                }
              }
              if (parent->GetPdgCode() == -3112) {
                FillHistogram("AntiNeutronMomentum_MCvsMC-recon_Plus", primP, -primP + recon.Mag());
                FillHistogram("MC_AntiSigmaPlus_5", parentPt);
                // CPV2 cut
                if (ph->IsCPV2OK()) {
                  FillHistogram("MC_AntiSigmaPlus_5_CPV2", parentPt);
                }
                // Disp2 cut
                if (ph->IsDisp2OK()) {
                  FillHistogram("MC_AntiSigmaPlus_5_Disp2", parentPt);
                  // CPV2 and Disp2 cut
                  if (ph->IsCPV2OK()) {
                    FillHistogram("MC_AntiSigmaPlus_5_CPV2_Disp2", parentPt);
                  }
                }
              }
            }
          }
        }
      }
    } else {
      // Clusters with imaginary reconstructed momentum
      FillHistogram("NegativeTime-", t, cluE);
      if (fIsMC) {
        Int_t primLabel = clu->GetLabelAt(0);
        if (primLabel > -1) {
          AliAODMCParticle* prim0 = (AliAODMCParticle*)fStack->At(primLabel);
          if (prim0->GetPdgCode() == -2112) {
            FillHistogram("AntineutronNegativeTime-", t, cluE);
          }
        }
      }
    }
  }

  // Tracks
  TIter nextEv(fCurrentMixedList);
  for (Int_t i = 0; i < multTracks; i++) {
    AliAODTrack* track = (AliAODTrack*)fEvent->GetTrack(i);
    if (track->GetFilterMap() != fTracksBits)
      continue;
    if (TMath::Abs(track->Eta()) > 0.8)
      continue;
    Double_t pxtr = track->Px();
    Double_t pytr = track->Py();
    Double_t pztr = track->Pz();
    const Double_t mpi = 0.13957039;
    Double_t etr = sqrt(mpi * mpi + (pxtr * pxtr + pytr * pytr + pztr * pztr));
    TLorentzVector lvtr(pxtr, pytr, pztr, etr);
    Double_t trackPt = track->Pt();
    Int_t primLabelTrack = track->GetLabel();
    // DCA
    Float_t b[2];
    track->GetImpactParameters(b[0], b[1]); // xy:-5.43674e-02+3.14319e-02/x z:-3.65702e-02e-02+5.38728e-02/x pt<0.5
    if (b[0] == -999 || b[1] == -999)
      continue;
    // Track PID in TPC sigmas
    AliVParticle* inEvHMain = dynamic_cast<AliVParticle*>(track);
    Bool_t pidPion = kFALSE, pidKaon = kFALSE, pidProton = kFALSE;
    Double_t nsigmaProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton));
    Double_t nsigmaKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon));
    Double_t nsigmaPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion));
    // guess the particle based on the smaller nsigma
    if ((nsigmaKaon < nsigmaPion) && (nsigmaKaon < nsigmaProton) && (nsigmaKaon < 3))
      pidKaon = kTRUE;
    if ((nsigmaPion < nsigmaKaon) && (nsigmaPion < nsigmaProton) && (nsigmaPion < 3))
      pidPion = kTRUE;
    if ((nsigmaProton < nsigmaKaon) && (nsigmaProton < nsigmaPion) && (nsigmaProton < 3))
      pidProton = kTRUE;
//     if (!pidPion && !pidKaon && !pidProton)
//       pidUndef = kTRUE;
    // kink daughters
    if (track->GetKinkIndex(0) != 0) {
      std::cout << "label: " << primLabelTrack << std::endl;
      if (primLabelTrack > -1) {
        AliAODMCParticle* primpion = (AliAODMCParticle*)fStack->At(primLabelTrack);
        std::cout << "pdg: " << primpion->GetPdgCode() << std::endl;
      }
      continue; // TODO?????
    }
    if (track->Charge() == 1) {
      FillHistogram("TrackPlus_1", trackPt);
      FillHistogram("DCA_XY_Plus", b[0], trackPt);
      FillHistogram("DCA_Z_Plus", b[1], trackPt);
    }
    if (track->Charge() == -1) {
      FillHistogram("TrackMinus_1", trackPt);
      FillHistogram("DCA_XY_Minus", b[0], trackPt);
      FillHistogram("DCA_Z_Minus", b[1], trackPt);
    }
    // Pions
    if (fIsMC) {
      if (primLabelTrack > -1) {
        AliAODMCParticle* primpion = (AliAODMCParticle*)fStack->At(primLabelTrack);
        if (primpion->GetPdgCode() == 211) {
          FillHistogram("PionPlus_1", trackPt);
        }
        if (primpion->GetPdgCode() == -211) {
          FillHistogram("PionMinus_1", trackPt);
        }
        if (primpion->GetMother() > -1) {
          AliAODMCParticle* pionparent = (AliAODMCParticle*)fStack->At(primpion->GetMother());
          if (pionparent->GetPdgCode() == -3112 && primpion->GetPdgCode() == 211) {
            FillHistogram("MC_AntiSigmaPlus_6", pionparent->Pt());
            FillHistogram("DCA_XY_AntiSigmaPlus", b[0], trackPt);
            FillHistogram("DCA_Z_AntiSigmaPlus", b[1], trackPt);
          }
          if (pionparent->GetPdgCode() == -3222 && primpion->GetPdgCode() == -211) {
            FillHistogram("MC_AntiSigmaMinus_6", pionparent->Pt());
            FillHistogram("DCA_XY_AntiSigmaMinus", b[0], trackPt);
            FillHistogram("DCA_Z_AntiSigmaMinus", b[1], trackPt);
          }
        }
      }
    }

    if (pidPion) {
      if (track->Charge() == 1) {
        FillHistogram("TrackPlus_2", trackPt);
      }
      if (track->Charge() == -1) {
        FillHistogram("TrackMinus_2", trackPt);
      }
      if (fIsMC) {
        if (primLabelTrack > -1) {
          AliAODMCParticle* primpion = (AliAODMCParticle*)fStack->At(primLabelTrack);
          if (primpion->GetPdgCode() == 211) {
            FillHistogram("PionPlus_2", trackPt);
          }
          if (primpion->GetPdgCode() == -211) {
            FillHistogram("PionMinus_2", trackPt);
          }
          if (primpion->GetMother() > -1) {
            AliAODMCParticle* pionparent = (AliAODMCParticle*)fStack->At(primpion->GetMother());
            if (pionparent->GetPdgCode() == -3112 && primpion->GetPdgCode() == 211) {
              // FillHistogram("DCA_XYvsPionPt_Parent_AntiSigmaPlus", b[0], trackPt);
              // FillHistogram("DCA_ZvsPionPt_Parent_AntiSigmaPlus", b[1], trackPt);
              FillHistogram("MC_AntiSigmaPlus_7", pionparent->Pt());
            }
            if (pionparent->GetPdgCode() == -3222 && primpion->GetPdgCode() == -211) {
              // FillHistogram("DCA_XYvsPionPt_Parent_AntiSigmaMinus", b[0], trackPt);
              // FillHistogram("DCA_ZvsPionPt_Parent_AntiSigmaMinus", b[1], trackPt);
              FillHistogram("MC_AntiSigmaMinus_7", pionparent->Pt());
            }
          }
        }
      }

      //=================Event mixing==============
      while (TClonesArray* event2 = static_cast<TClonesArray*>(nextEv())) {
        Int_t nPhotons2 = event2->GetEntriesFast();
        for (Int_t j = 0; j < nPhotons2; j++) {
          AliCaloPhoton* ph2 = static_cast<AliCaloPhoton*>(event2->At(j));
          //             Int_t primLabel2 = ph2->GetPrimaryAtVertex();
          Double_t m = (lvtr + *ph2).M();
          Double_t pt = (lvtr + *ph2).Pt();

          // Try to construct particle out of track and nbar
          Double_t bm = fEvent->GetMagneticField();

          AliESDtrack btrk((AliVTrack*)track);
          AliExternalTrackParam bt(btrk);

          Double_t cv[21];
          for (Int_t ii = 0; ii < 21; ii++)
            cv[ii] = 0.0;
          Double_t p[3] = { ph2->Px(), ph2->Py(), ph2->Pz() };
          AliExternalTrackParam lV0TrajObject(fvtx5, p, cv, +1), *nt = &lV0TrajObject;
          nt->ResetCovariance(1); // won't use

          Double_t xn, xp, dca;
          dca = nt->GetDCA(&bt, bm, xn, xp); //+: pt<0.25 dca>0.2, -:pt<0.25, dca>0.1
          nt->PropagateTo(xn, bm);
          bt.PropagateTo(xp, bm);
          AliESDv0 v0(*nt, nPhotons2, bt, multTracks);
          Float_t cpa = v0.GetV0CosineOfPointingAngle(fvtx5[0], fvtx5[1], fvtx5[2]); //>0.9 по модулю
          Double_t rad =
            sqrt((fvtx5[0] - v0.Xv()) * (fvtx5[0] - v0.Xv()) + (fvtx5[1] - v0.Yv()) * (fvtx5[1] - v0.Yv()) +
                 (fvtx5[2] - v0.Zv()) * (fvtx5[2] - v0.Zv())); //>0.2для +, 0.1

          FillHistogram(Form("Mixed_InvMass_Charge%d", track->Charge()), m, pt);
          if (ph2->IsCPV2OK()) {
            FillHistogram(Form("CPV2_Mixed_InvMass_Charge%d", track->Charge()), m, pt);
          }
          if (ph2->IsDisp2OK()) {
            FillHistogram(Form("Disp2_Mixed_InvMass_Charge%d", track->Charge()), m, pt);
            if (ph2->IsCPV2OK()) {
              FillHistogram(Form("All3_Mixed_InvMass_Charge%d", track->Charge()), m, pt);

              if ((TestDCAZ(trackPt, b[1])) && (TestDCAXY(trackPt, b[0])) && abs(cpa) > 0.9) {
                if (track->Charge() == 1 && rad > 0.2 && DCADaugPlus(trackPt, dca)) {
                  FillHistogram("All4_Mixed_InvMass_Charge1", m, pt);
                }
                if (track->Charge() == -1 && rad > 0.1 && DCADaugMinus(trackPt, dca)) {
                  FillHistogram("All4_Mixed_InvMass_Charge-1", m, pt);
                }
              }
            }
          }
        }
      }

      //======================Fill Real==========================
      for (Int_t j = 0; j < inPHOS; j++) {
        AliCaloPhoton* ph2 = static_cast<AliCaloPhoton*>(fGamma->At(j));
        Int_t primLabel2 = ph2->GetPrimaryAtVertex();
        Double_t m = (lvtr + *ph2).M();
        Double_t pt = (lvtr + *ph2).Pt();

        // Try to construct particle out of track and nbar
        Double_t bm = fEvent->GetMagneticField();

        AliESDtrack btrk((AliVTrack*)track);
        AliExternalTrackParam bt(btrk);

        Double_t cv[21];
        for (Int_t ii = 0; ii < 21; ii++)
          cv[ii] = 0.0;
        Double_t p[3] = { ph2->Px(), ph2->Py(), ph2->Pz() };
        AliExternalTrackParam lV0TrajObject(fvtx5, p, cv, +1), *nt = &lV0TrajObject;
        nt->ResetCovariance(1); // won't use

        Double_t xn, xp, dca;
        dca = nt->GetDCA(&bt, bm, xn, xp); //+: pt<0.25 dca>0.2, -:pt<0.25, dca>0.1
        nt->PropagateTo(xn, bm);
        bt.PropagateTo(xp, bm);
        AliESDv0 v0(*nt, inPHOS, bt, multTracks);
        Float_t cpa = v0.GetV0CosineOfPointingAngle(fvtx5[0], fvtx5[1], fvtx5[2]); //>0.9
        Double_t rad = sqrt((fvtx5[0] - v0.Xv()) * (fvtx5[0] - v0.Xv()) + (fvtx5[1] - v0.Yv()) * (fvtx5[1] - v0.Yv()) +
                            (fvtx5[2] - v0.Zv()) * (fvtx5[2] - v0.Zv())); //>0.2 for positive +, 0.1 for negative

        if (fIsMC) {
          if (primLabelTrack > -1) {
            AliAODMCParticle* primpion = (AliAODMCParticle*)fStack->At(primLabelTrack);
            if (primpion->GetMother() > -1) {
              AliAODMCParticle* pionparent = (AliAODMCParticle*)fStack->At(primpion->GetMother());
              if (pionparent->GetPdgCode() == -3112 && primpion->GetPdgCode() == 211) {
                FillHistogram("DCA_Daughters_AntiSigmaPlus", dca, trackPt);
                FillHistogram("CPA_Daughters_AntiSigmaPlus", cpa, trackPt);
                FillHistogram("RAD_Daughters_AntiSigmaPlus", rad, trackPt);
              }
              if (pionparent->GetPdgCode() == -3222 && primpion->GetPdgCode() == -211) {
                FillHistogram("DCA_Daughters_AntiSigmaMinus", dca, trackPt);
                FillHistogram("CPA_Daughters_AntiSigmaMinus", cpa, trackPt);
                FillHistogram("RAD_Daughters_AntiSigmaMinus", rad, trackPt);
              }
            }
          }
        }

        if (track->Charge() == 1) {
          FillHistogram("DCA_Daughters_Plus", dca, trackPt);
          FillHistogram("CPA_Daughters_Plus", cpa, trackPt);
          FillHistogram("RAD_Daughters_Plus", rad, trackPt);
        }
        if (track->Charge() == -1) {
          FillHistogram("DCA_Daughters_Minus", dca, trackPt);
          FillHistogram("CPA_Daughters_Minus", cpa, trackPt);
          FillHistogram("RAD_Daughters_Minus", rad, trackPt);
        }

        FillHistogram(Form("InvMass_Charge%d", track->Charge()), m, pt);
        if (ph2->IsCPV2OK() == kTRUE) {
          FillHistogram(Form("CPV2_InvMass_Charge%d", track->Charge()), m, pt);
        }
        if (ph2->IsDisp2OK()) {
          FillHistogram(Form("Disp2_InvMass_Charge%d", track->Charge()), m, pt);
          if (ph2->IsCPV2OK()) {
            FillHistogram(Form("All3_InvMass_Charge%d", track->Charge()), m, pt);

            if ((TestDCAZ(trackPt, b[1])) && (TestDCAXY(trackPt, b[0])) && abs(cpa) > 0.9) {
              if (track->Charge() == 1 && rad > 0.2 && DCADaugPlus(trackPt, dca)) {
                FillHistogram("All4_InvMass_Charge1", m, pt);
              }
              if (track->Charge() == -1 && rad > 0.1 && DCADaugMinus(trackPt, dca)) {
                FillHistogram("All4_InvMass_Charge-1", m, pt);
              }
            }
          }
        }

        if (fIsMC) {
          if (primLabelTrack > -1 && primLabel2 > -1) {
            // Track
            AliAODMCParticle* prim2 = (AliAODMCParticle*)fStack->At(primLabelTrack);
            Int_t primTrackPdg = prim2->GetPdgCode();
            TLorentzVector lv1;
            prim2->Momentum(lv1);
            //                 Double_t e1 = prim2->E();

            // Cluster
            AliAODMCParticle* prim3 = (AliAODMCParticle*)fStack->At(primLabel2);
            TLorentzVector lv2;
            prim3->Momentum(lv2);
            //                 Double_t e2 = prim2->E();

            Double_t mPrim = (lv1 + lv2).M();
            Double_t ptPrim = (lv1 + lv2).Pt();

            Int_t primLabelCluster1 = prim3->GetPdgCode();
            FillHistogram(Form("MC_InvMass_Charge%d", track->Charge()), mPrim, ptPrim);
            // CPV2 cut
            if (ph2->IsCPV2OK()) {
              FillHistogram(Form("CPV2_MC_InvMass_Charge%d", track->Charge()), mPrim, ptPrim);
            }
            if (ph2->IsDisp2OK()) {
              FillHistogram(Form("Disp2_MC_InvMass_Charge%d", track->Charge()), mPrim, ptPrim);
              if (ph2->IsCPV2OK()) {
                FillHistogram(Form("All3_MC_InvMass_Charge%d", track->Charge()), mPrim, ptPrim);
                if (TestDCAZ(trackPt, b[1]) && TestDCAXY(trackPt, b[0]) && abs(cpa) > 0.9) {
                  if (track->Charge() > 0 && rad > 0.2 && DCADaugPlus(trackPt, dca)) {
                    FillHistogram("All4_MC_InvMass_Charge1", mPrim, ptPrim);
                  }
                  if (track->Charge() < 0 && rad > 0.1 && DCADaugMinus(trackPt, dca)) {
                    FillHistogram("All4_MC_InvMass_Charge-1", mPrim, ptPrim);
                  }
                }
              }
            }
            // PDG pion and nbar
            if (primTrackPdg == 211 && primLabelCluster1 == -2112) {
              FillHistogram("MC_InvMass_AntiSigmaPlus", mPrim, ptPrim);
              FillHistogram("InvMass_AntiSigmaPlus", m, pt);
              // CPV2 cut
              if (ph2->IsCPV2OK()) {
                FillHistogram("CPV2_MC_InvMass_AntiSigmaPlus", mPrim, ptPrim);
                FillHistogram("CPV2_InvMass_AntiSigmaPlus", m, pt);
              }
              // Disp2 cut
              if (ph2->IsDisp2OK()) {
                FillHistogram("Disp2_MC_InvMass_AntiSigmaPlus", mPrim, ptPrim);
                FillHistogram("Disp2_InvMass_AntiSigmaPlus", m, pt);

                if (ph2->IsCPV2OK()) {
                  FillHistogram("All3_MC_InvMass_AntiSigmaPlus", mPrim, ptPrim);
                  FillHistogram("All3_InvMass_AntiSigmaPlus", m, pt);

                  if (TestDCAZ(trackPt, b[1]) && TestDCAXY(trackPt, b[0]) && abs(cpa) > 0.9 && rad > 0.2 &&
                      DCADaugPlus(trackPt, dca)) {
                    FillHistogram("All4_MC_InvMass_AntiSigmaPlus", mPrim, ptPrim);
                    FillHistogram("All4_InvMass_AntiSigmaPlus", m, pt);
                  }
                }
              }
              if (prim2->GetMother() > -1 && prim3->GetMother() > -1) {
                AliAODMCParticle* prim2parent = (AliAODMCParticle*)fStack->At(prim2->GetMother());
                AliAODMCParticle* prim3parent = (AliAODMCParticle*)fStack->At(prim3->GetMother());
                if (prim2parent->GetPdgCode() == -3112 && prim3parent->GetPdgCode() == -3112 &&
                    prim2parent == prim3parent) {
                  FillHistogram("MC_InvMass_AntiSigmaPlus_ParentCheck", mPrim, ptPrim);
                  FillHistogram("InvMass_AntiSigmaPlus_ParentCheck", m, pt);
                  FillHistogram("MC_APlus_4", lv2.Pt());
                  FillHistogram("Rec_APlus_4", ph2->Pt());
                  FillHistogram("PionPlus_3", trackPt);
                  // CPV2 cut
                  if (ph2->IsCPV2OK() == kTRUE) {
                    FillHistogram("CPV2_MC_InvMass_AntiSigmaPlus_ParentCheck", mPrim, ptPrim);
                    FillHistogram("CPV2_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
                    FillHistogram("MC_APlus_5", lv2.Pt());
                    FillHistogram("Rec_APlus_5", ph2->Pt());
                  }
                  // Disp2 cut
                  if (ph2->IsDisp2OK() == kTRUE) {
                    FillHistogram("Disp2_MC_InvMass_AntiSigmaPlus_ParentCheck", mPrim, ptPrim);
                    FillHistogram("Disp2_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
                    FillHistogram("MC_APlus_6", lv2.Pt());
                    FillHistogram("Rec_APlus_6", ph2->Pt());
                  }
                  // All3 cuts
                  if (ph2->IsCPV2OK() == kTRUE && ph2->IsDisp2OK() == kTRUE) {
                    FillHistogram("All3_MC_InvMass_AntiSigmaPlus_ParentCheck", mPrim, ptPrim);
                    FillHistogram("All3_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
                    FillHistogram("MC_APlus_7", lv2.Pt());
                    FillHistogram("Rec_APlus_7", ph2->Pt());
                    FillHistogram("PionPlus_4", trackPt);
                  }
                  // All4 cuts
                  if (ph2->IsCPV2OK() == kTRUE && (TestDCAZ(trackPt, b[1])) && (TestDCAXY(trackPt, b[0])) &&
                      ph2->IsDisp2OK() == kTRUE && abs(cpa) > 0.9 && rad > 0.2 && DCADaugPlus(trackPt, dca)) {
                    FillHistogram("All4_MC_InvMass_AntiSigmaPlus_ParentCheck", mPrim, ptPrim);
                    FillHistogram("All4_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
                    FillHistogram("MC_APlus_8", lv2.Pt());
                    FillHistogram("Rec_APlus_8", ph2->Pt());
                    FillHistogram("PionPlus_5", trackPt);
                  }
                }
              }
            }
            // PDG pion and nbar
            if (primTrackPdg == -211 && primLabelCluster1 == -2112) {
              FillHistogram("MC_InvMass_AntiSigmaMinus", mPrim, ptPrim);
              FillHistogram("InvMass_AntiSigmaMinus", m, pt);
              // CPV2 cut
              if (ph2->IsCPV2OK() == kTRUE) {
                FillHistogram("CPV2_MC_InvMass_AntiSigmaMinus", mPrim, ptPrim);
                FillHistogram("CPV2_InvMass_AntiSigmaMinus", m, pt);
              }
              // Disp2 cut
              if (ph2->IsDisp2OK() == kTRUE) {
                FillHistogram("Disp2_MC_InvMass_AntiSigmaMinus", mPrim, ptPrim);
                FillHistogram("Disp2_InvMass_AntiSigmaMinus", m, pt);
              }
              // All3 cuts
              if (ph2->IsCPV2OK() == kTRUE && ph2->IsDisp2OK() == kTRUE) {
                FillHistogram("All3_MC_InvMass_AntiSigmaMinus", mPrim, ptPrim);
                FillHistogram("All3_InvMass_AntiSigmaMinus", m, pt);
              }
              // All4 cuts
              if (ph2->IsCPV2OK() == kTRUE && (TestDCAZ(trackPt, b[1])) && (TestDCAXY(trackPt, b[0])) &&
                  ph2->IsDisp2OK() == kTRUE && abs(cpa) > 0.9 && rad > 0.1 && DCADaugMinus(trackPt, dca)) {
                FillHistogram("All4_MC_InvMass_AntiSigmaMinus", mPrim, ptPrim);
                FillHistogram("All4_InvMass_AntiSigmaMinus", m, pt);
              }
              // PDG Codes, pion and nbar
              if (prim2->GetMother() > -1 && prim3->GetMother() > -1) {
                AliAODMCParticle* prim2parent = (AliAODMCParticle*)fStack->At(prim2->GetMother());
                AliAODMCParticle* prim3parent = (AliAODMCParticle*)fStack->At(prim3->GetMother());
                if (prim2parent->GetPdgCode() == -3222 && prim3parent->GetPdgCode() == -3222 &&
                    prim2parent == prim3parent) {
                  FillHistogram("MC_InvMass_AntiSigmaMinus_ParentCheck", mPrim, ptPrim);
                  FillHistogram("InvMass_AntiSigmaMinus_ParentCheck", m, pt);
                  FillHistogram("MC_AMinus_4", lv2.Pt());
                  FillHistogram("Rec_AMinus_4", ph2->Pt());
                  FillHistogram("PionMinus_3", trackPt);
                  // CPV2 cut
                  if (ph2->IsCPV2OK()) {
                    FillHistogram("CPV2_MC_InvMass_AntiSigmaMinus_ParentCheck", mPrim, ptPrim);
                    FillHistogram("CPV2_InvMass_AntiSigmaMinus_ParentCheck", m, pt);
                    FillHistogram("MC_AMinus_5", lv2.Pt());
                    FillHistogram("Rec_AMinus_5", ph2->Pt());
                  }
                  // Disp2 cut
                  if (ph2->IsDisp2OK()) {
                    FillHistogram("Disp2_MC_InvMass_AntiSigmaMinus_ParentCheck", mPrim, ptPrim);
                    FillHistogram("Disp2_InvMass_AntiSigmaMinus_ParentCheck", m, pt);
                    FillHistogram("MC_AMinus_6", lv2.Pt());
                    FillHistogram("Rec_AMinus_6", ph2->Pt());
                  }
                  // All3 cuts
                  if (ph2->IsCPV2OK() && ph2->IsDisp2OK()) {
                    FillHistogram("All3_MC_InvMass_AntiSigmaMinus_ParentCheck", mPrim, ptPrim);
                    FillHistogram("All3_InvMass_AntiSigmaMinus_ParentCheck", m, pt);
                    FillHistogram("MC_AMinus_7", lv2.Pt());
                    FillHistogram("Rec_AMinus_7", ph2->Pt());
                    FillHistogram("PionMinus_4", trackPt);
                  }
                  // All4 cuts
                  if (ph2->IsCPV2OK() == kTRUE && (TestDCAZ(trackPt, b[1])) && (TestDCAXY(trackPt, b[0])) &&
                      ph2->IsDisp2OK() == kTRUE && abs(cpa) > 0.9 && rad > 0.1 && DCADaugMinus(trackPt, dca)) {
                    FillHistogram("All4_MC_InvMass_AntiSigmaMinus_ParentCheck", mPrim, ptPrim);
                    FillHistogram("All4_InvMass_AntiSigmaMinus_ParentCheck", m, pt);
                    FillHistogram("MC_AMinus_8", lv2.Pt());
                    FillHistogram("Rec_AMinus_8", ph2->Pt());
                    FillHistogram("PionMinus_5", trackPt);
                  }
                }
              }
            }
            if (primLabelCluster1 == -2112) {
              FillHistogram(Form("MC_InvMass_AntiNeutron_Charge%d", track->Charge()), mPrim, ptPrim);
              FillHistogram(Form("InvMass_AntiNeutron_Charge%d", track->Charge()), m, pt);
            }
            // Track- pion
            if (primTrackPdg == -211) {
              FillHistogram("MC_InvMass_PiMinus", mPrim, ptPrim);
              FillHistogram("InvMass_PiMinus", m, pt);
            }
            if (primTrackPdg == 211) {
              FillHistogram("MC_InvMass_PiPlus", mPrim, ptPrim);
              FillHistogram("InvMass_PiPlus", m, pt);
            }
          }
        }
      }
    }
  }
}

void AliAnalysisSigmaBarCharged::FillHistogram(const char* key, Double_t x) const
{
  TH1* hist = dynamic_cast<TH1*>(fOutputContainer->FindObject(key));
  if (hist) {
    hist->Fill(x);
  } else {
    AliError(Form("can not find histogram (of instance TH1) <%s> ", key));
  }
}

void AliAnalysisSigmaBarCharged::FillHistogram(const char* key, Double_t x, Double_t y) const
{
  TObject* obj = fOutputContainer->FindObject(key);
  TH1* th1 = dynamic_cast<TH1*>(obj);
  if (th1) {
    th1->Fill(x, y);
  } else {
    TH2* th2 = dynamic_cast<TH2*>(obj);
    if (th2) {
      th2->Fill(x, y);
    } else {
      AliError(Form("can not find histogram (of instance TH1) <%s> ", key));
    }
  }
}

void AliAnalysisSigmaBarCharged::FillHistogram(const char* key, Double_t x, Double_t y, Double_t z) const
{
  TObject* obj = fOutputContainer->FindObject(key);

  TH2* th2 = dynamic_cast<TH2*>(obj);
  if (th2) {
    th2->Fill(x, y, z);
    return;
  }

  TH3* th3 = dynamic_cast<TH3*>(obj);
  if (th3) {
    th3->Fill(x, y, z);
    return;
  }

  AliError(Form("can not findi histogram (of instance TH2) <%s> ", key));
}
void AliAnalysisSigmaBarCharged::Terminate(Option_t*) {}

Double_t AliAnalysisSigmaBarCharged::RealRes(Double_t x)
{
  // Simulates time resolution
  if (x <= 0.2) {
    x = 0.2;
  }
  if (x >= 10) {
    x = 10;
  }
  return sqrt(pow(9.217682 * TMath::Exp(-x / 3.575588e+11) / x - 0.02196338 * x + 0.02202962 * x * x, 2) -
              pow(0.5, 2)) *
         1e-9;
}

Bool_t AliAnalysisSigmaBarCharged::TestDCAXY(Double_t pt, Double_t dca)
{
  // Checks if track DCA far enough from primary vertex due to large Sigma lifetime
  if (pt > 0.5)
    return true;
  if (TMath::Abs(dca) >= -5.43674e-02 + 3.14319e-02 / pt) { // dca>=f(pt)
    return true;
  } else {
    return false;
  }
}

Bool_t AliAnalysisSigmaBarCharged::TestDCAZ(Double_t pt, Double_t dca)
{
  // Checks if track DCA in z direction far enough from primary vertex due to large Sigma lifetime
  if (pt > 0.5)
    return true;
  if (TMath::Abs(dca) >= -3.65702e-02 + 5.38728e-02 / pt) {
    return true;
  } else {
    return false;
  }
}

Bool_t AliAnalysisSigmaBarCharged::DCADaugPlus(Double_t pt, Double_t dca)
{
  // Daughter DCA cut
  if (pt < 0.25 && dca < 0.2) {
    return false;
  } else {
    return true;
  }
}

Bool_t AliAnalysisSigmaBarCharged::DCADaugMinus(Double_t pt, Double_t dca)
{
  if (pt < 0.25 && dca < 0.1) {
    return false;
  } else {
    return true;
  }
}
