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

// Analysis task for antineutron measurement in PHOS and mesurment of charged
// Sigma bar Authors: Pavel Gordeev, Dmitri Blau, Dmitri Peresunko
#include "AliAnalysisSigmaBarCharged.h"

#include <cstdio>
#include <iostream>

#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliCaloPhoton.h"
#include "AliCascadeVertexer.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliEventplane.h"
#include "AliExternalTrackParam.h"
#include "AliGenEventHeader.h"
#include "AliLog.h"
#include "AliMultSelection.h"
#include "AliOADBContainer.h"
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

AliAnalysisSigmaBarCharged::AliAnalysisSigmaBarCharged(const char *name)
    : AliAnalysisTaskSE(name), fOutputContainer(nullptr), fPIDResponse(nullptr),
      fUtils(nullptr), fPHOSgeom(nullptr), fRunNumber(0), fGamma(nullptr),
      fCurrentMixedList(nullptr), fCluTimeCut(150.e-9), fCluNbarMinE(0.6),
      fCPAplusCut(0.), fCPAminusCut(0.), fDCAdaugplusCut(0.06),
      fDCAdaugminusCut(0.06), fRADplusCut(0.25), fRADminusCut(0.15),
      fCPVCut(10.), fDispCut(4.), fNcellCut(7), fTracksBits(4), fTrackEta(0.8),
      fNTPCclusters(60), fTPCsigmas(3), fIsMC(false), fAdditionHist(false),
      fInvMassHist(false), fQAhist(false), fPHOSClusterTOFOption(0),
      fCentrality(0), fCentBin(0), fStack(nullptr), fEvent(nullptr) {
  for (Int_t i = 0; i < 10; i++) {
    for (Int_t j = 0; j < 10; j++) {
      fPHOSEvents[i][j] = 0x0; // Container for PHOS photons
    }
  }
  for (Int_t mod = 0; mod < 6; mod++)
    fPHOSBadMap[mod] = 0x0;

  // Output slots #0 write into a TH1 container
  DefineOutput(1, THashList::Class());
}

AliAnalysisSigmaBarCharged::AliAnalysisSigmaBarCharged(
    const AliAnalysisSigmaBarCharged &rh)
    : AliAnalysisTaskSE(rh.GetName()), fOutputContainer(nullptr),
      fPIDResponse(nullptr), fUtils(nullptr), fPHOSgeom(nullptr), fRunNumber(0),
      fGamma(nullptr), fCurrentMixedList(nullptr), fCluTimeCut(150.e-9),
      fCluNbarMinE(0.6), fCPAplusCut(0.), fCPAminusCut(0.),
      fDCAdaugplusCut(0.06), fDCAdaugminusCut(0.06), fRADplusCut(0.25),
      fRADminusCut(0.15), fCPVCut(10.), fDispCut(4.), fNcellCut(7),
      fTracksBits(4), fTrackEta(0.8), fNTPCclusters(60), fTPCsigmas(3),
      fIsMC(false), fAdditionHist(false), fInvMassHist(false), fQAhist(false),
      fPHOSClusterTOFOption(0), fCentrality(0), fCentBin(0), fStack(nullptr),
      fEvent(nullptr) {
  for (Int_t i = 0; i < 10; i++) {
    for (Int_t j = 0; j < 10; j++) {
      fPHOSEvents[i][j] = 0x0; // Container for PHOS photons
    }
  }
  for (Int_t mod = 0; mod < 6; mod++)
    fPHOSBadMap[mod] = 0x0;

  if (fOutputContainer)
    delete fOutputContainer;
  fOutputContainer = new THashList();
}

AliAnalysisSigmaBarCharged &
AliAnalysisSigmaBarCharged::operator=(const AliAnalysisSigmaBarCharged &rh) {
  this->~AliAnalysisSigmaBarCharged();
  new (this) AliAnalysisSigmaBarCharged(rh);
  return *this;
}

AliAnalysisSigmaBarCharged::~AliAnalysisSigmaBarCharged() {
  if (fOutputContainer) {
    delete fOutputContainer;
    fOutputContainer = 0x0;
  }
  for (Int_t i = 0; i < 10; i++) {
    for (Int_t j = 0; j < 10; j++) {
      if (fPHOSEvents[i][j]) {
        delete fPHOSEvents[i][j];
        fPHOSEvents[i][j] = 0x0;
      }
    }
  }
  for (Int_t mod = 0; mod < 6; mod++) {
    if (fPHOSBadMap[mod]) {
      delete fPHOSBadMap[mod];
      fPHOSBadMap[mod] = 0x0;
    }
  }
}

void AliAnalysisSigmaBarCharged::UserCreateOutputObjects() {
  if (fOutputContainer != nullptr) {
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  // Criteria
  char cPT[13][12];
  snprintf(cPT[0], 12, "Photon");
  snprintf(cPT[1], 12, "Electron");
  snprintf(cPT[2], 12, "Positron");
  snprintf(cPT[3], 12, "Proton");
  snprintf(cPT[4], 12, "AntiProton");
  snprintf(cPT[5], 12, "PiPlus");
  snprintf(cPT[6], 12, "PiMinus");
  snprintf(cPT[7], 12, "Neutron");
  snprintf(cPT[8], 12, "AntiNeutron");
  snprintf(cPT[9], 12, "KPlus");
  snprintf(cPT[10], 12, "KMinus");
  snprintf(cPT[11], 12, "KLong");
  snprintf(cPT[12], 12, "Other");

  // Binning
  int dispmin = 0, dispmax = 15, dispbins = 600; // dispersion: m02 and m20
  double tofmin = -250.e-9, tofmax = 250.e-9;
  int tofbins = 2000; // Time Of Flight
  double dcaxymin = -2.4, dcaxymax = 2.4, dcazmin = -3.2, dcazmax = 3.2;
  int dcabins = 800; // DCA XY and DCA Z
  double invmin = 1.0, invmax = 3.0;
  int invbins = 2000;                            // invariant mass
  int recmin = -20, recmax = 20, recbins = 1000; // reconstructed momentum
  int dcamin = 0, dcamax = 15;                   // dca between daughters
  int cpamin = -1, cpamax = 1, cpabins = 400;    // cosine of pointing angle
  int radmin = 0, radmax = 15,
      radbins = 800; // distance between primary and secondary vertex
  int ebins = 800;
  double emin = 0., emax = 20.;

  if (fIsMC) {
    fOutputContainer->Add(
        new TH1F("MC_AntiSigmaPlus_1",
                 "Spectrum of AntiSigmaPlus;#it{p}_{T} (GeV/#it{c});Counts",
                 ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("MC_AntiSigmaMinus_1",
                 "Spectrum of AntiSigmaMinus;#it{p}_{T} (GeV/#it{c});Counts",
                 ebins, emin, emax));
    if (fAdditionHist) {
      fOutputContainer->Add(new TH1F(
          "MC_Nbar_1", "Spectrum of nbar;#it{p}_{T} (GeV/#it{c});Counts", ebins,
          emin, emax));
      fOutputContainer->Add(new TH1F(
          "MC_Pbar_1", "Spectrum of pbar;#it{p}_{T} (GeV/#it{c});Counts", ebins,
          emin, emax));
      fOutputContainer->Add(new TH1F(
          "MC_Nbar_prim_1", "Spectrum of nbar;#it{p}_{T} (GeV/#it{c});Counts",
          ebins, emin, emax));
      fOutputContainer->Add(
          new TH1F("MC_Nbar_nonprim_1",
                   "Spectrum of nbar;#it{p}_{T} (GeV/#it{c});Counts", ebins,
                   emin, emax));
    }
    if (fQAhist) {
      fOutputContainer->Add(new TH1F(
          "MC_Pion_1", "Spectrum of pion;#it{p}_{T} (GeV/#it{c});Counts", ebins,
          emin, emax));
    }
  }

  if (fQAhist) {
    for (Int_t iModule = 1; iModule <= 4; iModule++) {
      fOutputContainer->Add(new TH2F(Form("OccupancyPlot_%dM", iModule),
                                     ";x, cells; z, cells", 64, 0, 64, 56, 0,
                                     56));
      fOutputContainer->Add(
          new TH2F(Form("ClusterTOFvsE_%dM", iModule),
                   "Cluster time vs energy;TOF (s);#it{E}_{clu} (GeV)", tofbins,
                   tofmin, tofmax, ebins, emin, emax));
    }
    fOutputContainer->Add(
        new TH2F("Acceptance_track", "Spec", 400, -10, 10, 800, -1, 7));
    fOutputContainer->Add(
        new TH2F("PhivsPt_track_1", "Spec", ebins, emin, emax, 800, -1, 7));
    fOutputContainer->Add(
        new TH2F("PhivsPt_track_2", "Spec", ebins, emin, emax, 800, -1, 7));
    fOutputContainer->Add(
        new TH1F("Spec_track_FB5_0", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spec_track_FB5_1", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spec_track_FB5_2", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spec_track_FB5_3", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spec_track_FB5_4", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spec_track_FB4_1", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spec_track_FB4_2", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spec_track_FB4_3", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spec_track_FB4_4", "Spec", ebins, emin, emax));
  }

  // Time of flight vs energy
  fOutputContainer->Add(new TH2F(
      "hClusterTOFvsE", "Cluster time vs energy;TOF (s);#it{E}_{clu} (GeV)",
      tofbins, tofmin, tofmax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F(
      "hClusterTOFvsE_1", "Cluster time vs energy;TOF (s);#it{E}_{clu} (GeV)",
      tofbins, tofmin, tofmax, ebins, emin, emax));
  fOutputContainer->Add(
      new TH2F("hClusterpTvsE_1",
               "Reconstructed #it{p}_{T} vs energy;#it{p}_{T, rec} "
               "(GeV/#it{c});#it{E}_{clu} (GeV)",
               ebins, emin, emax, ebins, emin, emax));
  if (fAdditionHist) {
    fOutputContainer->Add(
        new TH2F("hClusterTOFvsE_cut1",
                 "Cluster time vs energy;TOF (s);#it{E}_{clu} (GeV)", tofbins,
                 tofmin, tofmax, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("hClusterTOFvsE_cut2",
                 "Cluster time vs energy;TOF (s);#it{E}_{clu} (GeV)", tofbins,
                 tofmin, tofmax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F(
        "hClusterTOFvsE_2", "Cluster time vs energy;TOF (s);#it{E}_{clu} (GeV)",
        tofbins, tofmin, tofmax, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("hClusterpTvsE_2",
                 "Reconstructed #it{p}_{T} vs energy;#it{p}_{T, rec} "
                 "(GeV/#it{c});#it{E}_{clu} (GeV)",
                 ebins, emin, emax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("All_Disp",
                                   "Dispersion of all particles;M20 "
                                   "(cm^{2});M02 (cm^{2});#it{E}_{clu} (GeV)",
                                   300, 0, 12, 300, 0, 12));
    // fOutputContainer->Add(new TH2F("All_M02", "Dispersion of all
    // particles;M02 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax)); fOutputContainer->Add(new TH2F("All_M20", "Dispersion of all
    // particles;M20 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax));
    fOutputContainer->Add(new TH2F(
        "All_Sum",
        "Dispersion of all particles;M02 + M20 (cm^{2});#it{E}_{clu} (GeV)",
        600, 0, 24, ebins, emin, emax));
    // fOutputContainer->Add(new TH2F("Nbar_all_M02_1", "Dispersion of all
    // particles;M02 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax)); fOutputContainer->Add(new TH2F("Nbar_all_M20_1", "Dispersion of
    // all particles;M20 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax));
    fOutputContainer->Add(new TH2F(
        "Nbar_all_Sum_1",
        "Dispersion of all particles;M02 + M20 (cm^{2});#it{E}_{clu} (GeV)",
        600, 0, 24, ebins, emin, emax));
    // fOutputContainer->Add(new TH2F("Nbar_all_M02_2", "Dispersion of all
    // particles;M02 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax)); fOutputContainer->Add(new TH2F("Nbar_all_M20_2", "Dispersion of
    // all particles;M20 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax));
    fOutputContainer->Add(new TH2F(
        "Nbar_all_Sum_2",
        "Dispersion of all particles;M02 + M20 (cm^{2});#it{E}_{clu} (GeV)",
        600, 0, 24, ebins, emin, emax));
    // fOutputContainer->Add(new TH2F("Pbar_all_M02_1", "Dispersion of all
    // particles;M02 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax)); fOutputContainer->Add(new TH2F("Pbar_all_M20_1", "Dispersion of
    // all particles;M20 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax)); fOutputContainer->Add(new TH2F("Pbar_all_Sum_1", "Dispersion of
    // all particles;M02 + M20 (cm^{2});#it{E}_{clu} (GeV)", 600, 0, 24, ebins,
    // emin, emax)); fOutputContainer->Add(new TH2F("Pbar_all_M02_2",
    // "Dispersion of all particles;M02 (cm^{2});#it{E}_{clu} (GeV)", 300, 0,
    // 12, ebins, emin, emax)); fOutputContainer->Add(new TH2F("Pbar_all_M20_2",
    // "Dispersion of all particles;M20 (cm^{2});#it{E}_{clu} (GeV)", 300, 0,
    // 12, ebins, emin, emax)); fOutputContainer->Add(new TH2F("Pbar_all_Sum_2",
    // "Dispersion of all particles;M02 + M20 (cm^{2});#it{E}_{clu} (GeV)", 600,
    // 0, 24, ebins, emin, emax));

    // Spec of clu and nbar,pbar candidates
    fOutputContainer->Add(
        new TH1F("Spectrum_pt_allnbar_1", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_eclu_allnbar_1", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_pt_allnbar_2", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_eclu_allnbar_2", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_pt_allpbar_1", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_eclu_allpbar_1", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_pt_allpbar_2", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_eclu_allpbar_2", "Spec", ebins, emin, emax));

    // TPC PID
    fOutputContainer->Add(new TH2F(
        "TPC_clusters_track", "TPC clusters;Clusters;#it{p}_{T} (GeV/#it{c})",
        160, 0, 160, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("TPC_pid_track",
                                   "TPC n#sigma;n#sigma;#it{p} (GeV/#it{c})",
                                   80, -10, 10, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("TPC_pid_proton",
                                   "TPC n#sigma;n#sigma;#it{p} (GeV/#it{c})",
                                   80, -10, 10, ebins, emin, emax));
  }

  if (fIsMC && fAdditionHist) {
    fOutputContainer->Add(
        new TH2F("hClusterTOFvsE_nbar",
                 "Cluster time vs energy;TOF (s);#it{E}_{clu} (GeV)", tofbins,
                 tofmin, tofmax, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("hClusterTOFvsE_nbar_1",
                 "Cluster time vs energy;TOF (s);#it{E}_{clu} (GeV)", tofbins,
                 tofmin, tofmax, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("hClusterpTvsE_nbar_1",
                 "Reconstructed #it{p}_{T} vs energy;#it{p}_{T, rec} "
                 "(GeV/#it{c});#it{E}_{clu} (GeV)",
                 ebins, emin, emax, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("hClusterpTmcvsE_nbar_1",
                 "Reconstructed #it{p}_{T} vs energy;#it{p}_{T, rec} "
                 "(GeV/#it{c});#it{E}_{clu} (GeV)",
                 ebins, emin, emax, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("hClusterTOFvsE_nbar_2",
                 "Cluster time vs energy;TOF (s);#it{E}_{clu} (GeV)", tofbins,
                 tofmin, tofmax, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("hClusterpTvsE_nbar_2",
                 "Reconstructed #it{p}_{T} vs energy;#it{p}_{T, rec} "
                 "(GeV/#it{c});#it{E}_{clu} (GeV)",
                 ebins, emin, emax, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("hClusterpTmcvsE_nbar_2",
                 "Reconstructed #it{p}_{T} vs energy;#it{p}_{T, rec} "
                 "(GeV/#it{c});#it{E}_{clu} (GeV)",
                 ebins, emin, emax, ebins, emin, emax));

    fOutputContainer->Add(new TH2F(
        "AntiProton_Disp", "Dispersion of pbar;M20 (cm^{2});M02 (cm^{2})", 300,
        0, 12, 300, 0, 12));
    fOutputContainer->Add(new TH2F(
        "Photon_Disp", "Dispersion of photon;M20 (cm^{2});M02 (cm^{2})", 300, 0,
        12, 300, 0, 12));
    fOutputContainer->Add(new TH2F(
        "AntiNeutron_Disp",
        "Dispersion of nbar;M20 (cm^{2});M02 (cm^{2});#it{E}_{clu} (GeV)", 300,
        0, 12, 300, 0, 12));

    fOutputContainer->Add(
        new TH2F("EcluvsP_nbar",
                 "Response matrix;#it{E}_{clu} (GeV);#it{p}_{MC} (GeV/#it{c})",
                 ebins, emin, emax, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("EcluvsP_pbar",
                 "Response matrix;#it{E}_{clu} (GeV);#it{p}_{MC} (GeV/#it{c})",
                 ebins, emin, emax, ebins, emin, emax));

    fOutputContainer->Add(
        new TH2F("Nbar_pmcvsprec-pmc_prim",
                 "nbar_pmcvsprec-pmc;#it{p}_{MC} (GeV/#it{c});#it{p}_{rec} "
                 "#minus #it{p}_{MC} (GeV/#it{c})",
                 ebins, emin, emax, recbins, recmin, recmax));
    fOutputContainer->Add(
        new TH2F("Nbar_prec-pmcvseclu_prim",
                 "nbar_pmcvsprec-pmc;#it{E}_{clu} (GeV);#it{p}_{rec} #minus "
                 "#it{p}_{MC} (GeV/#it{c})",
                 ebins, emin, emax, recbins, recmin, recmax));
    fOutputContainer->Add(new TH2F(
        "Nbar_pmcvsprec_prim",
        "nbar_pmcvsprec-pmc;#it{p}_{MC} (GeV/#it{c});#it{p}_{rec} (GeV/#it{c})",
        ebins, emin, emax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F(
        "Nbar_emcvserec_prim",
        "nbar_emcvserec-emc;#it{E}_{MC} (GeV/#it{c});#it{E}_{rec} (GeV/#it{c})",
        ebins, emin, emax, ebins, emin, emax));

    fOutputContainer->Add(
        new TH2F("Nbar_pmcvsprec-pmc",
                 "nbar_pmcvsprec-pmc;#it{p}_{MC} (GeV/#it{c});#it{p}_{rec} "
                 "#minus #it{p}_{MC} (GeV/#it{c})",
                 ebins, emin, emax, recbins, recmin, recmax));
    fOutputContainer->Add(
        new TH2F("Nbar_prec-pmcvseclu",
                 "nbar_pmcvsprec-pmc;#it{E}_{clu} (GeV);#it{p}_{rec} #minus "
                 "#it{p}_{MC} (GeV/#it{c})",
                 ebins, emin, emax, recbins, recmin, recmax));
    fOutputContainer->Add(new TH2F(
        "Nbar_pmcvsprec",
        "nbar_pmcvsprec-pmc;#it{p}_{MC} (GeV/#it{c});#it{p}_{rec} (GeV/#it{c})",
        ebins, emin, emax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F(
        "Nbar_emcvserec",
        "nbar_emcvserec-emc;#it{E}_{MC} (GeV/#it{c});#it{E}_{rec} (GeV/#it{c})",
        ebins, emin, emax, ebins, emin, emax));

    fOutputContainer->Add(new TH1F(
        "Spectrum_Sum", "Spectrum of clusters;#it{E}_{clu} (GeV);Counts", ebins,
        emin, emax));
    fOutputContainer->Add(new TH1F(
        "Ncells_Sum", "Ncells of clusters;#it{N}_{cells};Counts", 50, 0, 50));
    fOutputContainer->Add(new TH1F(
        "Spectrum_Sum_1", "Spectrum of clusters;#it{E}_{clu} (GeV);Counts",
        ebins, emin, emax));
    fOutputContainer->Add(new TH1F(
        "Ncells_Sum_1", "Ncells of clusters;#it{N}_{cells};Counts", 50, 0, 50));
    for (int iPID = 0; iPID < 13; iPID++) {
      fOutputContainer->Add(new TH1F(
          Form("Spectrum_%s", cPT[iPID]),
          "Spectrum of clusters;#it{E}_{clu} (GeV);Counts", ebins, emin, emax));
      fOutputContainer->Add(
          new TH1F(Form("Ncells_%s", cPT[iPID]),
                   "Spectrum of clusters;#it{N}_{cells};Counts", 50, 0, 50));
      fOutputContainer->Add(new TH1F(
          Form("Spectrum_%s_1", cPT[iPID]),
          "Spectrum of clusters;#it{E}_{clu} (GeV);Counts", ebins, emin, emax));
      fOutputContainer->Add(
          new TH1F(Form("Ncells_%s_1", cPT[iPID]),
                   "Spectrum of clusters;#it{N}_{cells};Counts", 50, 0, 50));
    }
    // fOutputContainer->Add(new TH1F("PDG_Other", "pdg", 8000, -4000, 4000));

    // Spectre for ROC
    // fOutputContainer->Add(new TH1F("Spectrum_Nbar", "Spec", ebins, emin,
    // emax)); Clusters default cuts and variations
    fOutputContainer->Add(
        new TH1F("Spectrum_d_0_0_0", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_0_d_0_0", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_0_0_d_0", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_0_0_0_d", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_Nbar_d_0_0_0", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_Nbar_0_d_0_0", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_Nbar_0_0_d_0", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_Nbar_0_0_0_d", "Spec", ebins, emin, emax));
    for (Int_t i = 1; i < 14; i++) {
      fOutputContainer->Add(
          new TH1F(Form("Spectrum_%d_0_0_0", i), "Spec", ebins, emin, emax));
      fOutputContainer->Add(
          new TH1F(Form("Spectrum_0_%d_0_0", i), "Spec", ebins, emin, emax));
      fOutputContainer->Add(
          new TH1F(Form("Spectrum_0_0_%d_0", i), "Spec", ebins, emin, emax));
      fOutputContainer->Add(
          new TH1F(Form("Spectrum_0_0_0_%d", i), "Spec", ebins, emin, emax));
      fOutputContainer->Add(new TH1F(Form("Spectrum_Nbar_%d_0_0_0", i), "Spec",
                                     ebins, emin, emax));
      fOutputContainer->Add(new TH1F(Form("Spectrum_Nbar_0_%d_0_0", i), "Spec",
                                     ebins, emin, emax));
      fOutputContainer->Add(new TH1F(Form("Spectrum_Nbar_0_0_%d_0", i), "Spec",
                                     ebins, emin, emax));
      fOutputContainer->Add(new TH1F(Form("Spectrum_Nbar_0_0_0_%d", i), "Spec",
                                     ebins, emin, emax));
    }

    fOutputContainer->Add(
        new TH2F("Nbar_pmcvsprec-pmc_plus",
                 "nbar_pmcvsprec-pmc;#it{p}_{MC} (GeV/#it{c});#it{p}_{rec} "
                 "#minus #it{p}_{MC} (GeV/#it{c})",
                 ebins, emin, emax, recbins, recmin, recmax));
    fOutputContainer->Add(
        new TH2F("Nbar_prec-pmcvseclu_plus",
                 "nbar_pmcvsprec-pmc;#it{E}_{clu} (GeV);#it{p}_{rec} #minus "
                 "#it{p}_{MC} (GeV/#it{c})",
                 ebins, emin, emax, recbins, recmin, recmax));
    fOutputContainer->Add(new TH2F(
        "Nbar_pmcvsprec_plus",
        "nbar_pmcvsprec-pmc;#it{p}_{MC} (GeV/#it{c});#it{p}_{rec} (GeV/#it{c})",
        ebins, emin, emax, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("Nbar_pmcvsprec-pmc_minus",
                 "nbar_pmcvsprec-pmc;#it{p}_{MC} (GeV/#it{c});#it{p}_{rec} "
                 "#minus #it{p}_{MC} (GeV/#it{c})",
                 ebins, emin, emax, recbins, recmin, recmax));
    fOutputContainer->Add(
        new TH2F("Nbar_prec-pmcvseclu_minus",
                 "nbar_pmcvsprec-pmc;#it{E}_{clu} (GeV);#it{p}_{rec} #minus "
                 "#it{p}_{MC} (GeV/#it{c})",
                 ebins, emin, emax, recbins, recmin, recmax));
    fOutputContainer->Add(new TH2F(
        "Nbar_pmcvsprec_minus",
        "nbar_pmcvsprec-pmc;#it{p}_{MC} (GeV/#it{c});#it{p}_{rec} (GeV/#it{c})",
        ebins, emin, emax, ebins, emin, emax));

    fOutputContainer->Add(
        new TH2F("SV_plus", "SV", 1000, 0., 0.1, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("SV_minus", "SV", 1000, 0., 0.1, ebins, emin, emax));

    fOutputContainer->Add(
        new TH2F("t_compare_plus", "compare t", 500, -1, 1, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("t_compare_plus_1", "compare t", 1000, -0.05,
                                   0.05, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("t_compare_plus_2", "compare t", 1000, 0.,
                                   0.25, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("t_compare_plus_3", "compare t", 1000, -0.1,
                                   0.1, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("t_compare_minus", "compare t", 500, -1, 1,
                                   ebins, emin, emax));
    fOutputContainer->Add(new TH2F("t_compare_minus_1", "compare t", 1000,
                                   -0.05, 0.05, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("t_compare_minus_2", "compare t", 1000, 0.,
                                   0.25, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("t_compare_minus_3", "compare t", 1000, -0.1,
                                   0.1, ebins, emin, emax));

    fOutputContainer->Add(
        new TH2F("TPC_clusters_pionsigma",
                 "TPC clusters;Clusters;#it{p}_{T} (GeV/#it{c})", 160, 0, 160,
                 ebins, emin, emax));

    // Sigma cut varitions and default
    fOutputContainer->Add(
        new TH1F("Sigma_Plus_DCAdaug", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Sigma_Plus_RAD", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Sigma_Plus_CPA", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Track_Plus_DCAdaug", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Track_Plus_RAD", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Track_Plus_CPA", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Sigma_Minus_DCAdaug", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Sigma_Minus_RAD", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Sigma_Minus_CPA", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Track_Minus_DCAdaug", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Track_Minus_RAD", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Track_Minus_CPA", "Spec", ebins, emin, emax));
    for (Int_t i = 1; i < 11; i++) {
      fOutputContainer->Add(new TH1F(Form("Sigma_Plus_DCAdaug_%d", i), "Spec",
                                     ebins, emin, emax));
      fOutputContainer->Add(
          new TH1F(Form("Sigma_Plus_RAD_%d", i), "Spec", ebins, emin, emax));
      fOutputContainer->Add(
          new TH1F(Form("Sigma_Plus_CPA_%d", i), "Spec", ebins, emin, emax));
      fOutputContainer->Add(new TH1F(Form("Track_Plus_DCAdaug_%d", i), "Spec",
                                     ebins, emin, emax));
      fOutputContainer->Add(
          new TH1F(Form("Track_Plus_RAD_%d", i), "Spec", ebins, emin, emax));
      fOutputContainer->Add(
          new TH1F(Form("Track_Plus_CPA_%d", i), "Spec", ebins, emin, emax));
      fOutputContainer->Add(new TH1F(Form("Sigma_Minus_DCAdaug_%d", i), "Spec",
                                     ebins, emin, emax));
      fOutputContainer->Add(
          new TH1F(Form("Sigma_Minus_RAD_%d", i), "Spec", ebins, emin, emax));
      fOutputContainer->Add(
          new TH1F(Form("Sigma_Minus_CPA_%d", i), "Spec", ebins, emin, emax));
      fOutputContainer->Add(new TH1F(Form("Track_Minus_DCAdaug_%d", i), "Spec",
                                     ebins, emin, emax));
      fOutputContainer->Add(
          new TH1F(Form("Track_Minus_RAD_%d", i), "Spec", ebins, emin, emax));
      fOutputContainer->Add(
          new TH1F(Form("Track_Minus_CPA_%d", i), "Spec", ebins, emin, emax));
    }

    fOutputContainer->Add(new TH2F("DCA_Daughters_AntiSigmaPlus",
                                   "DCA between daughters", dcabins, dcamin,
                                   dcamax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("DCA_Daughters_AntiSigmaMinus",
                                   "DCA between daughters", dcabins, dcamin,
                                   dcamax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("CPA_Daughters_AntiSigmaPlus",
                                   "DCA between daughters", cpabins, cpamin,
                                   cpamax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("CPA_Daughters_AntiSigmaMinus",
                                   "DCA between daughters", cpabins, cpamin,
                                   cpamax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("RAD_Daughters_AntiSigmaPlus",
                                   "RAD between daughters", radbins, radmin,
                                   radmax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("RAD_Daughters_AntiSigmaMinus",
                                   "RAD between daughters", radbins, radmin,
                                   radmax, ebins, emin, emax));

    fOutputContainer->Add(new TH2F("DCA_XY_pionsigma", "DCA XY", 480, -2.4, 2.4,
                                   ebins, emin, emax));
    fOutputContainer->Add(new TH2F("DCA_Z_pionsigma", "DCA Z", 640, -3.2, 3.2,
                                   ebins, emin, emax));

    fOutputContainer->Add(
        new TH2F("SV_x_plus", "SV_x", 800, -20, 20, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("SV_y_plus", "SV_y", 800, -20, 20, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("SV_z_plus", "SV_z", 800, -20, 20, ebins, emin, emax));

    fOutputContainer->Add(
        new TH2F("SV_x_minus", "SV_x", 800, -20, 20, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("SV_y_minus", "SV_y", 800, -20, 20, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("SV_z_minus", "SV_z", 800, -20, 20, ebins, emin, emax));

    // Proton spec from clu and nbar;
    // Spec of clu and nbar,pbar candidates
    fOutputContainer->Add(
        new TH1F("Spectrum_pt_nbar_1", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_ptmc_nbar_1", "Spec", ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Response_nbar_1", "Spec", ebins, emin, emax,
                                   ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_eclu_nbar_1", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_pt_nbar_2", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_ptmc_nbar_2", "Spec", ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Response_nbar_2", "Spec", ebins, emin, emax,
                                   ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_eclu_nbar_2", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_pt_pbar_1", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_ptmc_pbar_1", "Spec", ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Response_pbar_1", "Spec", ebins, emin, emax,
                                   ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_eclu_pbar_1", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_pt_pbar_2", "Spec", ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_ptmc_pbar_2", "Spec", ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Response_pbar_2", "Spec", ebins, emin, emax,
                                   ebins, emin, emax));
    fOutputContainer->Add(
        new TH1F("Spectrum_eclu_pbar_2", "Spec", ebins, emin, emax));
    // fOutputContainer->Add(new TH2F("Nbar_M02_1", "Dispersion of all
    // particles;M02 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax)); fOutputContainer->Add(new TH2F("Nbar_M20_1", "Dispersion of all
    // particles;M20 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax));
    fOutputContainer->Add(new TH2F(
        "Nbar_Sum_1",
        "Dispersion of all particles;M02 + M20 (cm^{2});#it{E}_{clu} (GeV)",
        600, 0, 24, ebins, emin, emax));
    // fOutputContainer->Add(new TH2F("Nbar_M02_2", "Dispersion of all
    // particles;M02 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax)); fOutputContainer->Add(new TH2F("Nbar_M20_2", "Dispersion of all
    // particles;M20 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax));
    fOutputContainer->Add(new TH2F(
        "Nbar_Sum_2",
        "Dispersion of all particles;M02 + M20 (cm^{2});#it{E}_{clu} (GeV)",
        600, 0, 24, ebins, emin, emax));
    // fOutputContainer->Add(new TH2F("Pbar_M02_1", "Dispersion of all
    // particles;M02 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax)); fOutputContainer->Add(new TH2F("Pbar_M20_1", "Dispersion of all
    // particles;M20 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax)); fOutputContainer->Add(new TH2F("Pbar_Sum_1", "Dispersion of all
    // particles;M02 + M20 (cm^{2});#it{E}_{clu} (GeV)", 600, 0, 24, ebins,
    // emin, emax)); fOutputContainer->Add(new TH2F("Pbar_M02_2", "Dispersion of
    // all particles;M02 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax)); fOutputContainer->Add(new TH2F("Pbar_M20_2", "Dispersion of all
    // particles;M20 (cm^{2});#it{E}_{clu} (GeV)", 300, 0, 12, ebins, emin,
    // emax)); fOutputContainer->Add(new TH2F("Pbar_Sum_2", "Dispersion of all
    // particles;M02 + M20 (cm^{2});#it{E}_{clu} (GeV)", 600, 0, 24, ebins,
    // emin, emax));
  }

  // Invmas for MC truth AntiSigma with selections
  if (fIsMC) {
    fOutputContainer->Add(new TH2F("Cuts1_InvMass_AntiSigmaPlus_ParentCheck",
                                   "Invariant mass", invbins, invmin, invmax,
                                   ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Cuts1_InvMass_AntiSigmaMinus_ParentCheck",
                                   "Invariant mass", invbins, invmin, invmax,
                                   ebins, emin, emax));
    if (fQAhist) {
      fOutputContainer->Add(new TH2F("Ncell_InvMass_AntiSigmaPlus_ParentCheck",
                                     "Invariant mass", invbins, invmin, invmax,
                                     ebins, emin, emax));
      fOutputContainer->Add(new TH2F("Ncell_InvMass_AntiSigmaMinus_ParentCheck",
                                     "Invariant mass", invbins, invmin, invmax,
                                     ebins, emin, emax));
      fOutputContainer->Add(new TH2F("Disp_InvMass_AntiSigmaPlus_ParentCheck",
                                     "Invariant mass", invbins, invmin, invmax,
                                     ebins, emin, emax));
      fOutputContainer->Add(new TH2F("Disp_InvMass_AntiSigmaMinus_ParentCheck",
                                     "Invariant mass", invbins, invmin, invmax,
                                     ebins, emin, emax));
      fOutputContainer->Add(new TH2F("Wo_InvMass_AntiSigmaPlus_ParentCheck",
                                     "Invariant mass", invbins, invmin, invmax,
                                     ebins, emin, emax));
      fOutputContainer->Add(new TH2F("Wo_InvMass_AntiSigmaMinus_ParentCheck",
                                     "Invariant mass", invbins, invmin, invmax,
                                     ebins, emin, emax));
    }
    if (fInvMassHist) {
      for (Int_t i = 1; i < 34; i++) {
        fOutputContainer->Add(new TH2F(
            Form("Cuts%d_InvMass_AntiSigmaPlus_ParentCheck", i + 1),
            "Invariant mass", invbins, invmin, invmax, ebins, emin, emax));
        fOutputContainer->Add(new TH2F(
            Form("Cuts%d_InvMass_AntiSigmaMinus_ParentCheck", i + 1),
            "Invariant mass", invbins, invmin, invmax, ebins, emin, emax));
      }
    }
    // fOutputContainer->Add(new TH2F("Cuts1_InvMass_AntiSigmaPlus_Parent",
    // "Invariant mass", invbins, invmin, invmax, ebins, emin, emax));
    // fOutputContainer->Add(new TH2F("Cuts1_InvMass_AntiSigmaMinus_Parent",
    // "Invariant mass", invbins, invmin, invmax, ebins, emin, emax));
    // fOutputContainer->Add(new TH1F("Plus_Parent_PDG", "pdg", 8000, -4000,
    // 4000)); fOutputContainer->Add(new TH1F("Minus_Parent_PDG", "pdg", 8000,
    // -4000, 4000)); fOutputContainer->Add(new TH1F("Plus_Cluster_PDG", "pdg",
    // 8000, -4000, 4000)); fOutputContainer->Add(new TH1F("Minus_Cluster_PDG",
    // "pdg", 8000, -4000, 4000));
  }

  if (fAdditionHist) {
    fOutputContainer->Add(new TH2F("DCA_Daughters_Plus",
                                   "DCA between daughters", dcabins, dcamin,
                                   dcamax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("DCA_Daughters_Test",
                                   "DCA between daughters", dcabins, dcamin,
                                   dcamax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("DCA_Daughters_Minus",
                                   "DCA between daughters", dcabins, dcamin,
                                   dcamax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("CPA_Daughters_Plus",
                                   "DCA between daughters", cpabins, cpamin,
                                   cpamax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("CPA_Daughters_Test",
                                   "DCA between daughters", cpabins, cpamin,
                                   cpamax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("CPA_Daughters_Minus",
                                   "DCA between daughters", cpabins, cpamin,
                                   cpamax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("RAD_Daughters_Plus",
                                   "RAD between daughters", radbins, radmin,
                                   radmax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("RAD_Daughters_Test",
                                   "RAD between daughters", radbins, radmin,
                                   radmax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("RAD_Daughters_Minus",
                                   "RAD between daughters", radbins, radmin,
                                   radmax, ebins, emin, emax));

    fOutputContainer->Add(
        new TH2F("DCA_XY_pion", "DCA XY", 480, -2.4, 2.4, ebins, emin, emax));
    fOutputContainer->Add(
        new TH2F("DCA_Z_pion", "DCA XY", 640, -3.2, 3.2, ebins, emin, emax));
  }

  // Inv mass
  fOutputContainer->Add(new TH2F("Cuts1_Mixed_InvMass_Charge1",
                                 "Invariant mass", invbins, invmin, invmax,
                                 ebins, emin, emax));
  fOutputContainer->Add(new TH2F("Cuts1_Mixed_InvMass_Charge-1",
                                 "Invariant mass", invbins, invmin, invmax,
                                 ebins, emin, emax));
  fOutputContainer->Add(new TH2F("Cuts1_InvMass_Charge1", "Invariant mass",
                                 invbins, invmin, invmax, ebins, emin, emax));
  fOutputContainer->Add(new TH2F("Cuts1_InvMass_Charge-1", "Invariant mass",
                                 invbins, invmin, invmax, ebins, emin, emax));
  if (fQAhist) {
    fOutputContainer->Add(new TH2F("Ncell_Mixed_InvMass_Charge1",
                                   "Invariant mass", invbins, invmin, invmax,
                                   ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Ncell_Mixed_InvMass_Charge-1",
                                   "Invariant mass", invbins, invmin, invmax,
                                   ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Disp_Mixed_InvMass_Charge1",
                                   "Invariant mass", invbins, invmin, invmax,
                                   ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Disp_Mixed_InvMass_Charge-1",
                                   "Invariant mass", invbins, invmin, invmax,
                                   ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Wo_Mixed_InvMass_Charge1", "Invariant mass",
                                   invbins, invmin, invmax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Wo_Mixed_InvMass_Charge-1",
                                   "Invariant mass", invbins, invmin, invmax,
                                   ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Cuts1_Mixed_Pi0", "Invariant mass", invbins,
                                   0., 0.5, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Ncell_InvMass_Charge1", "Invariant mass",
                                   invbins, invmin, invmax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Ncell_InvMass_Charge-1", "Invariant mass",
                                   invbins, invmin, invmax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Disp_InvMass_Charge1", "Invariant mass",
                                   invbins, invmin, invmax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Disp_InvMass_Charge-1", "Invariant mass",
                                   invbins, invmin, invmax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Wo_InvMass_Charge1", "Invariant mass",
                                   invbins, invmin, invmax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Wo_InvMass_Charge-1", "Invariant mass",
                                   invbins, invmin, invmax, ebins, emin, emax));
    fOutputContainer->Add(new TH2F("Cuts1_Pi0", "Invariant mass", invbins, 0.,
                                   0.5, ebins, emin, emax));
  }

  if (fInvMassHist && !fIsMC) {
    for (Int_t i = 1; i < 34; i++) {
      fOutputContainer->Add(new TH2F(
          Form("Cuts%d_Mixed_InvMass_Charge1", i + 1), "Invariant mass",
          invbins, invmin, invmax, ebins, emin, emax));
      fOutputContainer->Add(new TH2F(
          Form("Cuts%d_Mixed_InvMass_Charge-1", i + 1), "Invariant mass",
          invbins, invmin, invmax, ebins, emin, emax));

      fOutputContainer->Add(new TH2F(Form("Cuts%d_InvMass_Charge1", i + 1),
                                     "Invariant mass", invbins, invmin, invmax,
                                     ebins, emin, emax));
      fOutputContainer->Add(new TH2F(Form("Cuts%d_InvMass_Charge-1", i + 1),
                                     "Invariant mass", invbins, invmin, invmax,
                                     ebins, emin, emax));
    }
  }

  //========================================//
  fOutputContainer->Add(new TH1F("hSelEvents", "Events selected", 14, 0., 14));
  fOutputContainer->Add(new TH1F("hCentralityV0M", "V0M", 105, 0, 105));
  fOutputContainer->Add(new TH1F("hZvertex", "Z vertex", 200, -50., +50.));
  fOutputContainer->Add(new TH1F("hNvertexTracks",
                                 "N of primary tracks from the primary vertex",
                                 150, 0., 150.));

  for (Int_t i = 0; i < 10; i++) {
    for (Int_t j = 0; j < 10; j++) {
      fPHOSEvents[i][j] = 0x0; // Container for PHOS photons
    }
  }
  for (Int_t mod = 0; mod < 6; mod++)
    fPHOSBadMap[mod] = 0x0;

  PostData(1, fOutputContainer);
}

void AliAnalysisSigmaBarCharged::UserExec(Option_t *) {
  // analyze one event
  FillHistogram("hSelEvents", 1);
  fEvent = dynamic_cast<AliAODEvent *>(InputEvent());
  if (!fEvent) {
    Printf("ERROR: Could not retrieve event");
    // PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents", 2);

  fRunNumber = fEvent->GetRunNumber();

  // read geometry if not read yet
  if (fPHOSgeom == 0) {
    InitGeometry();
  }
  if (!fPHOSgeom) {
    Printf("ERROR: PHOS geometry hasn't been connected!");
  }

  fUtils = new AliAnalysisUtils();
  if (!fUtils) {
    Printf("ERROR: Could not retrieve an. utils");
    // PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents", 3);

  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler *>(
      AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  // particle identification
  if (!fPIDResponse) {
    if (aodH) {
      fPIDResponse = aodH->GetPIDResponse();
    } else {
      Printf("ERROR: Could not get AODInputHandler");
      return;
    }
  }
  FillHistogram("hSelEvents", 4);

  // Minimum Bias events
  if (aodH && !fIsMC) {
    if (!(aodH->IsEventSelected() & AliVEvent::kINT7)) {
      return;
    }
  }
  FillHistogram("hSelEvents", 5);

  // Test vertex
  const AliAODVertex *esdVertex5 = fEvent->GetPrimaryVertex();
  if (!esdVertex5) {
    Printf("ERROR: No primary vertex in AOD!");
    return;
  }
  fvtx5[0] = esdVertex5->GetX();
  fvtx5[1] = esdVertex5->GetY();
  fvtx5[2] = esdVertex5->GetZ();

  FillHistogram("hNvertexTracks", esdVertex5->GetNContributors());
  if (esdVertex5->GetNContributors() <= 0) {
    // PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents", 6);

  FillHistogram("hZvertex", fvtx5[2]);
  if (TMath::Abs(fvtx5[2]) > 10.) {
    // PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents", 7);

  // Pileup selection
  if (fUtils->IsPileUpEvent(fEvent))
    return;
  FillHistogram("hSelEvents", 8);

  // Pileup in Monte-Carlo
  if (fIsMC) {
    AliAODMCHeader *aodMCHeader = (AliAODMCHeader *)fEvent->FindListObject(
        AliAODMCHeader::StdBranchName());
    if (aodMCHeader) {
      // find cocktail header
      Int_t nGenerators = aodMCHeader->GetNCocktailHeaders();
      if (nGenerators > 0) {
        for (Int_t igen = 0; igen < nGenerators; igen++) {
          AliGenEventHeader *eventHeaderGen =
              aodMCHeader->GetCocktailHeader(igen);
          TString genname = eventHeaderGen->ClassName();
          bool isPileUp =
              AliAnalysisUtils::IsPileupInGeneratedEvent(aodMCHeader, genname);
          if (isPileUp)
            return;
        }
      }
    }
  }
  FillHistogram("hSelEvents", 9);

  // Centrality
  AliMultSelection *multSelection =
      (AliMultSelection *)fEvent->FindListObject("MultSelection");
  if (!multSelection) {
    Printf("ERROR: No Cent. inf.!");
    return;
  }
  FillHistogram("hSelEvents", 10);

  fCentrality = multSelection->GetMultiplicityPercentile("V0M");
  // Centrality class bin
  fCentBin = int(fCentrality / 10.);
  if (fCentBin >= 10) {
    fCentBin = 9;
  }
  FillHistogram("hCentralityV0M", fCentrality);

  // Vtx class z-bin
  Int_t zvtx = (Int_t)((fvtx5[2] + 10.) / 2.);
  if (zvtx < 0)
    zvtx = 0;
  if (zvtx > 9)
    zvtx = 9;

  if (!fPHOSEvents[zvtx][fCentBin]) {
    fPHOSEvents[zvtx][fCentBin] = new TList();
  }
  fCurrentMixedList = fPHOSEvents[zvtx][fCentBin];

  SelectSigma();

  // Remove old events
  if (fGamma->GetEntriesFast() > 0) {
    fCurrentMixedList->AddFirst(fGamma);
    fGamma = 0x0;
    if (fCurrentMixedList->GetSize() > 100) {
      TClonesArray *tmp =
          static_cast<TClonesArray *>(fCurrentMixedList->Last());
      fCurrentMixedList->Remove(tmp);
      delete tmp;
    }
  }

  PostData(1, fOutputContainer);
}

void AliAnalysisSigmaBarCharged::SelectSigma() {
  const Double_t c = 29979245800.;    // speed of light in cm/sec
  const Double_t mbar = 0.939485;     // neutron mass
  const Double_t mpi = 0.13957039;    // pion mass
  const Double_t msigmap = 1.197449;  //
  const Double_t msigmam = 1.18937;   //
  const Double_t ltantip = 1.479e-10; // lifetime of Sigma-plus-bar
  const Double_t ltantim = 8.018e-11; // lifetime of Sigma-minus-bar
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
    fStack = (TClonesArray *)fEvent->FindListObject(
        AliAODMCParticle::StdBranchName());

    // Generated particles
    Int_t multMC = fStack->GetEntriesFast();
    for (Int_t i = 0; i < multMC; i++) {
      AliAODMCParticle *prim = (AliAODMCParticle *)fStack->At(i);
      if (prim->GetPdgCode() == -3222 && TMath::Abs(prim->Y()) < 0.5) {
        if (prim->GetDaughterFirst() == -1 || prim->GetDaughterLast() == -1)
          continue;
        AliAODMCParticle *daug1 =
            (AliAODMCParticle *)fStack->At(prim->GetDaughterFirst());
        AliAODMCParticle *daug2 =
            (AliAODMCParticle *)fStack->At(prim->GetDaughterLast());
        if ((daug1->GetPdgCode() == -2112 && daug2->GetPdgCode() == -211) ||
            (daug2->GetPdgCode() == -2112 && daug1->GetPdgCode() == -211)) {
          FillHistogram("MC_AntiSigmaMinus_1", prim->Pt());
        }
      }
      if (prim->GetPdgCode() == -3112 && TMath::Abs(prim->Y()) < 0.5) {
        if (prim->GetDaughterFirst() == -1 || prim->GetDaughterLast() == -1)
          continue;
        AliAODMCParticle *daug1 =
            (AliAODMCParticle *)fStack->At(prim->GetDaughterFirst());
        AliAODMCParticle *daug2 =
            (AliAODMCParticle *)fStack->At(prim->GetDaughterLast());
        if ((daug1->GetPdgCode() == -2112 && daug2->GetPdgCode() == 211) ||
            (daug2->GetPdgCode() == -2112 && daug1->GetPdgCode() == 211)) {
          FillHistogram("MC_AntiSigmaPlus_1", prim->Pt());
        }
      }
      if (fAdditionHist) {
        if (prim->GetPdgCode() == -2112 && TMath::Abs(prim->Y()) < 0.5) {
          FillHistogram("MC_Nbar_1", prim->Pt());
          if (prim->IsPhysicalPrimary() || prim->IsPrimary()) {
            FillHistogram("MC_Nbar_prim_1", prim->Pt());
          } else {
            FillHistogram("MC_Nbar_nonprim_1", prim->Pt());
          }
        }
        if (prim->GetPdgCode() == -2212 && TMath::Abs(prim->Y()) < 0.5) {
          FillHistogram("MC_Pbar_1", prim->Pt());
        }
      }
      if (fQAhist) {
        if ((prim->GetPdgCode() == 211 || prim->GetPdgCode() == -211) &&
            TMath::Abs(prim->Y()) < 0.5) {
          FillHistogram("MC_Pion_1", prim->Pt());
        }
      }
    }
  } // End of MC particles loop

  // Scan clusters
  for (Int_t i = 0; i < multClust; i++) {
    AliAODCaloCluster *clu = fEvent->GetCaloCluster(i);
    // Select cluster type
    if (clu->GetType() != AliVCluster::kPHOSNeutral)
      continue;
    Double_t cluE = clu->E();
    // Energy shift, see pi0 measurement in pPb@8TeV by D.Peresunko, additional
    // energy non-linearity in MC
    if (fIsMC) {
      cluE = 1.0245 * cluE * (1. - 0.013 * exp(-0.5 * cluE));
    }
    // Exotic or single cell cluster
    if (clu->GetM02() <= 0.2)
      continue;
    if (clu->GetNCells() < 3)
      continue;
    if (cluE < 0.3) { // Min energy
      continue;
    }

    // Reconstructed (photon) momentum
    TLorentzVector lvclu;
    clu->GetMomentum(lvclu, fvtx5);
    // Position in PHOS surface
    Float_t pos[3];
    clu->GetPosition(pos);

    TVector3 global1(pos);
    Int_t relId[4];
    fPHOSgeom->GlobalPos2RelId(global1, relId);
    Int_t mod = relId[0];
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3];

    Double_t t = clu->GetTOF();
    Double_t r = TMath::Sqrt((fvtx5[0] - pos[0]) * (fvtx5[0] - pos[0]) +
                             (fvtx5[1] - pos[1]) * (fvtx5[1] - pos[1]) +
                             (fvtx5[2] - pos[2]) * (fvtx5[2] - pos[2]));
    Double_t tgamma = r / c;
    // Real data calibrated wrt photon arrival
    if (!fIsMC) {
      t = t + tgamma;
    }
    // Time resolution (realistic)
    Double_t sigt = 0.;
    if (fIsMC && fPHOSClusterTOFOption == 1) {
      sigt = RealRes(cluE);
      t = gRandom->Gaus(t, sigt);
    } // append realistic resolution

    // Time resolution (realistic)
    if (fIsMC && fPHOSClusterTOFOption == 2) {
      t = t + RealResv2(cluE);
    } // append realistic resolution

    // TOF versus E_clu distribution
    FillHistogram("hClusterTOFvsE", t - tgamma, cluE);
    if (fQAhist) {
      FillHistogram(Form("OccupancyPlot_%dM", mod), cellX, cellZ);
      FillHistogram(Form("ClusterTOFvsE_%dM", mod), t - tgamma, cluE);
    }

    if (fAdditionHist) {
      if ((clu->GetM02() >= -1 * clu->GetM20() + fDispCut)) {
        FillHistogram("hClusterTOFvsE_cut1", t - tgamma, cluE);
      }
      if ((clu->GetM02() >= -1 * clu->GetM20() + fDispCut) &&
          clu->GetNCells() >= fNcellCut) {
        FillHistogram("hClusterTOFvsE_cut2", t - tgamma, cluE);
      }
    }

    Int_t primLabel0;
    if (fIsMC) {
      primLabel0 = clu->GetLabelAt(0);
      if (primLabel0 > -1 && fAdditionHist) {
        AliAODMCParticle *prim0 = (AliAODMCParticle *)fStack->At(primLabel0);
        if (prim0->GetPdgCode() == -2112) {
          FillHistogram("hClusterTOFvsE_nbar", t - tgamma, cluE);
        }
      }
    }

    // Time cut; exclude negative times
    if ((t - tgamma <= 0) || (t - tgamma >= fCluTimeCut)) {
      continue;
    }

    // Reconstructed momentum
    Double_t prec = mbar / sqrt(pow(t * c / r, 2) - 1);

    TVector3 recon;
    recon.SetMagThetaPhi(prec, lvclu.Theta(), lvclu.Phi());
    Double_t reconE = sqrt(mbar * mbar + prec * prec);
    TLorentzVector nbar(recon.Px(), recon.Py(), recon.Pz(), reconE);

    // Arrays for clu cut variations
    Double_t Eclucut[13] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                            1.,  1.1, 1.2, 1.3, 1.4, 1.5};
    Double_t Ncellcut[13] = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    Double_t Dispcut[13] = {1.,  1.5, 2.,  2.5, 3.,  3.5, 4.,
                            4.5, 5.,  5.5, 6.,  6.5, 7.0};
    Double_t CPVcut[13] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};

    FillHistogram("hClusterTOFvsE_1", t - tgamma, cluE);
    FillHistogram("hClusterpTvsE_1", nbar.Pt(), cluE);
    if (fAdditionHist) {
      // FillHistogram("Spectrum", cluE);
      FillHistogram("All_Disp", clu->GetM20(), clu->GetM02());
      // FillHistogram("All_M02", clu->GetM02(),cluE);
      // FillHistogram("All_M20", clu->GetM20(),cluE);
      FillHistogram("All_Sum", clu->GetM02() + clu->GetM20(), cluE);
    }

    if (fAdditionHist) {
      // Eclu variation
      if (fIsMC && (clu->GetNCells() >= fNcellCut) &&
          (clu->GetM02() >= -1 * clu->GetM20() + fDispCut) &&
          (clu->GetEmcCpvDistance() > fCPVCut)) {
        FillHistogram("Spectrum_d_0_0_0", cluE);
        for (Int_t i1 = 1; i1 < 14; i1++) {
          if (cluE >= Eclucut[i1 - 1]) {
            FillHistogram(Form("Spectrum_%d_0_0_0", i1), cluE);
          }
        }
      }
      // Ncells variation
      if (fIsMC && (cluE >= fCluNbarMinE) &&
          (clu->GetM02() >= -1 * clu->GetM20() + fDispCut) &&
          (clu->GetEmcCpvDistance() > fCPVCut)) {
        FillHistogram("Spectrum_0_d_0_0", cluE);
        for (Int_t i1 = 1; i1 < 14; i1++) {
          if (clu->GetNCells() >= Ncellcut[i1 - 1]) {
            FillHistogram(Form("Spectrum_0_%d_0_0", i1), cluE);
          }
        }
      }
      // Disp variation
      if (fIsMC && (cluE >= fCluNbarMinE) && (clu->GetNCells() >= fNcellCut) &&
          (clu->GetEmcCpvDistance() > fCPVCut)) {
        FillHistogram("Spectrum_0_0_d_0", cluE);
        for (Int_t i1 = 1; i1 < 14; i1++) {
          if (clu->GetM02() >= -1 * clu->GetM20() + Dispcut[i1 - 1]) {
            FillHistogram(Form("Spectrum_0_0_%d_0", i1), cluE);
          }
        }
      }
      // CPV variation
      if (fIsMC && (cluE >= fCluNbarMinE) && (clu->GetNCells() >= fNcellCut) &&
          (clu->GetM02() >= -1 * clu->GetM20() + fDispCut)) {
        FillHistogram("Spectrum_0_0_0_d", cluE);
        for (Int_t i1 = 1; i1 < 14; i1++) {
          if (clu->GetEmcCpvDistance() > CPVcut[i1 - 1]) {
            FillHistogram(Form("Spectrum_0_0_0_%d", i1), cluE);
          }
        }
      }
    }

    if (fIsMC && fAdditionHist) {
      if (primLabel0 > -1) {
        AliAODMCParticle *prim0 = (AliAODMCParticle *)fStack->At(primLabel0);
        if (prim0->GetPdgCode() == -2112) {
          Double_t pmc = prim0->P();
          FillHistogram("hClusterTOFvsE_nbar_1", t - tgamma, cluE);
          FillHistogram("hClusterpTvsE_nbar_1", nbar.Pt(), cluE);
          FillHistogram("hClusterpTmcvsE_nbar_1", prim0->Pt(), cluE);
          // FillHistogram("Spectrum_Nbar", cluE);
          FillHistogram("AntiNeutron_Disp", clu->GetM20(), clu->GetM02());
          FillHistogram("EcluvsP_nbar", cluE, pmc);
          // Eclu variation
          if ((clu->GetNCells() >= fNcellCut) &&
              (clu->GetM02() >= -1 * clu->GetM20() + fDispCut) &&
              (clu->GetEmcCpvDistance() > fCPVCut)) {
            FillHistogram("Spectrum_Nbar_d_0_0_0", cluE);
            for (Int_t i1 = 1; i1 < 14; i1++) {
              if (cluE >= Eclucut[i1 - 1]) {
                FillHistogram(Form("Spectrum_Nbar_%d_0_0_0", i1), cluE);
              }
            }
          }
          // Ncells variation
          if ((cluE >= fCluNbarMinE) &&
              (clu->GetM02() >= -1 * clu->GetM20() + fDispCut) &&
              (clu->GetEmcCpvDistance() > fCPVCut)) {
            FillHistogram("Spectrum_Nbar_0_d_0_0", cluE);
            for (Int_t i1 = 1; i1 < 14; i1++) {
              if (clu->GetNCells() >= Ncellcut[i1 - 1]) {
                FillHistogram(Form("Spectrum_Nbar_0_%d_0_0", i1), cluE);
              }
            }
          }
          // Disp variation
          if ((cluE >= fCluNbarMinE) && (clu->GetNCells() >= fNcellCut) &&
              (clu->GetEmcCpvDistance() > fCPVCut)) {
            FillHistogram("Spectrum_Nbar_0_0_d_0", cluE);
            for (Int_t i1 = 1; i1 < 14; i1++) {
              if (clu->GetM02() >= -1 * clu->GetM20() + Dispcut[i1 - 1]) {
                FillHistogram(Form("Spectrum_Nbar_0_0_%d_0", i1), cluE);
              }
            }
          }
          // CPV variation
          if ((cluE >= fCluNbarMinE) && (clu->GetNCells() >= fNcellCut) &&
              (clu->GetM02() >= -1 * clu->GetM20() + fDispCut)) {
            FillHistogram("Spectrum_Nbar_0_0_0_d", cluE);
            for (Int_t i1 = 1; i1 < 14; i1++) {
              if (clu->GetEmcCpvDistance() > CPVcut[i1 - 1]) {
                FillHistogram(Form("Spectrum_Nbar_0_0_0_%d", i1), cluE);
              }
            }
          }
        }
        if (prim0->GetPdgCode() == -2212) {
          FillHistogram("EcluvsP_pbar", cluE, prim0->P());
          FillHistogram("AntiProton_Disp", clu->GetM20(), clu->GetM02());
        }
        if (prim0->GetPdgCode() == 22) {
          FillHistogram("Photon_Disp", clu->GetM20(), clu->GetM02());
        }
        // Classify particle
        TString partName;
        switch (prim0->GetPdgCode()) {
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
        case 2112:
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
        }

        // Spectrum of particles in PHOS
        FillHistogram("Spectrum_Sum", cluE);
        FillHistogram("Ncells_Sum", clu->GetNCells());
        FillHistogram(Form("Spectrum_%s", partName.Data()), cluE);
        FillHistogram(Form("Ncells_%s", partName.Data()), clu->GetNCells());
      }
    }

    // Set bits; After nbar cuts

    // Fill the Nbar array
    if (inPHOS >= fGamma->GetSize()) {
      fGamma->Expand(inPHOS + 50);
    }

    // Create AliCaloPhoton
    AliCaloPhoton *ph = new ((*fGamma)[inPHOS++])
        AliCaloPhoton(nbar.Px(), nbar.Py(), nbar.Pz(), nbar.E());

    // Default Cuts
    // Energy cut for antineutrons
    ph->SetCPVBit(cluE >= fCluNbarMinE);
    // Ncells cut
    ph->SetCPV2Bit(clu->GetNCells() >= fNcellCut);
    // CPV cut
    ph->SetDispBit(clu->GetEmcCpvDistance() > fCPVCut);
    // Disp cut
    ph->SetDisp2Bit(clu->GetM02() >= -1 * clu->GetM20() + fDispCut);

    // Record first label of the cluster
    ph->SetPrimaryAtVertex(primLabel0);

    // Set variables for momentum reconstruction
    ph->SetEMCx(pos[0]);
    ph->SetEMCy(pos[1]);
    ph->SetEMCz(pos[2]);
    ph->SetTime(t - tgamma);
    ph->SetWeight(cluE);

    if (fAdditionHist) {
      // Compare spec. in MC and data
      FillHistogram("Spectrum_pt_allnbar_1", ph->Pt());
      FillHistogram("Spectrum_eclu_allnbar_1", cluE);
      // FillHistogram("Nbar_all_M02_1", clu->GetM02(),cluE);
      // FillHistogram("Nbar_all_M20_1", clu->GetM20(),cluE);
      FillHistogram("Nbar_all_Sum_1", clu->GetM02() + clu->GetM20(), cluE);

      // All cut
      if (ph->IsCPVOK() && ph->IsCPV2OK() && ph->IsDispOK() &&
          ph->IsDisp2OK()) {
        FillHistogram("hClusterTOFvsE_2", t - tgamma, cluE);
        FillHistogram("hClusterpTvsE_2", nbar.Pt(), cluE);
        FillHistogram("Spectrum_pt_allnbar_2", ph->Pt());
        FillHistogram("Spectrum_eclu_allnbar_2", cluE);
        // FillHistogram("Nbar_all_M02_2", clu->GetM02(),cluE);
        // FillHistogram("Nbar_all_M20_2", clu->GetM20(),cluE);
        FillHistogram("Nbar_all_Sum_2", clu->GetM02() + clu->GetM20(), cluE);
      }

      if (fIsMC) {
        Int_t primLabel = clu->GetLabelAt(0);
        if (primLabel > -1) {
          AliAODMCParticle *prim = (AliAODMCParticle *)fStack->At(primLabel);
          if (prim->GetPdgCode() == -2112) {
            // Compare spec. in MC and data
            FillHistogram("Spectrum_pt_nbar_1", ph->Pt());
            FillHistogram("Spectrum_ptmc_nbar_1", prim->Pt());
            FillHistogram("Response_nbar_1", ph->Pt(), prim->Pt());
            FillHistogram("Spectrum_eclu_nbar_1", cluE);
            // FillHistogram("Nbar_M02_1", clu->GetM02(),cluE);
            // FillHistogram("Nbar_M20_1", clu->GetM20(),cluE);
            FillHistogram("Nbar_Sum_1", clu->GetM02() + clu->GetM20(), cluE);
            // All cut
            if (ph->IsCPVOK() && ph->IsCPV2OK() && ph->IsDispOK() &&
                ph->IsDisp2OK()) {
              FillHistogram("Spectrum_pt_nbar_2", ph->Pt());
              FillHistogram("Spectrum_ptmc_nbar_2", prim->Pt());
              FillHistogram("Response_nbar_2", ph->Pt(), prim->Pt());
              FillHistogram("Spectrum_eclu_nbar_2", cluE);
              // FillHistogram("Nbar_M02_2", clu->GetM02(),cluE);
              // FillHistogram("Nbar_M20_2", clu->GetM20(),cluE);
              FillHistogram("Nbar_Sum_2", clu->GetM02() + clu->GetM20(), cluE);
            }
          }
        }
      }

      if (fPIDResponse && clu->GetNTracksMatched() > 0) {
        AliAODTrack *trackclu =
            dynamic_cast<AliAODTrack *>(clu->GetTrackMatched(0));
        Double_t nsigmaProtonTPC =
            fPIDResponse->NumberOfSigmasTPC(trackclu, AliPID::kProton);
        FillHistogram("TPC_pid_proton", nsigmaProtonTPC, trackclu->P());
        if (TMath::Abs(nsigmaProtonTPC) < 3. && trackclu->Charge() == -1) {
          // Compare spec. in MC and data
          FillHistogram("Spectrum_pt_allpbar_1", ph->Pt());
          FillHistogram("Spectrum_eclu_allpbar_1", cluE);
          // FillHistogram("Pbar_all_M02_1", clu->GetM02(),cluE);
          // FillHistogram("Pbar_all_M20_1", clu->GetM20(),cluE);
          // FillHistogram("Pbar_all_Sum_1", clu->GetM02()+clu->GetM20(),cluE);

          // All cut
          if (ph->IsCPVOK() && ph->IsCPV2OK() && ph->IsDisp2OK() &&
              clu->GetEmcCpvDistance() < 2.) {
            FillHistogram("Spectrum_pt_allpbar_2", ph->Pt());
            FillHistogram("Spectrum_eclu_allpbar_2", cluE);
            // FillHistogram("Pbar_all_M02_2", clu->GetM02(),cluE);
            // FillHistogram("Pbar_all_M20_2", clu->GetM20(),cluE);
            // FillHistogram("Pbar_all_Sum_2",
            // clu->GetM02()+clu->GetM20(),cluE);
          }
          if (fIsMC) {
            Int_t primLabel = clu->GetLabelAt(0);
            if (primLabel > -1) {
              AliAODMCParticle *prim =
                  (AliAODMCParticle *)fStack->At(primLabel);
              if (prim->GetPdgCode() == -2212) {
                // Compare spec. in MC and data
                FillHistogram("Spectrum_pt_pbar_1", ph->Pt());
                FillHistogram("Spectrum_ptmc_pbar_1", prim->Pt());
                FillHistogram("Response_pbar_1", ph->Pt(), prim->Pt());
                FillHistogram("Spectrum_eclu_pbar_1", cluE);
                // FillHistogram("Pbar_M02_1", clu->GetM02(),cluE);
                // FillHistogram("Pbar_M20_1", clu->GetM20(),cluE);
                // FillHistogram("Pbar_Sum_1",
                // clu->GetM02()+clu->GetM20(),cluE); All cut
                if (ph->IsCPVOK() && ph->IsCPV2OK() && ph->IsDisp2OK() &&
                    clu->GetEmcCpvDistance() < 2.) {
                  FillHistogram("Spectrum_pt_pbar_2", ph->Pt());
                  FillHistogram("Spectrum_ptmc_pbar_2", prim->Pt());
                  FillHistogram("Response_pbar_2", ph->Pt(), prim->Pt());
                  FillHistogram("Spectrum_eclu_pbar_2", cluE);
                  // FillHistogram("Pbar_M02_2", clu->GetM02(),cluE);
                  // FillHistogram("Pbar_M20_2", clu->GetM20(),cluE);
                  // FillHistogram("Pbar_Sum_2",
                  // clu->GetM02()+clu->GetM20(),cluE);
                }
              }
            }
          }
        }
      }
    }

    if (fIsMC && ph->IsCPVOK() && ph->IsCPV2OK() && ph->IsDispOK() &&
        ph->IsDisp2OK() && fAdditionHist) {
      Int_t primLabel = clu->GetLabelAt(0);
      if (primLabel > -1) {
        AliAODMCParticle *prim = (AliAODMCParticle *)fStack->At(primLabel);
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
        case 2112:
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
          // FillHistogram("PDG_Other", prim->GetPdgCode());
        }

        // Spectrum of particles in PHOS
        FillHistogram("Spectrum_Sum_1", cluE);
        FillHistogram("Ncells_Sum_1", clu->GetNCells());
        FillHistogram(Form("Spectrum_%s_1", partName.Data()), cluE);
        FillHistogram(Form("Ncells_%s_1", partName.Data()), clu->GetNCells());

        if (prim->GetPdgCode() == -2112) {
          FillHistogram("hClusterTOFvsE_nbar_2", t - tgamma, cluE);
          FillHistogram("hClusterpTvsE_nbar_2", nbar.Pt(), cluE);
          FillHistogram("hClusterpTmcvsE_nbar_2", prim->Pt(), cluE);

          Double_t vtxmc[3];
          prim->XvYvZv(vtxmc);
          Double_t pmc = prim->P();
          Double_t Emc = prim->E();

          if (prim->IsPhysicalPrimary() || prim->IsPrimary()) {
            FillHistogram("Nbar_pmcvsprec-pmc_prim", pmc, -pmc + prec);
            FillHistogram("Nbar_prec-pmcvseclu_prim", cluE, -pmc + prec);
            FillHistogram("Nbar_pmcvsprec_prim", pmc, prec);
            FillHistogram("Nbar_emcvserec_prim", Emc, reconE);
          } else {
            FillHistogram("Nbar_pmcvsprec-pmc", pmc, -pmc + prec);
            FillHistogram("Nbar_prec-pmcvseclu", cluE, -pmc + prec);
            FillHistogram("Nbar_pmcvsprec", pmc, prec);
            FillHistogram("Nbar_emcvserec", Emc, reconE);
          }

          Int_t primparentLabel = prim->GetMother();
          if (primparentLabel > -1) {
            AliAODMCParticle *primparent =
                (AliAODMCParticle *)fStack->At(primparentLabel);
            if (primparent->GetPdgCode() == -3222) {
              AliAODMCParticle *pion;
              if (primparent->GetDaughterLast() == primLabel) {
                pion = (AliAODMCParticle *)fStack->At(
                    primparent->GetDaughterFirst());
              } else if (primparent->GetDaughterFirst() == primLabel) {
                pion = (AliAODMCParticle *)fStack->At(
                    primparent->GetDaughterLast());
              }
              if (pion->GetPdgCode() == -211) {
                FillHistogram("Nbar_pmcvsprec-pmc_minus", pmc, -pmc + prec);
                FillHistogram("Nbar_prec-pmcvseclu_minus", cluE, -pmc + prec);
                FillHistogram("Nbar_pmcvsprec_minus", pmc, prec);

                Double_t vtxpion[3]; // Vertex of Sigma bar decay
                pion->XvYvZv(vtxpion);
                Double_t rtosecondary = TMath::Sqrt(
                    (vtxpion[0] - fvtx5[0]) * (vtxpion[0] - fvtx5[0]) +
                    (vtxpion[1] - fvtx5[1]) * (vtxpion[1] - fvtx5[1]) +
                    (vtxpion[2] - fvtx5[2]) *
                        (vtxpion[2] - fvtx5[2])); // PV to decay vertex
                Double_t rsectophs = TMath::Sqrt(
                    (vtxpion[0] - pos[0]) * (vtxpion[0] - pos[0]) +
                    (vtxpion[1] - pos[1]) * (vtxpion[1] - pos[1]) +
                    (vtxpion[2] - pos[2]) *
                        (vtxpion[2] - pos[2])); // Decay vertex to PHOS

                // AntiSigma tof compare
                Int_t parentSigmalabel = primparent->GetMother();
                Double_t tofparent = primparent->T();
                while (parentSigmalabel > -1) {
                  AliAODMCParticle *Sigmamother =
                      (AliAODMCParticle *)fStack->At(parentSigmalabel);
                  tofparent = Sigmamother->T();
                  parentSigmalabel = Sigmamother->GetMother();
                }

                Double_t tofsigma = prim->T() - tofparent;
                Double_t tofnbar = rtosecondary / (pmc * c / Emc);
                Double_t prec_corr =
                    mbar / sqrt(pow((t - tofsigma) * c / rsectophs, 2) - 1);

                FillHistogram("t_compare_minus",
                              (tofsigma - tofnbar) / tofsigma,
                              primparent->Pt());
                FillHistogram("t_compare_minus_1", (tofsigma - tofnbar) / t,
                              primparent->Pt());
                FillHistogram("t_compare_minus_2", tofsigma / t,
                              primparent->Pt());
                FillHistogram("SV_minus", (rtosecondary + rsectophs - r) / r,
                              primparent->Pt());

                FillHistogram("t_compare_minus_3", (-prec_corr + prec) / prec,
                              primparent->Pt());
              }
            }
            if (primparent->GetPdgCode() == -3112) {
              AliAODMCParticle *pion;
              if (primparent->GetDaughterLast() == primLabel) {
                pion = (AliAODMCParticle *)fStack->At(
                    primparent->GetDaughterFirst());
              } else if (primparent->GetDaughterFirst() == primLabel) {
                pion = (AliAODMCParticle *)fStack->At(
                    primparent->GetDaughterLast());
              }
              if (pion->GetPdgCode() == 211) {
                FillHistogram("Nbar_pmcvsprec-pmc_plus", pmc, -pmc + prec);
                FillHistogram("Nbar_prec-pmcvseclu_plus", cluE, -pmc + prec);
                FillHistogram("Nbar_pmcvsprec_plus", pmc, prec);

                Double_t vtxpion[3]; // Vertex of Sigma bar decay
                pion->XvYvZv(vtxpion);
                Double_t rtosecondary = TMath::Sqrt(
                    (vtxpion[0] - fvtx5[0]) * (vtxpion[0] - fvtx5[0]) +
                    (vtxpion[1] - fvtx5[1]) * (vtxpion[1] - fvtx5[1]) +
                    (vtxpion[2] - fvtx5[2]) *
                        (vtxpion[2] - fvtx5[2])); // PV to decay vertex
                Double_t rsectophs = TMath::Sqrt(
                    (vtxpion[0] - pos[0]) * (vtxpion[0] - pos[0]) +
                    (vtxpion[1] - pos[1]) * (vtxpion[1] - pos[1]) +
                    (vtxpion[2] - pos[2]) *
                        (vtxpion[2] - pos[2])); // Decay vertex to PHOS

                // AntiSigma tof compare
                Int_t parentSigmalabel = primparent->GetMother();
                Double_t tofparent = primparent->T();
                while (parentSigmalabel > -1) {
                  AliAODMCParticle *Sigmamother =
                      (AliAODMCParticle *)fStack->At(parentSigmalabel);
                  tofparent = Sigmamother->T();
                  parentSigmalabel = Sigmamother->GetMother();
                }

                Double_t tofsigma = prim->T() - tofparent;
                Double_t tofnbar = rtosecondary / (pmc * c / Emc);
                Double_t prec_corr =
                    mbar / sqrt(pow((t - tofsigma) * c / rsectophs, 2) - 1);

                FillHistogram("t_compare_plus", (tofsigma - tofnbar) / tofsigma,
                              primparent->Pt());
                FillHistogram("t_compare_plus_1", (tofsigma - tofnbar) / t,
                              primparent->Pt());
                FillHistogram("t_compare_plus_2", tofsigma / t,
                              primparent->Pt());
                FillHistogram("SV_plus", (rtosecondary + rsectophs - r) / r,
                              primparent->Pt());

                FillHistogram("t_compare_plus_3", (-prec_corr + prec) / prec,
                              primparent->Pt());
              }
            }
          }
        }
      }
    }
  }

  TIter nextEv(fCurrentMixedList);

  if (fQAhist) {
    for (Int_t i = 0; i < inPHOS; i++) {
      AliCaloPhoton *ph1 = static_cast<AliCaloPhoton *>(fGamma->At(i));
      Int_t primLabel1 = ph1->GetPrimaryAtVertex();

      Double_t cluE1 = ph1->GetWeight();

      TVector3 nbar1mom;
      nbar1mom.SetMagThetaPhi(cluE1, ph1->Theta(), ph1->Phi());
      TLorentzVector nbar1(nbar1mom.Px(), nbar1mom.Py(), nbar1mom.Pz(), cluE1);
      //=================Event mixing==============
      while (TClonesArray *event2 = static_cast<TClonesArray *>(nextEv())) {
        Int_t nPhotons2 = event2->GetEntriesFast();
        for (Int_t j = 0; j < nPhotons2; j++) {
          AliCaloPhoton *ph2 = static_cast<AliCaloPhoton *>(event2->At(j));
          Int_t primLabel2 = ph2->GetPrimaryAtVertex();

          Double_t cluE2 = ph2->GetWeight();

          TVector3 nbar2mom;
          nbar2mom.SetMagThetaPhi(cluE2, ph2->Theta(), ph2->Phi());
          TLorentzVector nbar2(nbar2mom.Px(), nbar2mom.Py(), nbar2mom.Pz(),
                               cluE2);

          Double_t m = (nbar1 + nbar2).M();
          Double_t pt = (nbar1 + nbar2).Pt();

          FillHistogram("Cuts1_Mixed_Pi0", m, pt);
        }
      }

      //======================Fill Real==========================
      for (Int_t j = i + 1; j < inPHOS; j++) {
        AliCaloPhoton *ph2 = static_cast<AliCaloPhoton *>(fGamma->At(j));
        Int_t primLabel2 = ph2->GetPrimaryAtVertex();

        Double_t cluE2 = ph2->GetWeight();

        TVector3 nbar2mom;
        nbar2mom.SetMagThetaPhi(cluE2, ph2->Theta(), ph2->Phi());
        TLorentzVector nbar2(nbar2mom.Px(), nbar2mom.Py(), nbar2mom.Pz(),
                             cluE2);

        Double_t m = (nbar1 + nbar2).M();
        Double_t pt = (nbar1 + nbar2).Pt();

        FillHistogram("Cuts1_Pi0", m, pt);
      }
    }
  }

  for (Int_t i = 0; i < multTracks; i++) {
    AliAODTrack *track = (AliAODTrack *)fEvent->GetTrack(i);
    Int_t primLabelTrack = track->GetLabel();
    Bool_t isMCpion = kFALSE;
    Bool_t isMCpionsigma = kFALSE;

    if (fIsMC) {
      if (primLabelTrack > -1) {
        AliAODMCParticle *primpion =
            (AliAODMCParticle *)fStack->At(primLabelTrack);
        if (primpion->GetPdgCode() == -211 || primpion->GetPdgCode() == 211) {
          isMCpion = kTRUE;
          if (primpion->GetMother() > -1) {
            AliAODMCParticle *primpionmother =
                (AliAODMCParticle *)fStack->At(primpion->GetMother());
            if (primpionmother->GetPdgCode() == -3112 ||
                primpionmother->GetPdgCode() == -3222) {
              isMCpionsigma = kTRUE;
            }
          }
        }
      }
    }

    // Check of entries of FB5 before and after FB4 selection
    if (track->TestFilterMask(BIT(5)) && fQAhist) {
      FillHistogram("Spec_track_FB5_0", track->Pt());
    }

    // Select FilterBit
    if (!track->TestFilterMask(BIT(fTracksBits))) {
      continue;
    }

    // kink daughters check
    if (track->GetKinkIndex(0) != 0) {
      std::cout << "label: " << primLabelTrack << std::endl;
      if (primLabelTrack > -1) {
        AliAODMCParticle *primpion =
            (AliAODMCParticle *)fStack->At(primLabelTrack);
        std::cout << "pdg: " << primpion->GetPdgCode() << std::endl;
      }
      continue; // TODO?????
    }

    // Spec of tracks with selected FB
    if (fQAhist) {
      FillHistogram("Spec_track_FB4_1", track->Pt());
      if (track->TestFilterMask(BIT(5))) {
        FillHistogram("Spec_track_FB5_1", track->Pt());
      }
    }

    Int_t TPCClust = track->GetTPCNcls();
    if (TMath::Abs(track->Eta()) <= fTrackEta && fAdditionHist) {
      FillHistogram("TPC_clusters_track", TPCClust, track->Pt());
      if (isMCpionsigma) {
        FillHistogram("TPC_clusters_pionsigma", TPCClust, track->Pt());
      }
    }

    Double_t pxtr = track->Px();
    Double_t pytr = track->Py();
    Double_t pztr = track->Pz();
    Double_t etr = sqrt(mpi * mpi + (pxtr * pxtr + pytr * pytr + pztr * pztr));
    Double_t trackPt = track->Pt();
    // DCA
    Float_t b[2];
    track->GetImpactParameters(b[0], b[1]);
    // if (b[0] == -999 || b[1] == -999)
    //   continue;

    // Track PID in TPC sigmas
    Double_t nsigmaPionTPC =
        TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion));
    Double_t nsigmaPionTPC_noabs =
        fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);

    // TPC dE/dx PID efficiency and purity check; TPCClusters
    if (fAdditionHist && TMath::Abs(track->Eta()) <= fTrackEta &&
        TPCClust >= fNTPCclusters) {
      FillHistogram("TPC_pid_track", nsigmaPionTPC_noabs, track->P());
    }
    if (fQAhist && TMath::Abs(track->Eta()) <= fTrackEta &&
        TPCClust >= fNTPCclusters) {
      FillHistogram("Spec_track_FB4_2", trackPt);
      FillHistogram("Acceptance_track", 1 / trackPt * track->Charge(),
                    track->Phi());
      FillHistogram("PhivsPt_track_1", trackPt, track->Phi());
      if (track->TestFilterMask(BIT(5))) {
        FillHistogram("Spec_track_FB5_2", trackPt);
      }
    }

    // DCA distributions
    if (fAdditionHist) {
      // DCA pion to PV
      if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= fNTPCclusters &&
          nsigmaPionTPC < fTPCsigmas) {
        FillHistogram("DCA_XY_pion", b[0], trackPt);
        FillHistogram("DCA_Z_pion", b[1], trackPt);
        if (isMCpionsigma) {
          FillHistogram("DCA_XY_pionsigma", b[0], trackPt);
          FillHistogram("DCA_Z_pionsigma", b[1], trackPt);
        }
      }
    }

    // QA for tracks; with cuts
    if (fQAhist) {
      if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= fNTPCclusters &&
          nsigmaPionTPC < fTPCsigmas) {
        FillHistogram("PhivsPt_track_2", trackPt, track->Phi());
        FillHistogram("Spec_track_FB4_3", trackPt);
        if (isMCpion) {
          FillHistogram("Spec_track_FB4_4", trackPt);
        }
        if (track->TestFilterMask(BIT(5))) {
          FillHistogram("Spec_track_FB5_3", trackPt);
          if (isMCpion) {
            FillHistogram("Spec_track_FB5_4", trackPt);
          }
        }
      }
    }

    //=================Event mixing==============
    while (TClonesArray *event2 = static_cast<TClonesArray *>(nextEv())) {
      Int_t nPhotons2 = event2->GetEntriesFast();
      for (Int_t j = 0; j < nPhotons2; j++) {
        AliCaloPhoton *ph2 = static_cast<AliCaloPhoton *>(event2->At(j));
        Int_t primLabel2 = ph2->GetPrimaryAtVertex();

        // Double_t t = ph2->GetTime();
        Double_t phoscoord[3] = {ph2->EMCx(), ph2->EMCy(), ph2->EMCz()};

        // Try to construct particle out of track and nbar
        Double_t bm = fEvent->GetMagneticField();

        AliESDtrack btrk((AliVTrack *)track);
        AliExternalTrackParam bt(btrk);

        Double_t decayparams[4];
        PropagateToDCACurvedBachelor(decayparams, ph2, fvtx5, &bt, bm);

        Double_t dca = decayparams[3];

        Double_t decayposition[3] = {decayparams[0], decayparams[1],
                                     decayparams[2]};

        pxtr = bt.Px();
        pytr = bt.Py();
        pztr = bt.Pz();
        etr = sqrt(mpi * mpi + (pxtr * pxtr + pytr * pytr + pztr * pztr));
        TLorentzVector lvtr(pxtr, pytr, pztr, etr);
        trackPt = bt.Pt();

        Double_t rad =
            sqrt((fvtx5[0] - decayparams[0]) * (fvtx5[0] - decayparams[0]) +
                 (fvtx5[1] - decayparams[1]) * (fvtx5[1] - decayparams[1]) +
                 (fvtx5[2] - decayparams[2]) *
                     (fvtx5[2] - decayparams[2])); // dist betw prim and v0

        Double_t sectophos =
            sqrt((phoscoord[0] - fvtx5[0]) * (phoscoord[0] - fvtx5[0]) +
                 (phoscoord[1] - fvtx5[1]) * (phoscoord[1] - fvtx5[1]) +
                 (phoscoord[2] - fvtx5[2]) *
                     (phoscoord[2] - fvtx5[2])); // dist betw prim and phos

        TLorentzVector nbar(ph2->Px(), ph2->Py(), ph2->Pz(), ph2->E());

        Double_t p[3]; // momentum of Sigmabar
        p[0] = nbar.Px() + pxtr;
        p[1] = nbar.Py() + pytr;
        p[2] = nbar.Pz() + pztr;

        Double_t cpa =
            GetCosineOfPointingAngle(p, decayposition, fvtx5[0], fvtx5[1],
                                     fvtx5[2]); // Cosine of Pointing Angle

        Double_t m = (lvtr + nbar).M();
        Double_t pt = (lvtr + nbar).Pt();

        // default
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
            TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
            track->Charge() == 1) {
          FillHistogram("Cuts1_Mixed_InvMass_Charge1", m, pt);
        }
        if (fQAhist) {
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
              TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
              ph2->IsDispOK() && track->Charge() == 1) {
            FillHistogram("Wo_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
              TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && track->Charge() == 1) {
            FillHistogram("Ncell_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
              TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Disp_Mixed_InvMass_Charge1", m, pt);
          }
        }
        if (fInvMassHist && !fIsMC) {
          if (TMath::Abs(track->Eta()) <= 0.7 && TPCClust >= fNTPCclusters &&
              nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
              TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
              track->Charge() == 1) {
            FillHistogram("Cuts2_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= 0.9 && TPCClust >= fNTPCclusters &&
              nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
              TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
              track->Charge() == 1) {
            FillHistogram("Cuts3_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= 50 &&
              nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
              TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
              track->Charge() == 1) {
            FillHistogram("Cuts4_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= 70 &&
              nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
              TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
              track->Charge() == 1) {
            FillHistogram("Cuts5_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < 2.5 &&
              TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
              TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
              track->Charge() == 1) {
            FillHistogram("Cuts6_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < 3.5 &&
              TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
              TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
              track->Charge() == 1) {
            FillHistogram("Cuts7_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0. &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts8_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0. &&
              TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts9_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0. &&
              TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts10_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts11_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts12_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts13_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts14_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts15_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts16_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0. &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts17_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0. &&
              TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts18_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0. &&
              TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts19_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts20_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts21_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts22_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts23_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts24_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0. &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts25_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0. &&
              TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts26_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0. &&
              TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts27_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts28_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts29_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts30_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts31_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts32_Mixed_InvMass_Charge1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
            FillHistogram("Cuts33_Mixed_InvMass_Charge1", m, pt);
          }
          if (ph2->GetTime() < 100.e-09) {
            if (TMath::Abs(track->Eta()) <= fTrackEta &&
                TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
                TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
                ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
                track->Charge() == 1) {
              FillHistogram("Cuts34_Mixed_InvMass_Charge1", m, pt);
            }
          }
        }
        // default
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
            TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
            track->Charge() == -1) {
          FillHistogram("Cuts1_Mixed_InvMass_Charge-1", m, pt);
        }
        if (fQAhist) {
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
              TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
              ph2->IsDispOK() && track->Charge() == -1) {
            FillHistogram("Wo_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
              TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && track->Charge() == -1) {
            FillHistogram("Ncell_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
              TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Disp_Mixed_InvMass_Charge-1", m, pt);
          }
        }
        if (fInvMassHist && !fIsMC) {
          if (TMath::Abs(track->Eta()) <= 0.7 && TPCClust >= fNTPCclusters &&
              nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
              TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
              track->Charge() == -1) {
            FillHistogram("Cuts2_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= 0.9 && TPCClust >= fNTPCclusters &&
              nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
              TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
              track->Charge() == -1) {
            FillHistogram("Cuts3_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= 50 &&
              nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
              TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
              track->Charge() == -1) {
            FillHistogram("Cuts4_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= 70 &&
              nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
              TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
              track->Charge() == -1) {
            FillHistogram("Cuts5_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < 2.5 &&
              TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
              TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
              track->Charge() == -1) {
            FillHistogram("Cuts6_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < 3.5 &&
              TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
              TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
              track->Charge() == -1) {
            FillHistogram("Cuts7_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0. &&
              TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts8_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0. &&
              TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts9_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0. &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts10_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts11_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts12_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts13_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts14_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts15_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts16_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0. &&
              TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts17_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0. &&
              TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts18_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0. &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts19_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts20_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts21_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts22_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts23_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts24_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0. &&
              TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts25_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0. &&
              TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts26_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0. &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts27_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts28_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts29_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts30_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts31_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts32_Mixed_InvMass_Charge-1", m, pt);
          }
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 &&
              TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
              ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
            FillHistogram("Cuts33_Mixed_InvMass_Charge-1", m, pt);
          }
          if (ph2->GetTime() < 100.e-09) {
            if (TMath::Abs(track->Eta()) <= fTrackEta &&
                TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
                TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
                ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
                track->Charge() == -1) {
              FillHistogram("Cuts34_Mixed_InvMass_Charge-1", m, pt);
            }
          }
        }
      }
    }

    //======================Fill Real==========================
    for (Int_t j = 0; j < inPHOS; j++) {
      AliCaloPhoton *ph2 = static_cast<AliCaloPhoton *>(fGamma->At(j));
      Int_t primLabel2 = ph2->GetPrimaryAtVertex();

      // Double_t t = ph2->GetTime();
      Double_t phoscoord[3] = {ph2->EMCx(), ph2->EMCy(), ph2->EMCz()};

      // Try to construct particle out of track and nbar
      Double_t bm = fEvent->GetMagneticField();

      AliESDtrack btrk((AliVTrack *)track);
      AliExternalTrackParam bt(btrk);

      Double_t decayparams[4];
      PropagateToDCACurvedBachelor(decayparams, ph2, fvtx5, &bt, bm);

      Double_t dca = decayparams[3];

      Double_t decayposition[3] = {decayparams[0], decayparams[1],
                                   decayparams[2]};

      pxtr = bt.Px();
      pytr = bt.Py();
      pztr = bt.Pz();
      etr = sqrt(mpi * mpi + (pxtr * pxtr + pytr * pytr + pztr * pztr));
      TLorentzVector lvtr(pxtr, pytr, pztr, etr);
      trackPt = bt.Pt();

      Double_t rad =
          sqrt((fvtx5[0] - decayparams[0]) * (fvtx5[0] - decayparams[0]) +
               (fvtx5[1] - decayparams[1]) * (fvtx5[1] - decayparams[1]) +
               (fvtx5[2] - decayparams[2]) *
                   (fvtx5[2] - decayparams[2])); // dist betw prim and v0

      Double_t sectophos =
          sqrt((phoscoord[0] - fvtx5[0]) * (phoscoord[0] - fvtx5[0]) +
               (phoscoord[1] - fvtx5[1]) * (phoscoord[1] - fvtx5[1]) +
               (phoscoord[2] - fvtx5[2]) *
                   (phoscoord[2] - fvtx5[2])); // dist betw prim and phos

      TLorentzVector nbar(ph2->Px(), ph2->Py(), ph2->Pz(), ph2->E());

      Double_t p[3]; // momentum of Sigmabar
      p[0] = nbar.Px() + pxtr;
      p[1] = nbar.Py() + pytr;
      p[2] = nbar.Pz() + pztr;

      Double_t cpa =
          GetCosineOfPointingAngle(p, decayposition, fvtx5[0], fvtx5[1],
                                   fvtx5[2]); // Cosine of Pointing Angle

      Double_t m = (lvtr + nbar).M();
      Double_t pt = (lvtr + nbar).Pt();

      Bool_t IsSigmaPlusv1 = kFALSE;
      Bool_t IsSigmaPlusv2 = kFALSE;
      Bool_t IsSigmaMinusv1 = kFALSE;
      Bool_t IsSigmaMinusv2 = kFALSE;
      Bool_t IsPairPlus = kFALSE;
      Bool_t IsPairMinus = kFALSE;
      Int_t PdgPairPlus = -1;
      Int_t PdgPairMinus = -1;

      Double_t DCAcut[10] = {0.01, 0.02, 0.03, 0.04, 0.05,
                             0.06, 0.07, 0.08, 0.09, 0.1};
      Double_t CPAcut[10] = {-0.9, -0.7, -0.5, -0.3, 0.,
                             0.3,  0.5,  0.7,  0.9,  0.95};
      Double_t RADcut[10] = {0.05, 0.1,  0.15, 0.2,  0.25,
                             0.3,  0.35, 0.4,  0.45, 0.5};

      if (fIsMC) {
        if (primLabelTrack > -1 && primLabel2 > -1) {
          AliAODMCParticle *primpion =
              (AliAODMCParticle *)fStack->At(primLabelTrack);
          AliAODMCParticle *primanti =
              (AliAODMCParticle *)fStack->At(primLabel2);
          if (primpion->GetMother() > -1 && primanti->GetMother() > -1) {
            AliAODMCParticle *pionparent =
                (AliAODMCParticle *)fStack->At(primpion->GetMother());
            AliAODMCParticle *antiparent =
                (AliAODMCParticle *)fStack->At(primanti->GetMother());
            // IsSigmaPlusv1 - first parent
            // if (pionparent->GetPdgCode() == -3112 && primpion->GetPdgCode()
            // == 211 && track->Charge()==1) {
            //   if (antiparent->GetPdgCode() == -3112 && pionparent ==
            //   antiparent) {
            //     IsSigmaPlusv1 = kTRUE;
            //   } else {
            //     Int_t dummy = antiparent->GetMother();
            //     while (dummy > -1) {
            //       AliAODMCParticle* antigrandparent =
            //       (AliAODMCParticle*)fStack->At(dummy); if
            //       (antigrandparent->GetPdgCode() == -3112 && pionparent ==
            //       antigrandparent) {
            //         IsSigmaPlusv1 = kTRUE;
            //         break;
            //       }
            //       dummy = antigrandparent->GetMother();
            //     }
            //   }
            // }
            // IsSigmaPlusv2 - check all parent
            if (track->Charge() == 1) {
              Int_t dummypion = primpion->GetMother();
              AliAODMCParticle *piongrandparent;
              while (dummypion > -1) {
                piongrandparent = (AliAODMCParticle *)fStack->At(dummypion);
                if (piongrandparent->GetPdgCode() == -3112) {
                  break;
                }
                dummypion = piongrandparent->GetMother();
              }
              if (piongrandparent) {
                Int_t dummyanti = primanti->GetMother();
                while (dummyanti > -1) {
                  AliAODMCParticle *antigrandparent =
                      (AliAODMCParticle *)fStack->At(dummyanti);
                  if (antigrandparent->GetPdgCode() == -3112 &&
                      piongrandparent == antigrandparent) {
                    IsSigmaPlusv2 = kTRUE;
                    break;
                  }
                  dummyanti = antigrandparent->GetMother();
                }
              }
              // Topological selections
              if (fAdditionHist && IsSigmaPlusv2 &&
                  TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
                  ph2->IsDisp2OK()) {
                FillHistogram("DCA_Daughters_AntiSigmaPlus", dca, pt);
                FillHistogram("CPA_Daughters_AntiSigmaPlus", cpa, pt);
                FillHistogram("RAD_Daughters_AntiSigmaPlus", rad, pt);

                // DCAdaug variation
                if (cpa >= fCPAplusCut && TestRAD(pt, rad, fRADplusCut) &&
                    (m >= msigmap - 0.02 && m <= msigmap + 0.02)) {
                  FillHistogram("Sigma_Plus_DCAdaug", pt);
                  for (Int_t i1 = 1; i1 < 11; i1++) {
                    if (TestDCADaug(pt, dca, DCAcut[i1 - 1])) {
                      FillHistogram(Form("Sigma_Plus_DCAdaug_%d", i1), pt);
                    }
                  }
                }
                // RAD variation
                if (cpa >= fCPAplusCut &&
                    TestDCADaug(pt, dca, fDCAdaugplusCut) &&
                    (m >= msigmap - 0.02 && m <= msigmap + 0.02)) {
                  FillHistogram("Sigma_Plus_RAD", pt);
                  for (Int_t i1 = 1; i1 < 11; i1++) {
                    if (TestRAD(pt, rad, RADcut[i1 - 1])) {
                      FillHistogram(Form("Sigma_Plus_RAD_%d", i1), pt);
                    }
                  }
                }
                // CPA variation
                if (TestDCADaug(pt, dca, fDCAdaugplusCut) &&
                    TestRAD(pt, rad, fRADplusCut) &&
                    (m >= msigmap - 0.02 && m <= msigmap + 0.02)) {
                  FillHistogram("Sigma_Plus_CPA", pt);
                  for (Int_t i1 = 1; i1 < 11; i1++) {
                    if (cpa >= CPAcut[i1 - 1]) {
                      FillHistogram(Form("Sigma_Plus_CPA_%d", i1), pt);
                    }
                  }
                }

                Double_t vtxpion[3];
                primpion->XvYvZv(vtxpion);
                // Double_t rsectophs = sqrt((phoscoord[0] - decayparams[0]) *
                // (phoscoord[0] - decayparams[0]) + (phoscoord[1] -
                // decayparams[1]) * (phoscoord[1] - decayparams[1]) +
                //                      (phoscoord[2] - decayparams[2]) *
                //                      (phoscoord[2] - decayparams[2])); //dist
                //                      betw v0 and phos

                FillHistogram("SV_x_plus", vtxpion[0] - decayparams[0],
                              pionparent->Pt());
                FillHistogram("SV_y_plus", vtxpion[1] - decayparams[1],
                              pionparent->Pt());
                FillHistogram("SV_z_plus", vtxpion[2] - decayparams[2],
                              pionparent->Pt());
              }
            }
            // IsSigmaMinusv1 - first parent
            // if (pionparent->GetPdgCode() == -3222 && primpion->GetPdgCode()
            // == -211 && track->Charge()==-1) {
            //   if (antiparent->GetPdgCode() == -3222 && pionparent ==
            //   antiparent) {
            //     IsSigmaMinusv1 = kTRUE;
            //   } else {
            //     Int_t dummy = antiparent->GetMother();
            //     while (dummy > -1) {
            //       AliAODMCParticle* antigrandparent =
            //       (AliAODMCParticle*)fStack->At(dummy); if
            //       (antigrandparent->GetPdgCode() == -3222 && pionparent ==
            //       antigrandparent) {
            //         IsSigmaMinusv1 = kTRUE;
            //         break;
            //       }
            //       dummy = antigrandparent->GetMother();
            //     }
            //   }
            // }
            // IsSigmaMinusv2 - all parent check
            if (track->Charge() == -1) {
              Int_t dummypion = primpion->GetMother();
              AliAODMCParticle *piongrandparent;
              while (dummypion > -1) {
                piongrandparent = (AliAODMCParticle *)fStack->At(dummypion);
                if (piongrandparent->GetPdgCode() == -3222) {
                  break;
                }
                dummypion = piongrandparent->GetMother();
              }
              if (piongrandparent) {
                Int_t dummyanti = primanti->GetMother();
                while (dummyanti > -1) {
                  AliAODMCParticle *antigrandparent =
                      (AliAODMCParticle *)fStack->At(dummyanti);
                  if (antigrandparent->GetPdgCode() == -3222 &&
                      piongrandparent == antigrandparent) {
                    IsSigmaMinusv2 = kTRUE;
                    break;
                  }
                  dummyanti = antigrandparent->GetMother();
                }
              }
              // topological selections
              if (fAdditionHist && IsSigmaMinusv2 &&
                  TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
                  ph2->IsDisp2OK()) {
                FillHistogram("DCA_Daughters_AntiSigmaMinus", dca, pt);
                FillHistogram("CPA_Daughters_AntiSigmaMinus", cpa, pt);
                FillHistogram("RAD_Daughters_AntiSigmaMinus", rad, pt);

                // DCAdaug variation
                if (cpa >= fCPAminusCut && TestRAD(pt, rad, fRADminusCut) &&
                    (m >= msigmam - 0.02 && m <= msigmam + 0.02)) {
                  FillHistogram("Sigma_Minus_DCAdaug", pt);
                  for (Int_t i1 = 1; i1 < 11; i1++) {
                    if (TestDCADaug(pt, dca, DCAcut[i1 - 1])) {
                      FillHistogram(Form("Sigma_Minus_DCAdaug_%d", i1), pt);
                    }
                  }
                }
                // RAD variation
                if (cpa >= fCPAminusCut &&
                    TestDCADaug(pt, dca, fDCAdaugminusCut) &&
                    (m >= msigmam - 0.02 && m <= msigmam + 0.02)) {
                  FillHistogram("Sigma_Minus_RAD", pt);
                  for (Int_t i1 = 1; i1 < 11; i1++) {
                    if (TestRAD(pt, rad, RADcut[i1 - 1])) {
                      FillHistogram(Form("Sigma_Minus_RAD_%d", i1), pt);
                    }
                  }
                }
                // CPA variation
                if (TestDCADaug(pt, dca, fDCAdaugminusCut) &&
                    TestRAD(pt, rad, fRADminusCut) &&
                    (m >= msigmam - 0.02 && m <= msigmam + 0.02)) {
                  FillHistogram("Sigma_Minus_CPA", pt);
                  for (Int_t i1 = 1; i1 < 11; i1++) {
                    if (cpa >= CPAcut[i1 - 1]) {
                      FillHistogram(Form("Sigma_Minus_CPA_%d", i1), pt);
                    }
                  }
                }

                Double_t vtxpion[3];
                primpion->XvYvZv(vtxpion);
                // Double_t rsectophs = sqrt((phoscoord[0] - decayparams[0]) *
                // (phoscoord[0] - decayparams[0]) + (phoscoord[1] -
                // decayparams[1]) * (phoscoord[1] - decayparams[1]) +
                //                      (phoscoord[2] - decayparams[2]) *
                //                      (phoscoord[2] - decayparams[2])); //dist
                //                      betw v0 and phos

                FillHistogram("SV_x_minus", vtxpion[0] - decayparams[0],
                              pionparent->Pt());
                FillHistogram("SV_y_minus", vtxpion[1] - decayparams[1],
                              pionparent->Pt());
                FillHistogram("SV_z_minus", vtxpion[2] - decayparams[2],
                              pionparent->Pt());
              }
            }
            // track and cluster parent, not AntiSigma
            // if (IsSigmaPlusv2 == kFALSE && track->Charge() == 1) {
            //   AliAODMCParticle* piongrandparent;
            //
            //   Int_t dummypion = primpion->GetMother();
            //   while (dummypion > -1) {
            //     piongrandparent = (AliAODMCParticle*)fStack->At(dummypion);
            //     if (piongrandparent->GetPdgCode() != -3112) {
            //       Int_t dummyanti = primanti->GetMother();
            //       while (dummyanti > -1) {
            //         AliAODMCParticle* antigrandparent =
            //         (AliAODMCParticle*)fStack->At(dummyanti); if
            //         (antigrandparent->GetPdgCode() != -3112 &&
            //         piongrandparent == antigrandparent) {
            //           IsPairPlus = kTRUE;
            //           break;
            //         }
            //         dummyanti = antigrandparent->GetMother();
            //       }
            //       if (IsPairPlus) {
            //         break;
            //       }
            //     }
            //     dummypion = piongrandparent->GetMother();
            //   }
            //
            //   if (IsPairPlus) {
            //     PdgPairPlus = piongrandparent->GetPdgCode();
            //   }
            // }
            // track and cluster parent, not AntiSigma
            // if (IsSigmaMinusv2 == kFALSE && track->Charge() == -1) {
            //   AliAODMCParticle* piongrandparent;
            //
            //   Int_t dummypion = primpion->GetMother();
            //   while (dummypion > -1) {
            //     piongrandparent = (AliAODMCParticle*)fStack->At(dummypion);
            //     if (piongrandparent->GetPdgCode() != -3222) {
            //       Int_t dummyanti = primanti->GetMother();
            //       while (dummyanti > -1) {
            //         AliAODMCParticle* antigrandparent =
            //         (AliAODMCParticle*)fStack->At(dummyanti); if
            //         (antigrandparent->GetPdgCode() != -3222 &&
            //         piongrandparent == antigrandparent) {
            //           IsPairMinus = kTRUE;
            //           break;
            //         }
            //         dummyanti = antigrandparent->GetMother();
            //       }
            //       if (IsPairMinus) {
            //         break;
            //       }
            //     }
            //     dummypion = piongrandparent->GetMother();
            //   }
            //
            //   if (IsPairMinus) {
            //     PdgPairMinus = pionparent->GetPdgCode();
            //   }
            // }
          }
        }
      }

      // Top selections
      if (fAdditionHist && track->Charge() == 1 &&
          TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= fNTPCclusters &&
          nsigmaPionTPC < fTPCsigmas && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
          ph2->IsDispOK() && ph2->IsDisp2OK()) {
        FillHistogram("DCA_Daughters_Plus", dca, pt);
        FillHistogram("CPA_Daughters_Plus", cpa, pt);
        FillHistogram("RAD_Daughters_Plus", rad, pt);

        // DCAdaug variation
        if (cpa >= fCPAplusCut && TestRAD(pt, rad, fRADplusCut) &&
            (m >= msigmap - 0.02 && m <= msigmap + 0.02) && fIsMC) {
          FillHistogram("Track_Plus_DCAdaug", pt);
          for (Int_t i1 = 1; i1 < 11; i1++) {
            if (TestDCADaug(pt, dca, DCAcut[i1 - 1])) {
              FillHistogram(Form("Track_Plus_DCAdaug_%d", i1), pt);
            }
          }
        }
        // RAD variation
        if (cpa >= fCPAplusCut && TestDCADaug(pt, dca, fDCAdaugplusCut) &&
            (m >= msigmap - 0.02 && m <= msigmap + 0.02) && fIsMC) {
          FillHistogram("Track_Plus_RAD", pt);
          for (Int_t i1 = 1; i1 < 11; i1++) {
            if (TestRAD(pt, rad, RADcut[i1 - 1])) {
              FillHistogram(Form("Track_Plus_RAD_%d", i1), pt);
            }
          }
        }
        // CPA variation
        if (TestDCADaug(pt, dca, fDCAdaugplusCut) &&
            TestRAD(pt, rad, fRADplusCut) &&
            (m >= msigmap - 0.02 && m <= msigmap + 0.02) && fIsMC) {
          FillHistogram("Track_Plus_CPA", pt);
          for (Int_t i1 = 1; i1 < 11; i1++) {
            if (cpa >= CPAcut[i1 - 1]) {
              FillHistogram(Form("Track_Plus_CPA_%d", i1), pt);
            }
          }
        }
        // Test of cuts
        if (TestDCADaug(pt, dca, fDCAdaugplusCut)) {
          FillHistogram("DCA_Daughters_Test", dca, pt);
        }
        if (cpa >= fCPAplusCut) {
          FillHistogram("CPA_Daughters_Test", cpa, pt);
        }
        if (TestRAD(pt, rad, fRADplusCut)) {
          FillHistogram("RAD_Daughters_Test", rad, pt);
        }
      }
      if (fAdditionHist && track->Charge() == -1 &&
          TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= fNTPCclusters &&
          nsigmaPionTPC < fTPCsigmas && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
          ph2->IsDispOK() && ph2->IsDisp2OK()) {
        FillHistogram("DCA_Daughters_Minus", dca, pt);
        FillHistogram("CPA_Daughters_Minus", cpa, pt);
        FillHistogram("RAD_Daughters_Minus", rad, pt);

        // DCAdaug variation
        if (cpa >= fCPAminusCut && TestRAD(pt, rad, fRADminusCut) &&
            (m >= msigmam - 0.02 && m <= msigmam + 0.02) && fIsMC) {
          FillHistogram("Track_Minus_DCAdaug", pt);
          for (Int_t i1 = 1; i1 < 11; i1++) {
            if (TestDCADaug(pt, dca, DCAcut[i1 - 1])) {
              FillHistogram(Form("Track_Minus_DCAdaug_%d", i1), pt);
            }
          }
        }
        // RAD variation
        if (cpa >= fCPAminusCut && TestDCADaug(pt, dca, fDCAdaugminusCut) &&
            (m >= msigmam - 0.02 && m <= msigmam + 0.02) && fIsMC) {
          FillHistogram("Track_Minus_RAD", pt);
          for (Int_t i1 = 1; i1 < 11; i1++) {
            if (TestRAD(pt, rad, RADcut[i1 - 1])) {
              FillHistogram(Form("Track_Minus_RAD_%d", i1), pt);
            }
          }
        }
        // CPA variation
        if (TestDCADaug(pt, dca, fDCAdaugminusCut) &&
            TestRAD(pt, rad, fRADminusCut) &&
            (m >= msigmam - 0.02 && m <= msigmam + 0.02) && fIsMC) {
          FillHistogram("Track_Minus_CPA", pt);
          for (Int_t i1 = 1; i1 < 11; i1++) {
            if (cpa >= CPAcut[i1 - 1]) {
              FillHistogram(Form("Track_Minus_CPA_%d", i1), pt);
            }
          }
        }
      }

      // default
      if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= fNTPCclusters &&
          nsigmaPionTPC < fTPCsigmas && TestDCADaug(pt, dca, fDCAdaugplusCut) &&
          cpa >= fCPAplusCut && TestRAD(pt, rad, fRADplusCut) &&
          ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
          ph2->IsDisp2OK() && track->Charge() == 1) {
        FillHistogram("Cuts1_InvMass_Charge1", m, pt);
      }
      if (fQAhist) {
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
            TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
            ph2->IsDispOK() && track->Charge() == 1) {
          FillHistogram("Wo_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
            TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && track->Charge() == 1) {
          FillHistogram("Ncell_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
            TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
            ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Disp_InvMass_Charge1", m, pt);
        }
      }
      if (fInvMassHist && !fIsMC) {
        if (TMath::Abs(track->Eta()) <= 0.7 && TPCClust >= fNTPCclusters &&
            nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
            TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
            track->Charge() == 1) {
          FillHistogram("Cuts2_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= 0.9 && TPCClust >= fNTPCclusters &&
            nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
            TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
            track->Charge() == 1) {
          FillHistogram("Cuts3_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= 50 &&
            nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
            TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
            track->Charge() == 1) {
          FillHistogram("Cuts4_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= 70 &&
            nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
            TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
            track->Charge() == 1) {
          FillHistogram("Cuts5_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < 2.5 &&
            TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
            TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
            track->Charge() == 1) {
          FillHistogram("Cuts6_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < 3.5 &&
            TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
            TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
            track->Charge() == 1) {
          FillHistogram("Cuts7_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0. && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts8_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0. && TestRAD(pt, rad, 0.25) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts9_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0. && TestRAD(pt, rad, 0.3) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts10_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts11_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 &&
            TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
            ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts12_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 && TestRAD(pt, rad, 0.3) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts13_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts14_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 &&
            TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
            ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts15_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 && TestRAD(pt, rad, 0.3) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts16_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0. && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts17_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0. && TestRAD(pt, rad, 0.25) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts18_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0. && TestRAD(pt, rad, 0.3) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts19_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0.3 && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts20_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0.3 && TestRAD(pt, rad, 0.3) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts21_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts22_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 &&
            TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
            ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts23_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 && TestRAD(pt, rad, 0.3) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts24_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0. && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts25_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0. && TestRAD(pt, rad, 0.25) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts26_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0. && TestRAD(pt, rad, 0.3) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts27_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts28_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 &&
            TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
            ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts29_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 && TestRAD(pt, rad, 0.3) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts30_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts31_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 &&
            TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
            ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts32_InvMass_Charge1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 && TestRAD(pt, rad, 0.3) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == 1) {
          FillHistogram("Cuts33_InvMass_Charge1", m, pt);
        }
        if (ph2->GetTime() < 100.e-09) {
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
              TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
              track->Charge() == 1) {
            FillHistogram("Cuts34_InvMass_Charge1", m, pt);
          }
        }
      }
      // default
      if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= fNTPCclusters &&
          nsigmaPionTPC < fTPCsigmas &&
          TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
          TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
          ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
        FillHistogram("Cuts1_InvMass_Charge-1", m, pt);
      }
      if (fQAhist) {
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
            TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
            ph2->IsDispOK() && track->Charge() == -1) {
          FillHistogram("Wo_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
            TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && track->Charge() == -1) {
          FillHistogram("Ncell_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
            TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
            ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Disp_InvMass_Charge-1", m, pt);
        }
      }
      if (fInvMassHist && !fIsMC) {
        if (TMath::Abs(track->Eta()) <= 0.7 && TPCClust >= fNTPCclusters &&
            nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
            TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
            track->Charge() == -1) {
          FillHistogram("Cuts2_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= 0.9 && TPCClust >= fNTPCclusters &&
            nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
            TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
            track->Charge() == -1) {
          FillHistogram("Cuts3_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= 50 &&
            nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
            TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
            track->Charge() == -1) {
          FillHistogram("Cuts4_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= 70 &&
            nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
            TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
            track->Charge() == -1) {
          FillHistogram("Cuts5_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < 2.5 &&
            TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
            TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
            track->Charge() == -1) {
          FillHistogram("Cuts6_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < 3.5 &&
            TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
            TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
            ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
            track->Charge() == -1) {
          FillHistogram("Cuts7_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0. && TestRAD(pt, rad, 0.1) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts8_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0. && TestRAD(pt, rad, 0.15) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts9_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0. && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts10_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 && TestRAD(pt, rad, 0.1) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts11_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 &&
            TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
            ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts12_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts13_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 && TestRAD(pt, rad, 0.1) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts14_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 &&
            TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
            ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts15_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts16_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0. && TestRAD(pt, rad, 0.1) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts17_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0. && TestRAD(pt, rad, 0.15) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts18_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0. && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts19_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0.3 && TestRAD(pt, rad, 0.1) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts20_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0.3 && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts21_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 && TestRAD(pt, rad, 0.1) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts22_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 &&
            TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
            ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts23_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts24_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0. && TestRAD(pt, rad, 0.1) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts25_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0. && TestRAD(pt, rad, 0.15) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts26_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0. && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts27_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 && TestRAD(pt, rad, 0.1) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts28_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 &&
            TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
            ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts29_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts30_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 && TestRAD(pt, rad, 0.1) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts31_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 &&
            TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
            ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts32_InvMass_Charge-1", m, pt);
        }
        if (TMath::Abs(track->Eta()) <= fTrackEta &&
            TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
            TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 && TestRAD(pt, rad, 0.2) &&
            ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
            ph2->IsDisp2OK() && track->Charge() == -1) {
          FillHistogram("Cuts33_InvMass_Charge-1", m, pt);
        }
        if (ph2->GetTime() < 100.e-09) {
          if (TMath::Abs(track->Eta()) <= fTrackEta &&
              TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
              TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
              TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
              ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
              track->Charge() == -1) {
            FillHistogram("Cuts34_InvMass_Charge-1", m, pt);
          }
        }
      }

      if (fIsMC) {
        if (primLabelTrack > -1 && primLabel2 > -1) {
          // Track
          AliAODMCParticle *prim2 =
              (AliAODMCParticle *)fStack->At(primLabelTrack);
          // Int_t primTrackPdg = prim2->GetPdgCode();
          // TLorentzVector lv1;
          // prim2->Momentum(lv1);

          // Cluster
          AliAODMCParticle *prim3 = (AliAODMCParticle *)fStack->At(primLabel2);
          // TLorentzVector lv2;
          // prim3->Momentum(lv2);

          // Double_t mPrim = (lv1 + lv2).M();
          // Double_t ptPrim = (lv1 + lv2).Pt();

          // Int_t primLabelCluster1 = prim3->GetPdgCode();

          // if (IsPairPlus) {
          //   if (TestDCADaug(pt,dca,fDCAdaugplusCut) && cpa>=fCPAplusCut &&
          //   TestRAD(pt,rad,fRADplusCut) && ph2->IsCPVOK() && ph2->IsCPV2OK()
          //   && ph2->IsDispOK() && ph2->IsDisp2OK() && pidPion &&
          //   track->Charge()==1) {
          //     FillHistogram("Cuts1_InvMass_AntiSigmaPlus_Parent", m, pt);
          //     FillHistogram("Plus_Parent_PDG", PdgPairPlus);
          //   }
          // }
          // if (IsPairMinus) {
          //   if (TestDCADaug(pt,dca,fDCAdaugminusCut) && cpa>=fCPAminusCut &&
          //   TestRAD(pt,rad,fRADminusCut) && ph2->IsCPVOK() && ph2->IsCPV2OK()
          //   && ph2->IsDispOK() && ph2->IsDisp2OK() && pidPion &&
          //   track->Charge()==-1) {
          //     FillHistogram("Cuts1_InvMass_AntiSigmaMinus_Parent", m, pt);
          //     FillHistogram("Minus_Parent_PDG", PdgPairMinus);
          //   }
          // }

          if (IsSigmaPlusv2) {
            if (TMath::Abs(track->Eta()) <= fTrackEta &&
                TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
                TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
                ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
                track->Charge() == 1) {
              FillHistogram("Cuts1_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
            }
            if (fQAhist) {
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
                  TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
                  ph2->IsDispOK() && track->Charge() == 1) {
                FillHistogram("Wo_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
                  TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
                  ph2->IsCPV2OK() && ph2->IsDispOK() && track->Charge() == 1) {
                FillHistogram("Ncell_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
                  TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Disp_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
              }
            }
            if (fInvMassHist) {
              if (TMath::Abs(track->Eta()) <= 0.7 &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
                  TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
                  ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == 1) {
                FillHistogram("Cuts2_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
              }
              if (TMath::Abs(track->Eta()) <= 0.9 &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
                  TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
                  ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == 1) {
                FillHistogram("Cuts3_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= 50 &&
                  nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
                  TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
                  ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == 1) {
                FillHistogram("Cuts4_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= 70 &&
                  nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
                  TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
                  ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == 1) {
                FillHistogram("Cuts5_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < 2.5 &&
                  TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
                  TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
                  ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == 1) {
                FillHistogram("Cuts6_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < 3.5 &&
                  TestDCADaug(pt, dca, fDCAdaugplusCut) && cpa >= fCPAplusCut &&
                  TestRAD(pt, rad, fRADplusCut) && ph2->IsCPVOK() &&
                  ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == 1) {
                FillHistogram("Cuts7_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts8_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts9_InvMass_AntiSigmaPlus_ParentCheck", m, pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts10_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts11_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts12_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts13_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts14_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts15_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts16_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts17_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts18_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts19_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts20_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts21_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts22_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts23_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts24_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts25_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts26_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts27_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts28_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts29_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts30_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts31_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.25) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts32_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.3) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() && track->Charge() == 1) {
                FillHistogram("Cuts33_InvMass_AntiSigmaPlus_ParentCheck", m,
                              pt);
              }
              if (ph2->GetTime() < 100.e-09) {
                if (TMath::Abs(track->Eta()) <= fTrackEta &&
                    TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                    TestDCADaug(pt, dca, fDCAdaugplusCut) &&
                    cpa >= fCPAplusCut && TestRAD(pt, rad, fRADplusCut) &&
                    ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
                    ph2->IsDisp2OK() && track->Charge() == 1) {
                  FillHistogram("Cuts34_InvMass_AntiSigmaPlus_ParentCheck", m,
                                pt);
                }
              }
            }
          }

          if (IsSigmaMinusv2) {
            if (TMath::Abs(track->Eta()) <= fTrackEta &&
                TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                TestDCADaug(pt, dca, fDCAdaugminusCut) && cpa >= fCPAminusCut &&
                TestRAD(pt, rad, fRADminusCut) && ph2->IsCPVOK() &&
                ph2->IsCPV2OK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
                track->Charge() == -1) {
              FillHistogram("Cuts1_InvMass_AntiSigmaMinus_ParentCheck", m, pt);
            }
            if (fQAhist) {
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, fDCAdaugminusCut) &&
                  cpa >= fCPAminusCut && TestRAD(pt, rad, fRADminusCut) &&
                  ph2->IsCPVOK() && ph2->IsDispOK() && track->Charge() == -1) {
                FillHistogram("Wo_InvMass_AntiSigmaMinus_ParentCheck", m, pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, fDCAdaugminusCut) &&
                  cpa >= fCPAminusCut && TestRAD(pt, rad, fRADminusCut) &&
                  ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
                  track->Charge() == -1) {
                FillHistogram("Ncell_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, fDCAdaugminusCut) &&
                  cpa >= fCPAminusCut && TestRAD(pt, rad, fRADminusCut) &&
                  ph2->IsCPVOK() && ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Disp_InvMass_AntiSigmaMinus_ParentCheck", m, pt);
              }
            }
            if (fInvMassHist) {
              if (TMath::Abs(track->Eta()) <= 0.7 &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, fDCAdaugminusCut) &&
                  cpa >= fCPAminusCut && TestRAD(pt, rad, fRADminusCut) &&
                  ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
                  ph2->IsDisp2OK() && track->Charge() == -1) {
                FillHistogram("Cuts2_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= 0.9 &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, fDCAdaugminusCut) &&
                  cpa >= fCPAminusCut && TestRAD(pt, rad, fRADminusCut) &&
                  ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
                  ph2->IsDisp2OK() && track->Charge() == -1) {
                FillHistogram("Cuts3_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= 50 &&
                  nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, fDCAdaugminusCut) &&
                  cpa >= fCPAminusCut && TestRAD(pt, rad, fRADminusCut) &&
                  ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
                  ph2->IsDisp2OK() && track->Charge() == -1) {
                FillHistogram("Cuts4_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta && TPCClust >= 70 &&
                  nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, fDCAdaugminusCut) &&
                  cpa >= fCPAminusCut && TestRAD(pt, rad, fRADminusCut) &&
                  ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
                  ph2->IsDisp2OK() && track->Charge() == -1) {
                FillHistogram("Cuts5_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < 2.5 &&
                  TestDCADaug(pt, dca, fDCAdaugminusCut) &&
                  cpa >= fCPAminusCut && TestRAD(pt, rad, fRADminusCut) &&
                  ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
                  ph2->IsDisp2OK() && track->Charge() == -1) {
                FillHistogram("Cuts6_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < 3.5 &&
                  TestDCADaug(pt, dca, fDCAdaugminusCut) &&
                  cpa >= fCPAminusCut && TestRAD(pt, rad, fRADminusCut) &&
                  ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
                  ph2->IsDisp2OK() && track->Charge() == -1) {
                FillHistogram("Cuts7_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts8_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts9_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts10_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts11_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts12_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts13_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts14_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts15_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.05) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts16_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts17_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts18_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts19_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts20_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts21_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts22_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts23_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.06) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts24_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts25_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts26_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0. &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts27_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts28_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts29_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0.3 &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts30_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.1) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts31_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.15) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts32_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (TMath::Abs(track->Eta()) <= fTrackEta &&
                  TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                  TestDCADaug(pt, dca, 0.07) && cpa >= 0.5 &&
                  TestRAD(pt, rad, 0.2) && ph2->IsCPVOK() && ph2->IsCPV2OK() &&
                  ph2->IsDispOK() && ph2->IsDisp2OK() &&
                  track->Charge() == -1) {
                FillHistogram("Cuts33_InvMass_AntiSigmaMinus_ParentCheck", m,
                              pt);
              }
              if (ph2->GetTime() < 100.e-09) {
                if (TMath::Abs(track->Eta()) <= fTrackEta &&
                    TPCClust >= fNTPCclusters && nsigmaPionTPC < fTPCsigmas &&
                    TestDCADaug(pt, dca, fDCAdaugminusCut) &&
                    cpa >= fCPAminusCut && TestRAD(pt, rad, fRADminusCut) &&
                    ph2->IsCPVOK() && ph2->IsCPV2OK() && ph2->IsDispOK() &&
                    ph2->IsDisp2OK() && track->Charge() == -1) {
                  FillHistogram("Cuts34_InvMass_AntiSigmaMinus_ParentCheck", m,
                                pt);
                }
              }
            }
          }
        }
      }
    }
  }
}

void AliAnalysisSigmaBarCharged::FillHistogram(const char *key,
                                               Double_t x) const {
  TH1 *hist = dynamic_cast<TH1 *>(fOutputContainer->FindObject(key));
  if (hist) {
    hist->Fill(x);
  } else {
    AliError(Form("can not find histogram (of instance TH1) <%s> ", key));
  }
}

void AliAnalysisSigmaBarCharged::FillHistogram(const char *key, Double_t x,
                                               Double_t y) const {
  TObject *obj = fOutputContainer->FindObject(key);
  TH1 *th1 = dynamic_cast<TH1 *>(obj);
  if (th1) {
    th1->Fill(x, y);
  } else {
    TH2 *th2 = dynamic_cast<TH2 *>(obj);
    if (th2) {
      th2->Fill(x, y);
    } else {
      AliError(Form("can not find histogram (of instance TH1) <%s> ", key));
    }
  }
}

void AliAnalysisSigmaBarCharged::FillHistogram(const char *key, Double_t x,
                                               Double_t y, Double_t z) const {
  TObject *obj = fOutputContainer->FindObject(key);

  TH2 *th2 = dynamic_cast<TH2 *>(obj);
  if (th2) {
    th2->Fill(x, y, z);
    return;
  }

  TH3 *th3 = dynamic_cast<TH3 *>(obj);
  if (th3) {
    th3->Fill(x, y, z);
    return;
  }

  AliError(Form("can not findi histogram (of instance TH2) <%s> ", key));
}
void AliAnalysisSigmaBarCharged::Terminate(Option_t *) {}

Double_t CBfunction(Double_t *x, Double_t *par) {
  // Parameterizatin of signal
  Double_t m = par[1];
  Double_t s = par[2];
  Double_t dx = (x[0] - m) / s;
  Double_t alpha = par[3];
  Double_t n = par[4];

  Double_t result;
  Double_t fact1TLessMinosAlphaL = alpha / n;
  Double_t fact2TLessMinosAlphaL = (n / alpha) - alpha - dx;
  Double_t fact1THihgerAlphaH = alpha / n;
  Double_t fact2THigherAlphaH = (n / alpha) - alpha + dx;

  if (-alpha <= dx && alpha >= dx) {
    result = exp(-0.5 * dx * dx);
  } else if (dx < -alpha) {
    result = exp(-0.5 * alpha * alpha) *
             pow(fact1TLessMinosAlphaL * fact2TLessMinosAlphaL, -n);
  } else if (dx > alpha) {
    result = exp(-0.5 * alpha * alpha) *
             pow(fact1THihgerAlphaH * fact2THigherAlphaH, -n);
  }
  return par[0] * result;
}

Double_t AliAnalysisSigmaBarCharged::RealRes(Double_t x) {
  // Simulates time resolution
  if (x <= 0.2) {
    x = 0.2;
  }
  if (x >= 10) {
    x = 10;
  }
  return sqrt(pow(9.217682 * TMath::Exp(-x / 3.575588e+11) / x -
                      0.02196338 * x + 0.02202962 * x * x,
                  2) -
              pow(0.5, 2)) *
         1e-9;
}

Double_t AliAnalysisSigmaBarCharged::RealResv2(Double_t x) {
  Double_t sigma, alpha, n;
  // Sigma
  if (x <= 0.4) {
    sigma = 9.27927e-09 * exp(0.4 * 1.24013e-01) / 0.4 - 1.36967e-09;
  } else if (x >= 8.) {
    sigma = 9.27927e-09 * exp(8. * 1.24013e-01) / 8. - 1.36967e-09;
  } else {
    sigma = 9.27927e-09 * exp(x * 1.24013e-01) / x - 1.36967e-09;
  }

  sigma = sqrt(sigma * sigma - 0.5 * 0.5 * 1e-9 * 1e-9);

  // Alpha
  if (x <= 0.4) {
    alpha = exp(-1.20237e+00 * 0.4 + 6.78364e-01) - 1.34970e-02 +
            5.21402e-01 * 0.4 + 1.01185e-01 * 0.4 * 0.4 -
            2.90393e-02 * 0.4 * 0.4 * 0.4 + 1.48951e-03 * 0.4 * 0.4 * 0.4 * 0.4;
  } else if (x >= 10.) {
    alpha = exp(-1.20237e+00 * 10. + 6.78364e-01) - 1.34970e-02 +
            5.21402e-01 * 10. + 1.01185e-01 * 10. * 10. -
            2.90393e-02 * 10. * 10. * 10. + 1.48951e-03 * 10. * 10. * 10. * 10.;
  } else {
    alpha = exp(-1.20237e+00 * x + 6.78364e-01) - 1.34970e-02 +
            5.21402e-01 * x + 1.01185e-01 * x * x - 2.90393e-02 * x * x * x +
            1.48951e-03 * x * x * x * x;
  }

  // N
  if (x <= 0.5) {
    n = exp(2.51256e+00 - 7.40150e-01 * 0.5) + 2.19880e-01;
  } else {
    n = exp(2.51256e+00 - 7.40150e-01 * x) + 2.19880e-01;
  }

  Int_t bin, binx, ibin, loop;
  Double_t r1, xres;

  TF1 *f1 = new TF1("CB", CBfunction, -250e-9, 250e-9, 5);
  f1->SetParameters(1, 0, sigma, alpha, n);

  TAxis *xAxis = new TAxis(2000, -250e-9, 250e-9);

  Int_t first = xAxis->GetFirst();
  Int_t last = xAxis->GetLast();
  Int_t nbinsx = last - first + 1;

  Double_t *integral = new Double_t[nbinsx + 1];
  integral[0] = 0;
  for (binx = 1; binx <= nbinsx; binx++) {
    Double_t fint = f1->Integral(xAxis->GetBinLowEdge(binx + first - 1),
                                 xAxis->GetBinUpEdge(binx + first - 1), 0.);
    integral[binx] = integral[binx - 1] + fint;
  }

  //   - Normalize integral to 1
  for (bin = 1; bin <= nbinsx; bin++)
    integral[bin] /= integral[nbinsx];

  r1 = gRandom->Rndm();
  ibin = TMath::BinarySearch(nbinsx, &integral[0], r1);
  xres = xAxis->GetBinLowEdge(ibin + first) +
         xAxis->GetBinWidth(ibin + first) * (r1 - integral[ibin]) /
             (integral[ibin + 1] - integral[ibin]);

  delete[] integral;
  delete xAxis;
  delete f1;

  return xres;
}

Bool_t AliAnalysisSigmaBarCharged::TestDCADaug(Double_t pt, Double_t dca,
                                               Double_t dcashift) {
  // Daughter DCA cut
  if (dca < dcashift + exp(-1.381 * pt - 2.232)) {
    return true;
  } else {
    return false;
  }
}

Bool_t AliAnalysisSigmaBarCharged::TestRAD(Double_t pt, Double_t rad,
                                           Double_t radshift) {
  // RAD Plus cut
  if (rad > 0.193 * pt + radshift) {
    return true;
  } else {
    return false;
  }
}

void AliAnalysisSigmaBarCharged::PropagateToDCACurvedBachelor(
    Double_t *decayparams, AliCaloPhoton *v, Double_t v0vtx[3],
    AliExternalTrackParam *t, Double_t b) {
  //--------------------------------------------------------------------
  // This function returns the DCA between the V0 and the track
  // assumes that bachelor track is not straight
  // algorithm based on AliExternalTrackParam::GetDCA with zero curvature track
  //--------------------------------------------------------------------

  // Double_t decayparams[4]; //decayparams[0]-decayparams[2] -- cascade vertex
  // (secondary vtx); decayparams[3] -- dca daughters

  // Double_t alpha=t->GetAlpha(), cs1=TMath::Cos(alpha), sn1=TMath::Sin(alpha);
  Double_t r[3];
  t->GetXYZ(r);
  // Double_t x1=r[0], y1=r[1], z1=r[2];
  Double_t p[3];
  t->GetPxPyPz(p);
  // Double_t px1=p[0], py1=p[1], pz1=p[2];

  Double_t x2 = v0vtx[0], y2 = v0vtx[1],
           z2 = v0vtx[2]; // position and momentum of V0
  Double_t px2 = v->Px(), py2 = v->Py(), pz2 = v->Pz();

  // v->GetXYZ(x2,y2,z2);
  // v->GetPxPyPz(px2,py2,pz2);

  decayparams[3] = 1e+33; // initial dca
  Double_t dy2 = 1e-10;
  Double_t dz2 = 1e-10;
  Double_t dx2 = 1e-10;

  // Create dummy V0 track
  // V0 properties to get started
  Double_t xyz[3] = {v0vtx[0], v0vtx[1], v0vtx[2]},
           pxpypz[3] = {v->Px(), v->Py(), v->Pz()}, cv[21];
  for (Int_t ii = 0; ii < 21; ii++)
    cv[ii] = 0.0; // something small

  // Mockup track for V0 trajectory (no covariance)
  // AliExternalTrackParam *hV0Traj = new
  // AliExternalTrackParam(xyz,pxpypz,cv,+1);
  AliExternalTrackParam lV0TrajObject(xyz, pxpypz, cv, +1),
      *hV0Traj = &lV0TrajObject;
  hV0Traj->ResetCovariance(1); // won't use

  // Re-acquire helix parameters for bachelor (necessary!)
  Double_t p1[8];
  t->GetHelixParameters(p1, b);
  p1[6] = TMath::Sin(p1[2]);
  p1[7] = TMath::Cos(p1[2]);

  Double_t p2[8];
  hV0Traj->GetHelixParameters(
      p2, 0.0); // p2[4]=0 -> no curvature (fine, predicted in Evaluate)
  p2[6] = TMath::Sin(p2[2]);
  p2[7] = TMath::Cos(p2[2]);

  Double_t r1[3], g1[3], gg1[3];
  Double_t t1 = 0.;
  Evaluate(p1, t1, r1, g1, gg1);
  Double_t r2[3], g2[3], gg2[3];
  Double_t t2 = 0.;
  Evaluate(p2, t2, r2, g2, gg2);

  Double_t dx = r2[0] - r1[0], dy = r2[1] - r1[1], dz = r2[2] - r1[2];
  Double_t dm = dx * dx / dx2 + dy * dy / dy2 + dz * dz / dz2;

  Int_t max = 27; // standard in AliExternalTrackParam::GetDCA, good performance
  while (max--) {
    Double_t gt1 = -(dx * g1[0] / dx2 + dy * g1[1] / dy2 + dz * g1[2] / dz2);
    Double_t gt2 = +(dx * g2[0] / dx2 + dy * g2[1] / dy2 + dz * g2[2] / dz2);
    Double_t h11 = (g1[0] * g1[0] - dx * gg1[0]) / dx2 +
                   (g1[1] * g1[1] - dy * gg1[1]) / dy2 +
                   (g1[2] * g1[2] - dz * gg1[2]) / dz2;
    Double_t h22 = (g2[0] * g2[0] + dx * gg2[0]) / dx2 +
                   (g2[1] * g2[1] + dy * gg2[1]) / dy2 +
                   (g2[2] * g2[2] + dz * gg2[2]) / dz2;
    Double_t h12 =
        -(g1[0] * g2[0] / dx2 + g1[1] * g2[1] / dy2 + g1[2] * g2[2] / dz2);

    Double_t det = h11 * h22 - h12 * h12;

    Double_t dt1, dt2;
    if (TMath::Abs(det) < 1.e-33) {
      //(quasi)singular Hessian
      dt1 = -gt1;
      dt2 = -gt2;
    } else {
      dt1 = -(gt1 * h22 - gt2 * h12) / det;
      dt2 = -(h11 * gt2 - h12 * gt1) / det;
    }

    if ((dt1 * gt1 + dt2 * gt2) > 0) {
      dt1 = -dt1;
      dt2 = -dt2;
    }

    // check delta(phase1) ?
    // check delta(phase2) ?

    if (TMath::Abs(dt1) / (TMath::Abs(t1) + 1.e-3) < 1.e-4)
      if (TMath::Abs(dt2) / (TMath::Abs(t2) + 1.e-3) < 1.e-4) {
        if ((gt1 * gt1 + gt2 * gt2) > 1.e-4 / dy2 / dy2) {
          AliDebug(1, " stopped at not a stationary point !");
        }
        Double_t lmb = h11 + h22;
        lmb = lmb - TMath::Sqrt(lmb * lmb - 4 * det);
        if (lmb < 0.) {
          AliDebug(1, " stopped at not a minimum !");
        }
        break;
      }

    Double_t dd = dm;
    for (Int_t div = 1;; div *= 2) {
      Evaluate(p1, t1 + dt1, r1, g1, gg1);
      Evaluate(p2, t2 + dt2, r2, g2, gg2);
      dx = r2[0] - r1[0];
      dy = r2[1] - r1[1];
      dz = r2[2] - r1[2];
      dd = dx * dx / dx2 + dy * dy / dy2 + dz * dz / dz2;
      if (dd < dm)
        break;
      dt1 *= 0.5;
      dt2 *= 0.5;
      if (div > 512) {
        AliDebug(1, " overshoot !");
        break;
      }
    }
    dm = dd;
    t1 += dt1;
    t2 += dt2;
  }
  if (max <= 0) {
    AliDebug(1, " too many iterations !");
  }
  Double_t cs = TMath::Cos(t->GetAlpha());
  Double_t sn = TMath::Sin(t->GetAlpha());
  Double_t xthis = r1[0] * cs + r1[1] * sn;

  // Propagate bachelor to the point of DCA
  if (!t->PropagateTo(xthis, b)) {
    AliDebug(1, " propagation failed !");
    decayparams[3] = 1e+33;
  }

  // V0 distance to bachelor: the desired distance
  Double_t rBachDCAPt[3];
  t->GetXYZ(rBachDCAPt);

  // dca between daughters
  decayparams[3] = GetDCASigmabar(xyz, pxpypz, rBachDCAPt); // dca is ok

  // reconstruction of secondary vertex
  // Double_t r[3]; t.GetXYZ(r);
  Double_t x1 = rBachDCAPt[0], y1 = rBachDCAPt[1],
           z1 = rBachDCAPt[2]; // position of the bachelor
  t->GetPxPyPz(p);
  Double_t px1 = p[0], py1 = p[1], pz1 = p[2]; // momentum of the bachelor track

  // Double_t x2,y2,z2;          // position of the V0
  // v.GetXYZ(x2,y2,z2);
  // Double_t px2,py2,pz2;       // momentum of V0
  // v.GetPxPyPz(px2,py2,pz2);

  Double_t a2 = ((x1 - x2) * px2 + (y1 - y2) * py2 + (z1 - z2) * pz2) /
                (px2 * px2 + py2 * py2 + pz2 * pz2);

  Double_t xm = x2 + a2 * px2;
  Double_t ym = y2 + a2 * py2;
  Double_t zm = z2 + a2 * pz2;

  // position of the cascade decay
  decayparams[0] = 0.5 * (x1 + xm);
  decayparams[1] = 0.5 * (y1 + ym);
  decayparams[2] = 0.5 * (z1 + zm);
}

Float_t AliAnalysisSigmaBarCharged::GetDCASigmabar(Double_t xyz[3],
                                                   Double_t pxpypz[3],
                                                   Double_t bachelor[3]) const {
  //--------------------------------------------------------------------
  // This function returns V0's impact parameter calculated in 3D
  //--------------------------------------------------------------------
  Double_t x = xyz[0], y = xyz[1], z = xyz[2];
  Double_t px = pxpypz[0];
  Double_t py = pxpypz[1];
  Double_t pz = pxpypz[2];

  Double_t dx = (bachelor[1] - y) * pz - (bachelor[2] - z) * py;
  Double_t dy = (bachelor[0] - x) * pz - (bachelor[2] - z) * px;
  Double_t dz = (bachelor[0] - x) * py - (bachelor[1] - y) * px;
  Double_t d = TMath::Sqrt((dx * dx + dy * dy + dz * dz) /
                           (px * px + py * py + pz * pz));
  return d;
}

void AliAnalysisSigmaBarCharged::Evaluate(const Double_t *h, Double_t t,
                                          Double_t r[3],  // radius vector
                                          Double_t g[3],  // first defivatives
                                          Double_t gg[3]) // second derivatives
{
  //--------------------------------------------------------------------
  // Calculate position of a point on a track and some derivatives
  //--------------------------------------------------------------------
  Double_t phase = h[4] * t + h[2];
  Double_t sn = TMath::Sin(phase), cs = TMath::Cos(phase);

  r[0] = h[5];
  r[1] = h[0];
  if (TMath::Abs(h[4]) > kAlmost0) {
    r[0] += (sn - h[6]) / h[4];
    r[1] -= (cs - h[7]) / h[4];
  } else {
    r[0] += t * cs;
    r[1] -= -t * sn;
  }
  r[2] = h[1] + h[3] * t;

  g[0] = cs;
  g[1] = sn;
  g[2] = h[3];

  gg[0] = -h[4] * sn;
  gg[1] = h[4] * cs;
  gg[2] = 0.;
}

Double_t AliAnalysisSigmaBarCharged::GetCosineOfPointingAngle(
    Double_t mom[3], Double_t pos[3], Double_t refPointX, Double_t refPointY,
    Double_t refPointZ) const {
  // calculates the pointing angle of the cascade wrt a reference point

  Double_t momCas[3] = {mom[0], mom[1], mom[2]}; // momentum of the cascade

  Double_t fPosXi[3] = {pos[0], pos[1], pos[2]}; // pos of cascade

  Double_t
      deltaPos[3]; // vector between the reference point and the cascade vertex
  deltaPos[0] = fPosXi[0] - refPointX;
  deltaPos[1] = fPosXi[1] - refPointY;
  deltaPos[2] = fPosXi[2] - refPointZ;

  Double_t momCas2 =
      momCas[0] * momCas[0] + momCas[1] * momCas[1] + momCas[2] * momCas[2];
  Double_t deltaPos2 = deltaPos[0] * deltaPos[0] + deltaPos[1] * deltaPos[1] +
                       deltaPos[2] * deltaPos[2];

  Double_t cosinePointingAngle =
      (deltaPos[0] * momCas[0] + deltaPos[1] * momCas[1] +
       deltaPos[2] * momCas[2]) /
      TMath::Sqrt(momCas2 * deltaPos2);

  return cosinePointingAngle;
}
void AliAnalysisSigmaBarCharged::InitGeometry() {
  // Rotation matrixes are set with Tender

  if (fPHOSgeom)
    return;

  fPHOSgeom = AliPHOSGeometry::GetInstance();

  if (!fPHOSgeom) {            // Geometry not yet constructed with Tender
    if (fRunNumber < 209122) { // Run1
      AliError("TaggedPhotons: Can not get geometry from TENDER, creating PHOS "
               "geometry for Run1\n");
      fPHOSgeom = AliPHOSGeometry::GetInstance("IHEP", "");
    } else {
      AliError("TaggedPhotons: Can not get geometry from TENDER, creating PHOS "
               "geometry for Run2\n");
      fPHOSgeom = AliPHOSGeometry::GetInstance("Run2", "");
    }
    AliOADBContainer geomContainer("phosGeo");
    geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root",
                               "PHOSRotationMatrixes");
    TObjArray *matrixes = (TObjArray *)geomContainer.GetObject(
        fRunNumber, "PHOSRotationMatrixes");
    for (Int_t mod = 0; mod < 5; mod++) {
      if (!matrixes->At(mod))
        continue;
      fPHOSgeom->SetMisalMatrix(((TGeoHMatrix *)matrixes->At(mod)), mod);
    }
  }

  // Read BadMap for MC simulations
  AliOADBContainer badmapContainer(Form("phosBadMap"));
  badmapContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSBadMaps.root",
                               "phosBadMap");
  TObjArray *maps =
      (TObjArray *)badmapContainer.GetObject(fRunNumber, "phosBadMap");
  if (!maps) {
    AliError("TaggedPhotons: Can not read Bad map\n");
  } else {
    AliInfo(Form("TaggedPhotons: Setting PHOS bad map with name %s \n",
                 maps->GetName()));
    for (Int_t mod = 0; mod < 5; mod++) {
      if (fPHOSBadMap[mod])
        delete fPHOSBadMap[mod];
      TH2I *h = (TH2I *)maps->At(mod);
      if (h)
        fPHOSBadMap[mod] = new TH2I(*h);
    }
  }
}
