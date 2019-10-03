/**************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
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

#include <TClonesArray.h>
#include <THashList.h>

#include <AliVEventHandler.h>
#include <AliAnalysisManager.h>
#include <AliClusterContainer.h>
#include <AliEMCALTriggerPatchInfo.h>

#include "AliEmcalTriggerSimQATask.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerSimQATask);
/// \endcond

/**
 * Dummy constructor
 */
AliEmcalTriggerSimQATask::AliEmcalTriggerSimQATask() : AliAnalysisTaskEmcal(),
fTriggerPatchesName("EmcalTriggers"),
fTriggerPatches(0),
fMinAmplitude(0),
fMinClusterEnergy(3.),
fPtBinWidth(0.5),
fMaxPt(120),
fEventTriggerBits(0x0),
fHistManager("AliEmcalTriggerSimQATask")
{}

/**
 * Named Constructor
 * \param name Name of the Trigger Simulation QA Task
 */
AliEmcalTriggerSimQATask::AliEmcalTriggerSimQATask(const char *name) : AliAnalysisTaskEmcal(name,kTRUE),
fTriggerPatchesName("EmcalTriggers"),
fTriggerPatches(0),
fMinAmplitude(0),
fMinClusterEnergy(3.),
fPtBinWidth(0.5),
fMaxPt(120),
fEventTriggerBits(0x0),
fHistManager(name)
{
  SetMakeGeneralHistograms(kTRUE);

  // Initializing array in CINT-compatible way
  //EventEMCALTriggerType_t fTriggerTypesValues[kNTriggerTypes] = {kMB, kEL0, kEG1, kEG2, kEJ1, kEJ2};
  EventEMCALTriggerType_t fTriggerTypesValues[kNTriggerTypes] = {kMB, kEL0, kDL0, kEG1, kDG1, kEG2, kDG2, kEJ1, kDJ1, kEJ2, kDJ2};
  memcpy (fTriggerTypes,fTriggerTypesValues,sizeof(fTriggerTypes));
}


/**
 * Destructor
 */
AliEmcalTriggerSimQATask::~AliEmcalTriggerSimQATask() {

}


/**
 * Init the analysis
 */
void AliEmcalTriggerSimQATask::ExecOnce() {
  AliAnalysisTaskEmcal::ExecOnce();
}

/**
 * Create objects, histograms
 */
void AliEmcalTriggerSimQATask::UserCreateOutputObjects() {
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  TString histname;
  TString title;

  Int_t nPtBins = TMath::CeilNint(fMaxPt/fPtBinWidth);

  Double_t fMinEta = -0.7;
  Double_t fMaxEta = 0.7;
  Double_t fMinDeltaEta = -1.4;
  Double_t fMaxDeltaEta = 1.4;
  Double_t fMinPhi = 0;
  Double_t fMaxPhi = 2*TMath::Pi();
  Double_t fMinDeltaPhi = -TMath::Pi();
  Double_t fMaxDeltaPhi = TMath::Pi();
  Int_t nEtaBins = 96;
  Int_t nPhiBins = 208;
  Int_t nDeltaEtaBins = 192;
  Int_t nDeltaPhiBins = 416;

  Double_t fMinCol = 0;
  Double_t fMaxCol = 48;
  Double_t fMinRow = 0;
  Double_t fMaxRow = 104;
  Int_t nCol = 48;
  Int_t nRow = 104;

  Double_t fMaxADCValue = 1024;
  Int_t nADCBins = 512;

  Int_t nSmallADCBins = 128;
  Int_t nSmallEBins   = 100;

  // Hist for counting clusters
  fHistManager.CreateTH1("NClusters","NClusters;N_{cluster}; counts",300,0,300);

  // loop over trigger types
  for (Int_t i = 0; i < kNTriggerTypes; i++) {
    // Clusters in the EMCal (not including DCal)
    histname = TString::Format("fHistEMCalClusEnergy_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#it{E}_{cluster} (GeV); counts";
    fHistManager.CreateTH1(histname.Data(),title.Data(),nPtBins,0,fMaxPt);
    // Clusters in the DCal
    histname = TString::Format("fHistDCalClusEnergy_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#it{E}_{cluster} (GeV); counts";
    fHistManager.CreateTH1(histname.Data(),title.Data(),nPtBins,0,fMaxPt);
  }

  for (Int_t i = 0; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistClusEtaPhi_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#eta_{cluster};#phi_{cluster};counts";
    fHistManager.CreateTH2(histname.Data(),title.Data(),nEtaBins,fMinEta,fMaxEta,nPhiBins,fMinPhi,fMaxPhi);
  }

  // Patch Energy Spectra
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistPatchEnergy_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#it{E}_{patch} (GeV); counts";
    fHistManager.CreateTH1(histname.Data(),title.Data(),nPtBins,0,fMaxPt);
  }
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistPatchOnlineADCAmp_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#it{E}_{patch} (ADC); counts";
    fHistManager.CreateTH1(histname.Data(),title.Data(),nADCBins,0,fMaxADCValue);
  }

  // Patch Energy Spectra vs Col
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistPatchColEnergy_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";Col;#it{E}_{patch} (ADC); counts";
    fHistManager.CreateTH2(histname.Data(),title.Data(),nCol,fMinCol,fMaxCol,nPtBins,0,fMaxPt);
  }
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistPatchColOnlineADCAmp_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";Col;#it{E}_{patch} (ADC); counts";
    fHistManager.CreateTH2(histname.Data(),title.Data(),nCol,fMinCol,fMaxCol,nADCBins,0,fMaxADCValue);
  }

  // Patch Energy Spectra vs Row
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistPatchRowEnergy_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";Row;#it{E}_{patch} (ADC); counts";
    fHistManager.CreateTH2(histname.Data(),title.Data(),nRow,fMinRow,fMaxRow,nPtBins,0,fMaxPt);
  }
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistPatchRowOnlineADCAmp_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";Row;#it{E}_{patch} (ADC); counts";
    fHistManager.CreateTH2(histname.Data(),title.Data(),nRow,fMinRow,fMaxRow,nADCBins,0,fMaxADCValue);
  }

  // Patch Corner Col,Row
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistPatchColRow_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";Col;Row;counts";
    fHistManager.CreateTH2(histname.Data(),title.Data(),nCol,fMinCol,fMaxCol,nRow,fMinRow,fMaxRow);
  }

  // Patch Geometric Centers
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistPatchEtaPhiGeo_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#eta_{patch};#phi_{patch};counts";
    fHistManager.CreateTH2(histname.Data(),title.Data(),nEtaBins,fMinEta,fMaxEta,nPhiBins,fMinPhi,fMaxPhi);
  }

  // Patch Cluster Matching Histograms
  Int_t nTRU = 0;
  if (!fGeom) {
    AliError("Geometry not found. Initializing...");
    fGeom = AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
    if (!fGeom) {
      AliError("Geometry instance could not be found. Assuming 52 TRU.");
      nTRU = 52;
    } else {
      nTRU = fGeom->GetNTotalTRU();
    }
  } else {
    nTRU = fGeom->GetNTotalTRU();
  }

  // Patch-Cluster Delta Phi vs Cluster Phi
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistClusPhiClusPatchDPhi_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#phi_{cluster};#Delta#phi_{cluster-patch};counts";
    fHistManager.CreateTH2(histname.Data(),title.Data(),nPhiBins,fMinPhi,fMaxPhi,nDeltaPhiBins,fMinDeltaPhi,fMaxDeltaPhi);
  }
  // Patch-Cluster Delta Phi vs Cluster Eta
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistClusEtaClusPatchDPhi_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#eta_{cluster};#Delta#phi_{cluster-patch};counts";
    fHistManager.CreateTH2(histname.Data(),title.Data(),nEtaBins,fMinEta,fMaxEta,nDeltaPhiBins,fMinDeltaPhi,fMaxDeltaPhi);
  }
  // Patch-Cluster Delta Eta vs Cluster Phi
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistClusPhiClusPatchDEta_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#phi_{cluster};#Delta#eta_{cluster-patch};counts";
    fHistManager.CreateTH2(histname.Data(),title.Data(),nPhiBins,fMinPhi,fMaxPhi,nDeltaEtaBins,fMinDeltaEta,fMaxDeltaEta);
  }
  // Patch-Cluster Delta Eta vs Cluster Eta
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistClusEtaClusPatchDEta_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#eta_{cluster};#Delta#eta_{cluster-patch};counts";
    fHistManager.CreateTH2(histname.Data(),title.Data(),nEtaBins,fMinEta,fMaxEta,nDeltaEtaBins,fMinDeltaEta,fMaxDeltaEta);
  }

  // Patch-Cluster Delta Eta, per TRU
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistClusPatchDEtaTRU_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#Delta#eta_{cluster-patch};TRU;counts";
    fHistManager.CreateTH2(histname.Data(),title.Data(),nDeltaEtaBins,fMinDeltaEta,fMaxDeltaEta,nTRU,0,nTRU);
  }
  // Patch-Cluster Delta Phi, per TRU
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistClusPatchDPhiTRU_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#Delta#phi_{cluster-patch};TRU;counts";
    fHistManager.CreateTH2(histname.Data(),title.Data(),nDeltaPhiBins,fMinDeltaPhi,fMaxDeltaPhi,nTRU,0,nTRU);
  }

  // Patch-Cluster Energy, per TRU
  for (Int_t i = 1; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistClusEnergyPatchAmpTRU_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";E_{clus} (GeV);E_{patch} (ADC);TRU;counts";
    fHistManager.CreateTH3(histname.Data(),title.Data(),nSmallEBins,0,fMaxPt,nSmallADCBins,0,fMaxADCValue,nTRU,0,nTRU);
  }


  fOutput->Add(fHistManager.GetListOfHistograms());
}

/**
 * Run analysis
 * \return Always true.
 */
Bool_t AliEmcalTriggerSimQATask::Run() {
  if (!fCaloClusters)
  {
    fCaloClusters = (TClonesArray*)GetClusterContainer(0);
  }

  return kTRUE;
}

/**
 * Fill Histograms
 * \return Always true.
 */
Bool_t AliEmcalTriggerSimQATask::FillHistograms() {
  // Loop over patches, identify which trigger conditions are met
  DoPatchLoop();

  // Loop over clusters
  DoClusterLoop();
  return kTRUE;
}

/**
 * Create a new instance, adds it to analysis manager
 * \return a pointer to the new instance of the class
 */
AliEmcalTriggerSimQATask * AliEmcalTriggerSimQATask::AddTaskEmcalTriggerSimQA() {
  TString ClusterContName = "caloClusters";

  // Get the pointer to the existing analysis manager via the static access method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
		::Error("AliEmcalTriggerSimQATask", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    ::Error("AliEmcalTriggerSimQATask", "This task requires an input event handler");
    return 0;
  }

	// Init the task and set settings
	TString taskName("AliEmcalTriggerSimQATask");
	AliEmcalTriggerSimQATask * eTask = new AliEmcalTriggerSimQATask(taskName);
  eTask->AddClusterContainer(ClusterContName);

	mgr->AddTask(eTask);
	// Create containers for input/output
	AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
	mgr->ConnectInput(eTask, 0, cinput1);

  TString commonoutput(Form("%s", AliAnalysisManager::GetCommonFileName()));

  TString contOutName(Form("%s_histos", taskName.Data()));
  mgr->ConnectOutput(eTask, 1, mgr->CreateContainer(contOutName, TList::Class(), AliAnalysisManager::kOutputContainer, commonoutput.Data()));

  return eTask;
}

Int_t AliEmcalTriggerSimQATask::GeneratePatchTriggerBits(AliEMCALTriggerPatchInfo * patch) {

  Int_t fPatchTriggerBits = 0x0;

  Bool_t fIsDCal = patch->IsDCalPHOS();
  // In Trigger Types array, DCal Triggers are 1 more than the EMCal Trigger type
  Int_t fTriggerInt = 0;

  if (patch->IsLevel0()) {
    fTriggerInt = kEL0 + fIsDCal;
    fPatchTriggerBits |= 0x1 << fTriggerInt;
  }
  if (patch->IsGammaHigh()) {
    fTriggerInt = kEG1 + fIsDCal;
    fPatchTriggerBits |= 0x1 << fTriggerInt;
  }
  if (patch->IsGammaLow()) {
    fTriggerInt = kEG2 + fIsDCal;
    fPatchTriggerBits |= 0x1 << fTriggerInt;
  }
  if (patch->IsJetHigh()) {
    fTriggerInt = kEJ1 + fIsDCal;
    fPatchTriggerBits |= 0x1 << fTriggerInt;
  }
  if (patch->IsJetLow()) {
    fTriggerInt = kEJ2 + fIsDCal;
    fPatchTriggerBits |= 0x1 << fTriggerInt;
  }
  return fPatchTriggerBits;
}

void AliEmcalTriggerSimQATask::DoPatchLoop() {
  fEventTriggerBits = 0x0; // Reset
  fTriggerPatches = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTriggerPatchesName));

  if (!fTriggerPatches) return;

  Int_t nPatches = fTriggerPatches->GetEntries();

  for (Int_t i = 0; i < nPatches; i++) {
    AliEMCALTriggerPatchInfo* patch = static_cast<AliEMCALTriggerPatchInfo*>(fTriggerPatches->At(i));
    if (!patch) continue;
    if (patch->GetADCAmp() < fMinAmplitude) continue;
    Bool_t fIsDCal = patch->IsDCalPHOS();
    // In Trigger Types array, DCal Triggers are 1 more than the EMCal Trigger type
    Int_t fTriggerInt = 0;
    if (patch->IsLevel0()) {
      fTriggerInt = kEL0 + fIsDCal;
      fEventTriggerBits |= 0x1 << fTriggerInt;
      FillPatchHistograms(patch,fTriggerInt);
    }
    if (patch->IsGammaHigh()) {
      fTriggerInt = kEG1 + fIsDCal;
      fEventTriggerBits |= 0x1 << fTriggerInt;
      FillPatchHistograms(patch,fTriggerInt);
    }
    if (patch->IsGammaLow()) {
      fTriggerInt = kEG2 + fIsDCal;
      fEventTriggerBits |= 0x1 << fTriggerInt;
      FillPatchHistograms(patch,fTriggerInt);
    }
    if (patch->IsJetHigh()) {
      fTriggerInt = kEJ1 + fIsDCal;
      fEventTriggerBits |= 0x1 << fTriggerInt;
      FillPatchHistograms(patch,fTriggerInt);
    }
    if (patch->IsJetLow()) {
      fTriggerInt = kEJ2 + fIsDCal;
      fEventTriggerBits |= 0x1 << fTriggerInt;
      FillPatchHistograms(patch,fTriggerInt);
    }
  }
}

void AliEmcalTriggerSimQATask::FillPatchHistograms(AliEMCALTriggerPatchInfo * patch, Int_t i) {
  fHistManager.FillTH2(Form("fHistPatchColEnergy_Trig_%s",fTriggerNames[i+1].Data()),patch->GetColStart(),patch->GetPatchE());
  fHistManager.FillTH2(Form("fHistPatchColOnlineADCAmp_Trig_%s",fTriggerNames[i+1].Data()),patch->GetColStart(),patch->GetADCAmp());
  fHistManager.FillTH2(Form("fHistPatchRowEnergy_Trig_%s",fTriggerNames[i+1].Data()),patch->GetRowStart(),patch->GetPatchE());
  fHistManager.FillTH2(Form("fHistPatchRowOnlineADCAmp_Trig_%s",fTriggerNames[i+1].Data()),patch->GetRowStart(),patch->GetADCAmp());
  fHistManager.FillTH1(Form("fHistPatchEnergy_Trig_%s",fTriggerNames[i+1].Data()),patch->GetPatchE());
  fHistManager.FillTH1(Form("fHistPatchOnlineADCAmp_Trig_%s",fTriggerNames[i+1].Data()),patch->GetADCAmp());
  fHistManager.FillTH2(Form("fHistPatchColRow_Trig_%s",fTriggerNames[i+1].Data()),patch->GetColStart(),patch->GetRowStart());
  fHistManager.FillTH2(Form("fHistPatchEtaPhiGeo_Trig_%s",fTriggerNames[i+1].Data()),patch->GetEtaGeo(),patch->GetPhiGeo());
}


void AliEmcalTriggerSimQATask::DoClusterLoop() {
  TString histname;

  AliClusterContainer * clusters = GetClusterContainer(0);
  if (!clusters) {
    AliError("Cluster Container Not Found");
    return ;
  }
  Int_t nClusters = clusters->GetNClusters();
  fHistManager.FillTH1("NClusters",nClusters);

  // Cycle over clusters
  for (Int_t i = 0; i < nClusters; i++) {
    AliVCluster * cluster = (AliVCluster *) clusters->GetAcceptCluster(i);
    if (!cluster) continue;
    MatchClusterToPatches(cluster);
    for (Int_t j = 0; j < kNTriggerTypes; j++) {
      // Check if this event had this trigger
      if (fTriggerTypes[j] < 0) {
        // This is Minimum Bias, so accept all events
      }
      else if (!(fEventTriggerBits & (0x1 << fTriggerTypes[j]))) {
        continue;
      }

      // Get Cluster vector
      TLorentzVector vCluster;
      clusters->GetMomentum(vCluster,cluster);
      Double_t fPhi = vCluster.Phi();
      if (fPhi < 0) fPhi+=2*TMath::Pi();

      // Split Cluster spectra into EMCal, DCal
      Bool_t isDCal = (fPhi > 4.); // Lazy check
      if (isDCal) {
        histname = TString::Format("fHistDCalClusEnergy_Trig_%s",fTriggerNames[j].Data());
      } else {
        histname = TString::Format("fHistEMCalClusEnergy_Trig_%s",fTriggerNames[j].Data());
      }
      fHistManager.FillTH1(histname,cluster->E());
//      histname = TString::Format("fHistClusEnergy_Trig_%s",fTriggerNames[j].Data());
//      fHistManager.FillTH1(histname,cluster->E());
      histname = TString::Format("fHistClusEtaPhi_Trig_%s",fTriggerNames[j].Data());
      fHistManager.FillTH1(histname,vCluster.Eta(),fPhi);
    }
  }
}

void AliEmcalTriggerSimQATask::MatchClusterToPatches(AliVCluster * cluster) {
  if (!fTriggerPatches) return;
  if (!fGeom) return; // Can't run this without geometry
  AliClusterContainer * clusters = GetClusterContainer(0);
  if (!clusters) {
    AliError("Cluster Container Not Found");
    return ;
  }

  Int_t nPatches = fTriggerPatches->GetEntries();
  for (Int_t i = 0; i < nPatches; i++) {
    AliEMCALTriggerPatchInfo* patch = static_cast<AliEMCALTriggerPatchInfo*>(fTriggerPatches->At(i));
    if (!patch) continue;
    if (patch->GetADCAmp() < fMinAmplitude) continue;
    if (cluster->E() < fMinClusterEnergy) continue;

    TLorentzVector vCluster; // Cluster Vector
    clusters->GetMomentum(vCluster,cluster);

    // Cluster, Patch centers
    Float_t ClusterPhi = vCluster.Phi();
    Float_t ClusterEta = vCluster.Eta();
    Float_t PatchPhi   = patch->GetPhiGeo();
    Float_t PatchEta   = patch->GetEtaGeo();

    Float_t DeltaPhi = fmod(ClusterPhi - PatchPhi,2*TMath::Pi());
    if (DeltaPhi >= TMath::Pi()) DeltaPhi -= 2*TMath::Pi();
    Float_t DeltaEta = ClusterEta - PatchEta;

    // Getting TRU info based on patch corner
    Int_t PatchRow = patch->GetRowStart();
    Int_t PatchCol = patch->GetColStart();

    Int_t fAbsId = -1;
    if (!fGeom->GetAbsFastORIndexFromPositionInEMCAL(PatchCol,PatchRow, fAbsId)) {
      AliError("Could not find FastOR index from position");
    }
    Int_t iTRU = -1;
    Int_t iADC = -1;
    if (!fGeom->GetTRUFromAbsFastORIndex(fAbsId, iTRU, iADC)) {
      AliError(Form("Mapping mismatch: could not find TRU from FastOR Id %d",fAbsId));
    }

    Int_t iTriggerClass = GeneratePatchTriggerBits(patch);

    if (iTriggerClass != 0) {
      for (Int_t i = 1; i < kNTriggerTypes; i++) {
        if (!(iTriggerClass & (0x1 << fTriggerTypes[i]))) continue;

        // Fill These histograms for all pairs:
        TString histName = TString::Format("fHistClusPatchDEtaTRU_Trig_%s",fTriggerNames[i].Data());
        fHistManager.FillTH2(histName,DeltaEta,iTRU);

        histName = TString::Format("fHistClusPatchDPhiTRU_Trig_%s",fTriggerNames[i].Data());
        fHistManager.FillTH2(histName,DeltaPhi,iTRU);

        histName = TString::Format("fHistClusPhiClusPatchDPhi_Trig_%s",fTriggerNames[i].Data());
        fHistManager.FillTH2(histName,ClusterPhi,DeltaPhi);
        histName = TString::Format("fHistClusEtaClusPatchDPhi_Trig_%s",fTriggerNames[i].Data());
        fHistManager.FillTH2(histName,ClusterEta,DeltaPhi);
        histName = TString::Format("fHistClusPhiClusPatchDEta_Trig_%s",fTriggerNames[i].Data());
        fHistManager.FillTH2(histName,ClusterPhi,DeltaEta);
        histName = TString::Format("fHistClusEtaClusPatchDEta_Trig_%s",fTriggerNames[i].Data());
        fHistManager.FillTH2(histName,ClusterEta,DeltaEta);
      }
    }

    // Match Method 1
    // Check Delta R between patch cluster

    // Match Method 2
    // Check if Cluster Center is in Patch
    Float_t PatchPhiMin = patch->GetPhiMin();
    Float_t PatchPhiMax = patch->GetPhiMax();
    Float_t PatchEtaMin = patch->GetEtaMin();
    Float_t PatchEtaMax = patch->GetEtaMax();

    if ((PatchPhiMin <= ClusterPhi) && (ClusterPhi <= PatchPhiMax)) {
      if ((PatchEtaMin <= ClusterEta) && (ClusterEta <= PatchEtaMax)) {
        // Patch Cluster Matched
        if (iTriggerClass != 0) {
          for (Int_t i = 1; i < kNTriggerTypes; i++) {
            if (!(iTriggerClass & (0x1 << fTriggerTypes[i]))) continue;
            TString sHistName = TString::Format("fHistClusEnergyPatchAmpTRU_Trig_%s",fTriggerNames[i].Data());
            fHistManager.FillTH3(sHistName,cluster->E(),patch->GetADCAmp(),iTRU);
          }
        }
      }
    }
  }
}


