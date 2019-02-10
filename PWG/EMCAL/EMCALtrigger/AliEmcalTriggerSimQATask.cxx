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
fPtBinWidth(0.5),
fMaxPt(120),
fEventTriggerBits(0x0),
fHistManager(name)
{
  SetMakeGeneralHistograms(kTRUE);

  // Initializing array in CINT-compatible way
  EventEMCALTriggerType_t fTriggerTypesValues[kNTriggerTypes] = {kMB, kEL0, kEG1, kEG2, kEJ1, kEJ2};
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
  Double_t fMinPhi = 0;
  Double_t fMaxPhi = 2*TMath::Pi();
  Int_t nEtaBins = 96;
  Int_t nPhiBins = 208;

  Double_t fMinCol = 0;
  Double_t fMaxCol = 48;
  Double_t fMinRow = 0;
  Double_t fMaxRow = 104;
  Int_t nCol = 48;
  Int_t nRow = 104;

  Double_t fMaxADCValue = 1024;
  Double_t nADCBins = 512;


  // Hist for counting clusters
  fHistManager.CreateTH1("NClusters","NClusters;N_{cluster}; counts",300,0,300);

  // loop over trigger types
  for (Int_t i = 0; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistClusEnergy_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#it{E}_{cluster} (GeV); counts";
    fHistManager.CreateTH1(histname.Data(),title.Data(),nPtBins,0,fMaxPt);
  }

  for (Int_t i = 0; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistClusEtaPhi_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#eta_{cluster};#phi_{cluster};counts";
    fHistManager.CreateTH2(histname.Data(),title.Data(),nEtaBins,fMinEta,fMaxEta,nPhiBins,fMinPhi,fMaxPhi);
  }

  // Only make patch histgrams for L1 patches:
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

void AliEmcalTriggerSimQATask::DoPatchLoop() {
  fEventTriggerBits = 0x0; // Reset
  fTriggerPatches = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTriggerPatchesName));

  if (!fTriggerPatches) return;

  Int_t nPatches = fTriggerPatches->GetEntries();

  for (Int_t i = 0; i < nPatches; i++) {
    AliEMCALTriggerPatchInfo* patch = static_cast<AliEMCALTriggerPatchInfo*>(fTriggerPatches->At(i));
    if (!patch) continue;
    if (patch->GetADCAmp() < fMinAmplitude) continue;
    if (patch->IsLevel0()) {
      fEventTriggerBits |= 0x1 << kL0;
      FillPatchHistograms(patch,1);
    }
    if (patch->IsGammaLow()) {
      fEventTriggerBits |= 0x1 << kEG1;
      FillPatchHistograms(patch,2);
    }
    if (patch->IsGammaHigh()) {
      fEventTriggerBits |= 0x1 << kEG2;
      FillPatchHistograms(patch,3);
    }
    if (patch->IsJetLow()) {
      fEventTriggerBits |= 0x1 << kEJ1;
      FillPatchHistograms(patch,4);
    }
    if (patch->IsJetHigh()) {
      fEventTriggerBits |= 0x1 << kEJ2;
      FillPatchHistograms(patch,5);
    }
  }
}

void AliEmcalTriggerSimQATask::FillPatchHistograms(AliEMCALTriggerPatchInfo * patch, Int_t i) {
  fHistManager.FillTH1(Form("fHistPatchEnergy_Trig_%s",fTriggerNames[i].Data()),patch->GetPatchE());
  fHistManager.FillTH1(Form("fHistPatchOnlineADCAmp_Trig_%s",fTriggerNames[i].Data()),patch->GetADCAmp());
  fHistManager.FillTH2(Form("fHistPatchColRow_Trig_%s",fTriggerNames[i].Data()),patch->GetColStart(),patch->GetRowStart());
  fHistManager.FillTH2(Form("fHistPatchEtaPhiGeo_Trig_%s",fTriggerNames[i].Data()),patch->GetEtaGeo(),patch->GetPhiGeo());
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

      histname = TString::Format("fHistClusEnergy_Trig_%s",fTriggerNames[j].Data());
      fHistManager.FillTH1(histname,cluster->E());
      histname = TString::Format("fHistClusEtaPhi_Trig_%s",fTriggerNames[j].Data());
      fHistManager.FillTH1(histname,vCluster.Eta(),fPhi);
    }
  }
}

