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
  EventEMCALTriggerType_t fTriggerTypesValues[kNTriggerTypes] = {kNTr, kEL0, kEG1, kEG2, kEJ1, kEJ2};
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

  // Hist for counting clusters
  fHistManager.CreateTH1("NClusters","NClusters;N_{cluster}; counts",300,0,300);

  // loop over trigger types
  for (Int_t i = 0; i < kNTriggerTypes; i++) {
    histname = TString::Format("fHistClusEnergy_Trig_%s",fTriggerNames[i].Data());
    title = histname + ";#it{E}_{cluster} (GeV); counts";
    fHistManager.CreateTH1(histname.Data(),title.Data(),nPtBins,0,fMaxPt);
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
    if (patch->IsLevel0()) fEventTriggerBits |= 0x1 << kL0;
    if (patch->IsGammaLow()) fEventTriggerBits |= 0x1 << kEG1;
    if (patch->IsGammaHigh()) fEventTriggerBits |= 0x1 << kEG2;
    if (patch->IsJetLow()) fEventTriggerBits |= 0x1 << kEJ1;
    if (patch->IsJetHigh()) fEventTriggerBits |= 0x1 << kEJ2;
  }
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
        if (fEventTriggerBits != 0) continue; // Classify no trigger as mininum bias
      }
      else if (!(fEventTriggerBits & (0x1 << fTriggerTypes[j]))) {
        continue;
      }

      histname = TString::Format("fHistClusEnergy_Trig_%s",fTriggerNames[j].Data());
      fHistManager.FillTH1(histname,cluster->E());
    }
  }
}

