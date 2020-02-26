#if !defined(__CINT__) || defined(__MAKE_CINT__) || defined(__CLING__)

#include <TString.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliVEventHandler.h"
#include "AliITSPidParams.h"
#include "AliAnalysisTaskSEITSsaSpectra.h"

#endif

AliAnalysisTaskSEITSsaSpectra* AddTaskITSsaSpectra(Int_t    pidMethod, // 0:kNSigCut, 1:kMeanCut, 2:kLanGaus
                                                   Bool_t   isMC         = kFALSE, //
                                                   Bool_t   defPriors    = kTRUE,
                                                   Bool_t   optNtuple    = kFALSE,
                                                   Bool_t   useUnfolding = kFALSE,
                                                   const char* unfpath   = "",
                                                   const char* suffix    = "")
{
  // Creates, configures and attaches to the train the task for pi, K , p spectra
  // with ITS standalone tracks
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    ::Error("AddTaskITSsaBayes", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if(!mgr->GetInputEventHandler()) {
    ::Error("AddTaskITSsaBayes", "This task requires an input event handler");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD")) {
    ::Error("AddUserTask", "This task requires to run on ESD");
    return NULL;
  }

  const char* pidName[3] = {"_NSigCut", "_MeanCut", "_LanGauMP"};
  TString kContSuffix(pidName[pidMethod]);

  // Create and configure the task
  AliAnalysisTaskSEITSsaSpectra* taskits = new AliAnalysisTaskSEITSsaSpectra(defPriors,optNtuple);

  taskits->SetPidTech(pidMethod);
  taskits->SetIsMC(isMC);
  if (pidMethod==2) {
    AliITSPidParams* aliParams = new AliITSPidParams(isMC);
    taskits->SetITSPidParams(aliParams);
  }

  const int nMultBins=14;
  double mult[nMultBins+1] = {-5.f,0.f,1.f,5.f,10.f,15.f,20.f,30.f,40.f,50.f,60.f,70.f,80.f,90.f,100.f};
  taskits->SetCentBins(nMultBins, mult);

  const int nPtBins=24;
  double ptBins[nPtBins+1] = {
   0.00f,0.05f,0.08f,0.10f,0.12f,0.14f,0.16f,0.18f,0.20f,0.25f,
   0.30f,0.35f,0.40f,0.45f,0.50f,0.55f,0.60f,0.65f,0.70f,0.75f,
   0.80f,0.85f,0.90f,0.95f,1.00f
  };
  taskits->SetPtBins(nPtBins, ptBins);

  //set unfolded probability matrices (if enabled)
  if(useUnfolding){
    taskits->SetUseUnfolding(useUnfolding);
    taskits->SetUnfoldingProb(unfpath);
  }

  taskits->Initialization();

  kContSuffix += suffix;
  mgr->AddTask(taskits);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput(taskits,  0, mgr->GetCommonInputContainer());

  TString outputFileName = AliAnalysisManager::GetCommonFileName();

  TString kMainContName("cListSpectraPiKaPrITSsa");
  kMainContName += kContSuffix;

  AliAnalysisDataContainer* kMainDataCont = NULL;
  kMainDataCont = mgr->CreateContainer(kMainContName.Data(),
                                       TList::Class(),
                                       AliAnalysisManager::kOutputContainer,
                                       outputFileName);
  mgr->ConnectOutput(taskits, 1, kMainDataCont);

  TString kDCAcutContName("cListDCAcutFunction");
  kDCAcutContName += kContSuffix;

  AliAnalysisDataContainer* kDCAcutDataCont = NULL;
  kDCAcutDataCont = mgr->CreateContainer(kDCAcutContName.Data(),
                                         TList::Class(),
                                         AliAnalysisManager::kParamContainer,
                                         outputFileName);
  mgr->ConnectOutput(taskits, 2, kDCAcutDataCont);

  TString kNtupleContName("cListTreeInfo");
  kNtupleContName += kContSuffix;

  AliAnalysisDataContainer* kNtupleDataCont = NULL;
  if (optNtuple){
        kNtupleDataCont = mgr->CreateContainer(kNtupleContName.Data(),
                                           TList::Class(),
                                           AliAnalysisManager::kOutputContainer,
                                           outputFileName);
    mgr->ConnectOutput(taskits, 3, kNtupleDataCont);
  }
  return taskits;
}
