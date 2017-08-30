AliAnalysisTaskSEITSsaSpectra* AddTaskITSsaSpectra(Int_t    pidMethod, // 0:kNSigCut, 1:kMeanCut, 2:kLanGaus
                                                   Bool_t   isMC       = kFALSE, //
                                                   Bool_t   optNtuple  = kFALSE,
                                                   const char* suffix  = "")
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

  char* pidName[3] = {"_NSigCut", "_MeanCut", "_LanGauMP"};
  TString kContSuffix(pidName[pidMethod]);

  // Create and configure the task
  AliAnalysisTaskSEITSsaSpectra* taskits = new AliAnalysisTaskSEITSsaSpectra();
  taskits->SetIsMC(isMC);
  taskits->SetFillNtuple(optNtuple);
  taskits->SetPidTech(pidMethod);

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

  if(optNtuple) {
    TString kNtupleContName("cListTreeInfo");
    kNtupleContName += kContSuffix;

    AliAnalysisDataContainer* kNtupleDataCont;
    kNtupleDataCont = mgr->CreateContainer(kContSuffix.Data(),
                                           TTree::Class(),
                                           AliAnalysisManager::kOutputContainer,
                                           outputFileName);
    mgr->ConnectOutput(taskits, 3, kNtupleDataCont);
  }

  return taskits;
}
