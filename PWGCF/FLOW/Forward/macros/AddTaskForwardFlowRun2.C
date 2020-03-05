/**
 * @file   FTAddMyTask.C
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 *
 * @brief  Add Q-cummulant forward task to train
 *
 *
 * @ingroup pwglf_forward_scripts_tasks
 */
/**
 * @defgroup pwglf_forward_flow Flow
 *
 * Code to deal with flow
 *
 * @ingroup pwglf_forward_topical
 */
/**
 * Add Flow task to train
 *
 * @ingroup pwglf_forward_flow
 */



AliAnalysisTaskSE* AddTaskForwardFlowRun2(bool mc,  bool esd, 
                                          TString nua_file, TString nue_file, 
                                          UInt_t tracktype, TString centrality,
                                          TString sec_file_fwd,
                                          TString sec_file_cent,
                                          TString suffix)
{
  std::cout << "______________________________________________________________________________" << std::endl;

  std::cout << "AddTaskForwardFlowRun2" << std::endl;

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    Fatal("","No analysis manager to connect to.");

  AliForwardFlowRun2Task* task = new AliForwardFlowRun2Task(suffix + "_task");

  task->fSettings.mc = mc;
  task->fSettings.esd = esd;
  task->fSettings.nua_file = nua_file;
  task->fSettings.nue_file = nue_file;
  task->fSettings.sec_file = sec_file_fwd;
  task->fSettings.sec_cent_file = sec_file_cent;
  std::cout << "Using tracktype = " << tracktype << std::endl;
  if (tracktype == 0){
    task->fSettings.useSPD = kTRUE;
  }
  else{
    if (tracktype == 768){
      task->fSettings.tracktype = AliForwardSettings::kHybrid;
    }
    else if (tracktype == 128){
      task->fSettings.tracktype = AliForwardSettings::kTPCOnly;
    }
    else if (tracktype == 32){
      task->fSettings.tracktype = AliForwardSettings::kGlobal;
    }
    else if (tracktype == 64){
      task->fSettings.tracktype = AliForwardSettings::kGlobalLoose;
    }
    else if (tracktype == 96){
      task->fSettings.tracktype = AliForwardSettings::kGlobalComb;
    }
    else{
      std::cout << "INVALID TRACK TYPE FOR TPC" << std::endl;
    }
  }
  
  task->fSettings.centrality_estimator = centrality; // "V0M";// RefMult08; // "V0M" // "SPDTracklets";

  TString combName = suffix;

  std::cout << "Container name: " << combName << std::endl;
  std::cout << "______________________________________________________________________________" << std::endl;

  task->fSettings.fileName = suffix;
  mgr->AddTask(task);

  AliAnalysisDataContainer *coutput_recon =
  mgr->CreateContainer(suffix,
   AliForwardFlowResultStorage::Class(),
   AliAnalysisManager::kOutputContainer,
   mgr->GetCommonFileName());
  mgr->ConnectOutput(task, 1, coutput_recon);
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());


  AliAnalysisDataContainer* valid = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("event_selection_xchange");
  mgr->ConnectInput(task,1,valid);

  return task;
}
/*
 * EOF
 *
 */
