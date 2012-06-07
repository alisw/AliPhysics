/**
 * @file   AddTaskForwardFlow.C
 * @author Alexander Hansen alexander.hansen@cern.ch 
 * @date   Wed Sep 07 12:14:17 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */
/** 
 * Add FMD event plane task to train 
 * 
 * @param mc 
 *
 * @ingroup pwglf_forward_flow
 */
void AddTaskFMDEventPlane(Bool_t mc = kFALSE)
{
  // --- Get the analysis manager ------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) Fatal("","No analysis manager to connect to.");

  // --- Check that we have an AOD input handler ---------------------
  UShort_t aodInput = 0;
  if (!(aodInput = AliForwardUtil::CheckForAOD())) 
    Fatal("","Cannot proceed without and AOD handler");
  if (aodInput == 2 &&
      !AliForwardUtil::CheckForTask("AliForwardMultiplicityBase")) 
    Fatal("","The relevant task wasn't added to the train");


  // --- Make the event plane task -----------------------------------
  AliFMDEventPlaneTask* task = new AliFMDEventPlaneTask("FMDEventPlane");
  task->GetEventPlaneFinder().SetUsePhiWeights(false);
  if (mc) task->SetMCInput(true);
  mgr->AddTask(task);

  // --- Create containers for output --------------------------------
  AliAnalysisDataContainer* sums = 
    mgr->CreateContainer("FMDEventPlaneSums", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer* output = 
    mgr->CreateContainer("FMDEventPlaneResults", TList::Class(), 
			 AliAnalysisManager::kParamContainer, 
			 AliAnalysisManager::GetCommonFileName());
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, sums);
  mgr->ConnectOutput(task, 2, output);

  return;
}
/*
 * EOF
 */
