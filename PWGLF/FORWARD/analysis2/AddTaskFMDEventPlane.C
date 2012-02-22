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
 * @defgroup pwglf_forward_flow Flow 
 * @ingroup pwglf_forward_topical
 */
/** 
 * Add Flow task to train 
 * 
 * @param type 
 * @param etabins 
 * @param mc 
 * @param addFlow 
 * @param addFType 
 * @param addFOrder 
 *
 * @ingroup pwglf_forward_flow
 */
void AddTaskFMDEventPlane(Bool_t mc = kFALSE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFMDEventPlane", "No analysis manager to connect to.");
    return NULL;
  }   

  AliAODInputHandler* aodInput = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
   
  Bool_t aod = kFALSE;
  if (aodInput) aod = kTRUE;
  if (!aod) {
    Error("AddTaskFMDEventPlane", "No analysis manager to connect to.");
    return NULL;
  }

  // --- Create containers for output --- //
  AliAnalysisDataContainer* sums = 
    mgr->CreateContainer("FMDEventPlaneSums", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer* output = 
    mgr->CreateContainer("FMDEventPlaneResults", TList::Class(), 
			 AliAnalysisManager::kParamContainer, 
			 AliAnalysisManager::GetCommonFileName());

  // --- For the selected flow tasks the input and output is set --- //
  
  AliFMDEventPlaneTask* task = new AliFMDEventPlaneTask("FMDEventPlane");
  task->GetEventPlaneFinder().SetUsePhiWeights(false);

  if (mc) task->SetMCInput(true);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, sums);
  mgr->ConnectOutput(task, 2, output);

  return;
}
