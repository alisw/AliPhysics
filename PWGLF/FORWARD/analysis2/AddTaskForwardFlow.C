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
void AddTaskForwardFlow(TString type = "", 
                        Bool_t mc = kFALSE,
                        TString addFlow = "",
                        Int_t addFType = 0,
                        Int_t addFOrder = 0)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddFMDFlowTask", "No analysis manager to connect to.");
    return NULL;
  }   

  AliAODInputHandler* aodInput = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
   
  Bool_t aod = kFALSE;
  if (aodInput) aod = kTRUE;
  if (!aod) {
    Error("AddTaskForwardFlow", "No analysis manager to connect to.");
    return NULL;
  }

  // --- Check which harmonics to calculate --- //

  Bool_t v1 = kTRUE;
  Bool_t v2 = kTRUE;
  Bool_t v3 = kTRUE;
  Bool_t v4 = kTRUE;
  Bool_t v5 = kTRUE;
  Bool_t v6 = kTRUE;

  if (type.Length() > 0) {
    if (!type.Contains("1")) v1 = kFALSE;
    if (!type.Contains("2")) v2 = kFALSE;
    if (!type.Contains("3")) v3 = kFALSE;
    if (!type.Contains("4")) v4 = kFALSE;
    if (!type.Contains("5")) v5 = kFALSE;
    if (!type.Contains("6")) v6 = kFALSE;
  }

  // --- Create containers for output --- //
  AliAnalysisDataContainer* sums = 
    mgr->CreateContainer("FlowQCSums", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer* output = 
    mgr->CreateContainer("FlowQCResults", TList::Class(), 
			 AliAnalysisManager::kParamContainer, 
			 AliAnalysisManager::GetCommonFileName());

  // --- For the selected flow tasks the input and output is set --- //
  
  AliForwardFlowTaskQC* task = 0;
  if (mc) 
    task = new AliForwardMCFlowTaskQC("QCumulants");
  else
    task = new AliForwardFlowTaskQC("QCumulants");
  mgr->AddTask(task); 

  // Set which harmonics to do
  task->SetDoHarmonics(v1, v2, v3, v4, v5, v6);
  // Set non-default axis for vertices
  TAxis* a = new TAxis(20, -10, 10);
  task->SetVertexAxis(a);
  // Set debug flag
  task->SetDebugLevel(0);
  // Set up adding flow to MC input
  if (mc) {
    AliForwardMCFlowTaskQC* mcTask = 
      static_cast<AliForwardMCFlowTaskQC*>task;
    mcTask->SetUseImpactParameter(true);
    mcTask->AddFlow(addFlow);
    mcTask->AddFlowType(addFType);
    mcTask->AddFlowOrder(addFOrder);
  }
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, sums);
  mgr->ConnectOutput(task, 2, output);

  return;
}
