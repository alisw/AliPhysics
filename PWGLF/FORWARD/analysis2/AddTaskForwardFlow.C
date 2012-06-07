/**
 * @file   AddTaskForwardFlow.C
 * @author Alexander Hansen alexander.hansen@cern.ch 
 * 
 * @brief  
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
 * @param type          Which moments to do 
 * @param mc            Monte-carlo input
 * @param dispVtx       Use satellite interactions 
 * @param outlierCutFMD Cut to remove events with outliers 
 * @param outlierCutSPD Cut to remove events with outliers 
 * @param addFlow       Afterburn what (MC only)
 * @param addFType      Afterburner parameterization
 * @param addFOrder     Afterburder order 
 *
 * @ingroup pwglf_forward_flow
 */
void AddTaskForwardFlow(TString  type          = "", 
                        Bool_t   mc            = kFALSE,
			Bool_t   dispVtx       = kFALSE,
			Double_t outlierCutFMD = 4.1, 
			Double_t outlierCutSPD = 4.1,
                        TString  addFlow       = "",
                        Int_t    addFType      = 0,
                        Int_t    addFOrder     = 0)
{
  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    Fatal("","No analysis manager to connect to.");

  // --- Check that we have an AOD input handler ---------------------
  UShort_t aodInput = 0;
  if (!(aodInput = AliForwardUtil::CheckForAOD())) 
    Fatal("","Cannot proceed without and AOD handler");
  if (aodInput == 2 && 
      (!AliForwardUtil::CheckForTask("AliForwardMultiplicityBase") ||
       !AliForwardUtil::CheckForTask("AliCentralMultiplicityTask")))
    Fatal("","The relevant tasks wasn't added to the train");

  // --- For the selected flow tasks the input and output is set -----
  AliForwardFlowTaskQC* task = 0;
  if (mc) {
    AliForwardMCFlowTaskQC* mcTask = new AliForwardMCFlowTaskQC("QCumulants");
    // --- Set up adding flow to MC input ----------------------------
    mcTask->SetUseImpactParameter(true);
    mcTask->AddFlow(addFlow);
    mcTask->AddFlowType(addFType);
    mcTask->AddFlowOrder(addFOrder);
    task = mcTask;
  }
  else 
    task = new AliForwardFlowTaskQC("QCumulants");
  mgr->AddTask(task); 

  // --- Check which harmonics to calculate --------------------------
  Bool_t v1 = type.Contains("1");
  Bool_t v2 = type.Contains("2");
  Bool_t v3 = type.Contains("3");
  Bool_t v4 = type.Contains("4");
  Bool_t v5 = type.Contains("5");
  Bool_t v6 = type.Contains("6");
  task->SetDoHarmonics(v1, v2, v3, v4, v5, v6);

  // --- Set non-default axis for vertices ---------------------------
  TAxis* a = 0;
  if (dispVtx) {
    AliForwardFlowTaskQC::fgDispVtx = true;
    a = new TAxis(6, 93.75, 318.75);
  }
  else 
    a = new TAxis(20, -10, 10);
  task->SetVertexAxis(a);

  // --- Set sigma cuts for outliers ---------------------------------
  task->SetDetectorCuts(outlierCutFMD, outlierCutSPD);

  // --- Create containers for output --------------------------------
  AliAnalysisDataContainer* sums = 
    mgr->CreateContainer("FlowQCSums", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer* output = 
    mgr->CreateContainer("FlowQCResults", TList::Class(), 
			 AliAnalysisManager::kParamContainer, 
			 AliAnalysisManager::GetCommonFileName());
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, sums);
  mgr->ConnectOutput(task, 2, output);
}
/*
 * EOF
 */
