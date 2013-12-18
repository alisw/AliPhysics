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
 * @param useEtaGap     Whehter to use @f$\eta@f$ gaps
 * @param etaGap        Size of @f$\eta@f$ gap
 * @param useCent       Whether to use centrality 
 * @param mc            Monte-carlo input
 * @param satVtx        Use satellite interactions 
 * @param outlierCutFMD Cut to remove events with outliers 
 * @param outlierCutSPD Cut to remove events with outliers 
 * @param addFlow       Afterburn what (MC only)
 * @param addFType      Afterburner parameterization
 * @param addFOrder     Afterburder order 
 *
 * @ingroup pwglf_forward_flow
 */
void AddTaskForwardFlow(TString  type          = "234", 
                        Bool_t   useEtaGap     = kFALSE,
                        Bool_t   mc            = kFALSE,
			Double_t outlierCutFMD = 4.0, 
			Double_t outlierCutSPD = 0,
			Double_t etaGap        = 2.0,
			Bool_t   useCent       = kFALSE,
			Bool_t   satVtx        = kFALSE,
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
    Fatal("","The relevant tasks weren't added to the train");

  // --- For the selected flow tasks the input and output is set -----
  const char* name = (useEtaGap ? "ForwardQCumulantsEtaGap", "ForwardQCumulants");
  AliForwardFlowTaskQC* task = 0;
  // --- Set up adding flow to MC input ----------------------------
  if (mc) {
    AliForwardMCFlowTaskQC* mcTask = new AliForwardMCFlowTaskQC(name);
    mcTask->SetUseImpactParameter(!useCent);
    if (addFlow.Data()[0] != '\0') {
      mcTask->AddFlow(addFlow);
      mcTask->AddFlowType(addFType);
      mcTask->AddFlowOrder(addFOrder);
    }
    task = mcTask;
  }
  // --- Or use normal task ----------------------------------------
  else 
    task = new AliForwardFlowTaskQC(name);
  
  mgr->AddTask(task); 

  // --- Set flow flags --------------------------------------------
  UShort_t flags = AliForwardFlowTaskQC::kSymEta;
  if (useEtaGap)           flags |= AliForwardFlowTaskQC::kEtaGap;
  if (satVtx)              flags |= AliForwardFlowTaskQC::kSatVtx;
  if (useEtaGap || satVtx) flags ^= AliForwardFlowTaskQC::kSymEta;
  task->SetFlowFlags(flags);
  
  // --- Set eta gap value -----------------------------------------
  task->SetEtaGapValue(etaGap);

  // --- Check which harmonics to calculate --------------------------
  const char* harm = type.Data();
  Int_t i = 0;
  std::cout << "Type string: " << type.Data();
  std::cout << "\t harm string: " << harm << std::endl;
  while (i < type.Length()) {
    char c = harm[i];
    std::cout << "Adding moment: " << c << std::endl;
    Short_t n = atoi(&c);
    std::cout << "Adding moment: " << n << std::endl;
    task->AddFlowMoment(n);
    i++;
  }

  // --- Set non-default axis for vertices ---------------------------
  TAxis* a = 0;
  if (satVtx) {
    a = new TAxis(6, 93.75, 318.75);
  }
  else 
    a = new TAxis(20, -10, 10);
  task->SetVertexAxis(a);

  // --- Set sigma cuts for outliers ---------------------------------
  task->SetDetectorCuts(outlierCutFMD, outlierCutSPD);

  // --- Create containers for output --------------------------------
  const char* sumName = (useEtaGap ? "FlowQCSumsEtaGap" : "FlowQCSums");
  const char* resName = (useEtaGap ? "FlowQCResultsEtaGap" : "FlowQCResults");
  AliAnalysisDataContainer* sums = 
    mgr->CreateContainer(sumName, TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer* output = 
    mgr->CreateContainer(resName, TList::Class(), 
			 AliAnalysisManager::kParamContainer, 
			 AliAnalysisManager::GetCommonFileName());
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, sums);
  mgr->ConnectOutput(task, 2, output);
}
/*
 * EOF
 */
