/**
 * @file   AddTaskForwardFlowQC.C
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
 * @param maxMom        Max moment to do 
 * @param useEtaGap     Whehter to use @f$\eta@f$ gaps
 * @param mc            Monte-carlo input
 * @param outlierCutFMD Cut to remove events with outliers 
 * @param outlierCutSPD Cut to remove events with outliers 
 * @param etaGap        Size of @f$\eta@f$ gap
 * @param useTPCForRef  Use TPC tracks for reference flow
 * @param useCent       Whether to use centrality or impact parameter for MC 
 * @param useMCVtx      Whether to use vertex info from MC header
 * @param satVtx        Use satellite interactions 
 * @param addFlow       Afterburn what (MC only)
 * @param addFType      Afterburner parameterization
 * @param addFOrder     Afterburder order 
 *
 * @ingroup pwglf_forward_flow
 */
void AddTaskForwardFlowQC(Int_t    maxMom        = 5,
                          TString  fwdDet        = "FMD",
                          Bool_t   useEtaGap     = kFALSE,
                          Bool_t   use3cor       = kFALSE,
                          Bool_t   mc            = kFALSE,
			  Double_t outlierCutFMD = 4.0, 
			  Double_t outlierCutSPD = 4.0,
			  Double_t etaGap        = 2.0,
			  Bool_t   useTPCForRef  = kFALSE,
			  Bool_t   useCent       = kFALSE,
			  Bool_t   useMCVtx      = kFALSE,
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
  fwdDet.ToUpper();
  
  // --- Set flow flags --------------------------------------------
  if (useEtaGap && use3cor) 
    Fatal("", "You're doing it wrong! Cannot do both eta-gap and 3-sub");
  UShort_t flags = AliForwardFlowTaskQC::kStdQC|AliForwardFlowTaskQC::kSymEta;
  if (useEtaGap)              flags = AliForwardFlowTaskQC::kEtaGap;
  if (use3cor)                flags = AliForwardFlowTaskQC::k3Cor;
  if (satVtx)                 flags |= AliForwardFlowTaskQC::kSatVtx;
  if (fwdDet.Contains("FMD")) flags |= AliForwardFlowTaskQC::kNUAcorr;
  if      (fwdDet.Contains("FMD"))   flags |= AliForwardFlowTaskQC::kFMD;
  else if (fwdDet.Contains("VZERO")) flags |= AliForwardFlowTaskQC::kVZERO;
  if (useTPCForRef) flags |= AliForwardFlowTaskQC::kTPC;

  const char* name = Form("ForwardFlowQC%s%s", fwdDet.Data(), AliForwardFlowTaskQC::GetQCType(flags, false));
  AliForwardFlowTaskQC* task = 0;
  // --- Set up adding flow to MC input ----------------------------
  if (mc) {
    AliForwardMCFlowTaskQC* mcTask = new AliForwardMCFlowTaskQC(name);
    mcTask->SetUseImpactParameter(!useCent);
    mcTask->SetUseMCHeaderVertex(useMCVtx);
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
//  mgr->SetSkipTerminate(true);
//  task->SelectCollisionCandidates(AliVEvent::kCentral);
  task->SetFlowFlags(flags);
  
  // --- Set eta gap value -----------------------------------------
  if (useEtaGap) {
    if (useTPCForRef) task->SetEtaGapValue(0.4);
    else              task->SetEtaGapValue(etaGap);
  }
  else if (useTPCForRef && fwdDet.Contains("FMD")) task->SetEtaGapValue(0.1);

  // --- Check which harmonics to calculate --------------------------
  task->SetMaxFlowMoment(maxMom);
  
  // --- Set non-default axis for vertices ---------------------------
  TAxis* a = 0;
  if (satVtx) {
    a = new TAxis(6, 93.75, 318.75);
  }
  else 
    a = new TAxis(20, -10, 10);
//    a = new TAxis(10, -5, 5);
  task->SetVertexAxis(a);

  // --- Set non-default axis for centrality -------------------------
//  Double_t cent[] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100 };
  Double_t cent[] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 95, 100 };
  Int_t nBins = sizeof(cent)/sizeof(Double_t) -1;
  TAxis* centAxis = new TAxis(nBins, cent);
  task->SetCentralityAxis(centAxis);

  // --- Set sigma cuts for outliers ---------------------------------
  task->SetDetectorCuts(outlierCutFMD, outlierCutSPD);

  // --- Create containers for output --------------------------------
  const char* sumName = Form("FlowQCSums%s%s", fwdDet.Data(), AliForwardFlowTaskQC::GetQCType(flags, false));
  const char* resName = Form("FlowQCResults%s%s", fwdDet.Data(), AliForwardFlowTaskQC::GetQCType(flags, false));
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
