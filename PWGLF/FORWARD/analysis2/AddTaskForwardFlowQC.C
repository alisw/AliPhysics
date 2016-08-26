/**
 * @file   AddTaskForwardFlowQC.C
 * @author Alexander Hansen <alexander.hansen@cern.ch>
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
 * @param maxMom        Max moment to do 
 * @param fwdDet        Which forward detector to use (FMD/VZERO)
 * @param useEtaGap     Whehter to use @f$\eta@f$ gaps
 * @param use3cor       Wheter to use the 2 subevent method
 * @param mc            Monte-carlo input
 * @param outlierCutFMD Cut to remove events with outliers 
 * @param outlierCutSPD Cut to remove events with outliers 
 * @param etaGap        Size of @f$\eta@f$ gap
 * @param useTracksForRef  Which tracks to use for reference flow
 * @param useCent       Whether to use centrality or impact parameter for MC 
 * @param ispA          True for pPb 
 * @param useMCVtx      Whether to use vertex info from MC header
 * @param satVtx        Use satellite interactions 
 * @param addFlow       Add flow afterburner (only MC)
 *
 * @ingroup pwglf_forward_flow
 */
void AddTaskForwardFlowQC(Int_t    maxMom          = 5,
                          TString  fwdDet          = "FMD",
                          Bool_t   useEtaGap       = kFALSE,
                          Bool_t   use3cor         = kFALSE,
                          Bool_t   mc              = kFALSE,
			  Double_t outlierCutFMD   = 4.0, 
			  Double_t outlierCutSPD   = 4.0,
			  Double_t etaGap          = 2.0,
			  UInt_t   useTracksForRef = 0,
			  Bool_t   useCent         = kFALSE,
			  Bool_t   ispA            = kFALSE,
			  Bool_t   useMCVtx        = kFALSE,
			  Bool_t   satVtx          = kFALSE,
			  Bool_t   addFlow         = kFALSE)
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
  if (useTracksForRef == 1) flags |= AliForwardFlowTaskQC::kTPC;
  if (useTracksForRef == 2) flags |= AliForwardFlowTaskQC::kHybrid;

  const char* name = Form("ForwardFlowQC%s%s", fwdDet.Data(), AliForwardFlowTaskQC::GetQCType(flags, false));
  AliForwardFlowTaskQC* task = 0;
  // --- Set up adding flow to MC input ----------------------------
  if (mc) {
    AliForwardMCFlowTaskQC* mcTask = new AliForwardMCFlowTaskQC(name);
    mcTask->SetUseImpactParameter(!useCent);
    mcTask->SetUseMCHeaderVertex(useMCVtx);
    if (addFlow) {
      mcTask->SetUseFlowWeights(true);
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
    if (useTracksForRef) task->SetEtaGapValue(0.4);
    else              task->SetEtaGapValue(etaGap);
  }
  else if (useTracksForRef && fwdDet.Contains("FMD")) task->SetEtaGapValue(0.0);

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
  TAxis* centAxis = 0;
  if (ispA) {
    Double_t cent[] = {0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100});
    Int_t nBins = sizeof(cent)/sizeof(Double_t) -1;
    centAxis = new TAxis(nBins, cent);
  } else {
    Double_t cent[] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 95, 100 };
    Int_t nBins = sizeof(cent)/sizeof(Double_t) -1;
    centAxis = new TAxis(nBins, cent);
  }
  task->SetCentralityAxis(centAxis);

  // --- Set sigma cuts for outliers ---------------------------------
  task->SetDetectorCuts(outlierCutFMD, outlierCutSPD);

  // --- Setup track cuts --------------------------------------------
  if (useTracksForRef > 0) {
    AliAnalysisFilter* trackCuts = new AliAnalysisFilter("trackFilter");
    if (!ispA) {
      if (useTracksForRef == 1) { // tpc tracks
	AliESDtrackCuts* tpcTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
	tpcTrackCuts->SetPtRange(0.2, 5.0);
	tpcTrackCuts->SetEtaRange(-0.8, 0.8);
	tpcTrackCuts->SetMinNClustersTPC(70);
	trackCuts->AddCuts(tpcTrackCuts);
	task->SetTrackCuts(trackCuts);
      } else if (useTracksForRef == 2) { // hybrid tracks
	// Basic cuts for global tracks - working for 10h!
	AliESDtrackCuts* baseCuts = new AliESDtrackCuts("base");
	TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
	baseCuts->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
	baseCuts->SetMinNClustersTPC(70);
	baseCuts->SetMaxChi2PerClusterTPC(4);
	baseCuts->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
	baseCuts->SetAcceptKinkDaughters(kFALSE);
	baseCuts->SetRequireTPCRefit(kTRUE);
	baseCuts->SetMaxFractionSharedTPCClusters(0.4);
	// ITS
	baseCuts->SetRequireITSRefit(kTRUE);
	//accept secondaries
	baseCuts->SetMaxDCAToVertexXY(2.4);
	baseCuts->SetMaxDCAToVertexZ(3.2);
	baseCuts->SetDCAToVertex2D(kTRUE);
	//reject fakes
	baseCuts->SetMaxChi2PerClusterITS(36);
	baseCuts->SetMaxChi2TPCConstrainedGlobal(36);
	baseCuts->SetRequireSigmaToVertex(kFALSE);
  //      baseCuts->SetEtaRange(-0.9,0.9);
  //      baseCuts->SetPtRange(0.15, 1E+15.);
	baseCuts->SetEtaRange(-0.8, 0.8);
	baseCuts->SetPtRange(0.2, 5.);
    
	// Good global tracks
	AliESDtrackCuts* globalCuts = baseCuts->Clone("global");
	globalCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
	
	// tracks with vertex constraint
	AliESDtrackCuts* constrainedCuts = baseCuts->Clone("vertexConstrained");
	constrainedCuts->SetRequireITSRefit(kFALSE);

	// Add
	trackCuts->AddCuts(globalCuts);
	trackCuts->AddCuts(constrainedCuts);
	task->SetTrackCuts(trackCuts);
      } // end of hybrid
    } // end of Pb-Pb
    if (ispA) {
      if (useTracksForRef == 1) { // tpc tracks
	AliESDtrackCuts* tpcTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
	tpcTrackCuts->SetPtRange(0.2, 5.0);
	tpcTrackCuts->SetEtaRange(-0.8, 0.8);
	tpcTrackCuts->SetMinNClustersTPC(70);
	trackCuts->AddCuts(tpcTrackCuts);
	task->SetTrackCuts(trackCuts);
      } else if (useTracksForRef == 2) { // hybrid tracks
        AliESDtrackCuts* baseCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
	baseCuts->SetMaxDCAToVertexXY(2.4);
	baseCuts->SetMaxDCAToVertexZ(3.2);
	baseCuts->SetDCAToVertex2D(kTRUE);
	baseCuts->SetMaxChi2TPCConstrainedGlobal(36);
	baseCuts->SetMaxFractionSharedTPCClusters(0.4);
	// extra
	baseCuts->SetEtaRange(-0.8, 0.8);
	baseCuts->SetPtRange(0.2, 5.);
        
        // Global good tracks 
	AliESDtrackCuts* globalCuts = baseCuts->Clone("global");

	// tracks with vertex constraint
	AliESDtrackCuts* constrainedCuts = baseCuts->Clone("vertexConstrained");
	constrainedCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
	constrainedCuts->SetRequireITSRefit(kTRUE);

        // Add
	trackCuts->AddCuts(globalCuts);
	trackCuts->AddCuts(constrainedCuts);
	task->SetTrackCuts(trackCuts);
      } // end of hybrid
    } // end of pA
  } // end of tracks

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
