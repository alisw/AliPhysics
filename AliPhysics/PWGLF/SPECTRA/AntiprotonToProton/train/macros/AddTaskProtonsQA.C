AliAnalysisTaskProtonQA *AddTaskProtonsQA(const Char_t * addname="", Bool_t lCollidingSystems=kTRUE,Bool_t lDelegateSelection = kTRUE,Bool_t fixDCA=kTRUE,Bool_t onDCAz=kFALSE){

	//--- get the current analysis manager ---//
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTask_Proton__QA", "No analysis manager found.");
		return 0;
	}

	// -- check for ESD and MC ---//
	Bool_t hasESD=kFALSE,hasMC=kFALSE;
	AliESDInputHandler *esdH = static_cast<AliESDInputHandler*>(mgr->GetInputEventHandler());
	if (esdH) hasESD=kTRUE; cout<<"ESD: "<<hasESD<<endl;
	if(!hasESD) return NULL;
  	AliMCEventHandler *mc = new AliMCEventHandler();
  	mgr->SetMCtruthEventHandler(mc);

	//========= Add task to the ANALYSIS manager =====
	TString taskname = "ProtonSysCheckQA";
	TString outname  = "ProtonSysCheckQA";
		
		if(lCollidingSystems) {
      	 taskname += "_cent";
       	 outname  += "_cent";
   		}
	else {
		taskname +="_MB";
		outname  +="_MB"; 
		}

		if(fixDCA) {
       		taskname += "_fixDCA";
        	outname  += "_fixDCA";
   		}
	else {
		taskname +="_PtDep";
		outname  +="_PtDep"; 
		}

	if(onDCAz) {
	 taskname += "_DCAz";
       	 outname  += "_DCAz";
	}
 taskname += "_";
 outname  += "_";

 taskname += addname;
 outname  += addname;

AliAnalysisTaskProtonQA *taskcheckQA = new AliAnalysisTaskProtonQA(taskname);

//--- analysis modes ---//
taskcheckQA->SetCollidingSystems(lCollidingSystems);
taskcheckQA->SetUsePhysicsSelection(lDelegateSelection);
taskcheckQA->SetMultiplicityMode(lCollidingSystems);
taskcheckQA->SetCentralityWindow(0,100);

//--- Acceptance  cuts ---//
taskcheckQA->SetPtSpace(6, 0.45, 1.05);
taskcheckQA->SetMaxPrimaryVtxPosZ(10.);
taskcheckQA->SetMinTPCClusters(80);
taskcheckQA->SetMaxChi2PerTPCCluster(3.5);
taskcheckQA->SetMaxChi2PerITSCluster(36.);
taskcheckQA->SetMinITSClusters(2);

//--- Primaries cuts ---//
if (onDCAz) taskcheckQA->SetMaxDCAZ(1.);
if (fixDCA) taskcheckQA->SetMaxDCAXY(10.);
else taskcheckQA->SetPtDependentDCAxy(7,0.0026,0.0050,1.01);

//--- PID ---//
taskcheckQA->SetPIDMode(AliAnalysisTaskProton::kSigma,3,3,0.7);

mgr->AddTask(taskcheckQA);

	//================================================
	//              data containers
	//================================================
	//            find input container
	AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(outname,
                                                            TList::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            AliAnalysisManager::GetCommonFileName());

//--- connect containers ---//
  mgr->ConnectInput(taskcheckQA, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskcheckQA,1,coutput1);

// set debug level
mgr->SetDebugLevel(0);

	return taskcheckQA;
}
