AliAnalysisTaskProtonAbsorbtion *AddTaskProtonsAbs(const Char_t * addname="", Bool_t lCollidingSystems=kTRUE,Bool_t lDelegateSelection = kTRUE,Bool_t fixDCA=kTRUE,Bool_t onDCAz=kFALSE){

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
	TString taskname = "ProtonSysCheckAbs";
	TString outname  = "ProtonSysCheckAbs";
		
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

AliAnalysisTaskProtonAbsorbtion *taskcheckAbs = new AliAnalysisTaskProtonAbsorbtion(taskname);

//--- analysis modes ---//
taskcheckAbs->SetCollidingSystems(lCollidingSystems);
taskcheckAbs->SetUsePhysicsSelection(lDelegateSelection);
taskcheckAbs->SetMultiplicityMode(lCollidingSystems);
taskcheckAbs->SetCentralityWindow(0,100);

//--- Acceptance  cuts ---//
taskcheckAbs->SetPtSpace(6, 0.45, 1.05);
taskcheckAbs->SetMaxPrimaryVtxPosZ(10.);
taskcheckAbs->SetMinTPCClusters(80);
taskcheckAbs->SetMaxChi2PerTPCCluster(3.5);
taskcheckAbs->SetMaxChi2PerITSCluster(36.);
taskcheckAbs->SetMinITSClusters(2);

//--- Primaries cuts ---//
if (onDCAz) taskcheckAbs->SetMaxDCAZ(1.);
if (fixDCA) taskcheckAbs->SetMaxDCAXY(0.2);
else taskcheckAbs->SetPtDependentDCAxy(7,0.0026,0.0050,1.01);

//--- PID ---//
taskcheckAbs->SetPIDMode(AliAnalysisTaskProton::kSigma,3,3,0.7);

mgr->AddTask(taskcheckAbs);

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
  mgr->ConnectInput(taskcheckAbs, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskcheckAbs,1,coutput1);

// set debug level
mgr->SetDebugLevel(0);

	return taskcheckAbs;
}
