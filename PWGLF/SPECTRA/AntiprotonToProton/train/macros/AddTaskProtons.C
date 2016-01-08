AliAnalysisTaskProton *AddTaskProtons(const Char_t * addname="", Bool_t lCollidingSystems=kTRUE,Bool_t lDelegateSelection = kTRUE,Bool_t fixDCA=kTRUE,Bool_t onDCAz=kFALSE){

	//--- get the current analysis manager ---//
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTask_Proton", "No analysis manager found.");
		return 0;
	}

	// -- check for ESD and MC ---//
	Bool_t hasESD=kFALSE,hasMC=kFALSE;
	AliESDInputHandler *esdH = static_cast<AliESDInputHandler*>(mgr->GetInputEventHandler());
	if (esdH) hasESD=kTRUE; cout<<"ESD: "<<hasESD<<endl;
	if(!hasESD) return NULL;

	//========= Add task to the ANALYSIS manager =====
	TString taskname = "ProtonSysCheck";
	TString outname  = "ProtonSysCheck";
		
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

AliAnalysisTaskProton *taskcheck = new AliAnalysisTaskProton(taskname);

//--- analysis modes ---//
taskcheck->SetCollidingSystems(lCollidingSystems);
taskcheck->SetUsePhysicsSelection(lDelegateSelection);
taskcheck->SetMultiplicityMode(lCollidingSystems);
taskcheck->SetCentralityWindow(0,100);

//--- Acceptance  cuts ---//
taskcheck->SetPtSpace(6, 0.45, 1.05);
taskcheck->SetMaxPrimaryVtxPosZ(10.);
taskcheck->SetMinTPCClusters(80);
taskcheck->SetMaxChi2PerTPCCluster(3.5);
taskcheck->SetMaxChi2PerITSCluster(36.);
taskcheck->SetMinITSClusters(2);

//--- Primaries cuts ---//
if (onDCAz) taskcheck->SetMaxDCAZ(2.);
if (fixDCA) taskcheck->SetMaxDCAXY(0.2);
else taskcheck->SetPtDependentDCAxy(7,0.0026,0.0050,1.01);

//--- PID ---//
taskcheck->SetPIDMode(AliAnalysisTaskProton::kSigma,3,3,0.7);

mgr->AddTask(taskcheck);

//Char_t outFileName[256]={0};
//sprintf(outFileName,"%s_Tree.root",taskname);

	//================================================
	//              data containers
	//================================================
	//            find input container
	AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(outname,
                                                            TList::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            AliAnalysisManager::GetCommonFileName());
outname  +="_QA";
	 AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(outname,
                                                            TList::Class(),
							    AliAnalysisManager::kOutputContainer,
                                                            AliAnalysisManager::GetCommonFileName());
//outname  +="_Syst";
// AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(outname,
//                                                           TList::Class(),
//							    AliAnalysisManager::kOutputContainer,
// 
//                                                           Form("%s:mmeres", AliAnalysisManager::GetCommonFileName()));

//--- connect containers ---//
  mgr->ConnectInput(taskcheck, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskcheck,1,coutput1);
  mgr->ConnectOutput(taskcheck,2,coutput2);
//  mgr->ConnectOutput(taskcheck,3,coutput3);

// set debug level
mgr->SetDebugLevel(0);

	return taskcheck;
}
