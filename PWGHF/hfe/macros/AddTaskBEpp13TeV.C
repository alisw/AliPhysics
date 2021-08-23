AliAnalysisTaskBEpp13TeV* AddTaskBEpp13TeV( ///-> to run locally
			TString uniqueID        = "",
			bool 	isMC 			= kFALSE,
			int minTPCNclusters =	100,
			int minTPCclustersPID =	80,
			double maxTPCchi2 = 4.,
			double minTPCclusterRatio =	0.6,
			int  minITSNcluster = 4,
			TString itsLayer = "kBoth",
			double eta = 0.8,
			double ptMin = 0.3,
			double ptMax = 30.,
			float dcaxy = 0.1,
			float dcaz = 0.2,
			double tpcPIDlow = -1.,
			double tpcPIDhigh = 3.,
			double tofPID = 3
)           
{
	
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	
	if (!mgr) {
		::Error("AddTaskBEpp13TeV", "No analysis manager to connect to.");
		return NULL;
	}
	
	if (!mgr->GetInputEventHandler()) {
		::Error("AddTaskBEpp13TeV", "This task requires an input event handler");
		return NULL;
	}

	///Task config
	AliAnalysisTaskBEpp13TeV *task = new AliAnalysisTaskBEpp13TeV(uniqueID.Data());
	printf("task ------------------------ %p\n ", task);
	task->SetMinTPCNcls(minTPCNclusters);
	task->SetMinTPCNclsPID(minTPCclustersPID);
	task->SetMaxTPCchi2(maxTPCchi2);
	task->SetMinTPCclsRatio(minTPCclusterRatio);
	task->SetMinITSNcls(minITSNcluster);
	task->SetITSlayer(itsLayer);
	task->SetEtaRange(eta);
	task->SetPtRange(ptMin, ptMax);
	task->SetDCARange(dcaxy, dcaz);

	task->SetPIDCuts(tpcPIDlow, tpcPIDhigh, tofPID);
  
	//Trigger
	task->SelectCollisionCandidates(AliVEvent::kINT7); //Selecting Minumum Bias events (selected randomlly)

	mgr->AddTask(task);
	
	//added to run on train
	TString fileName = AliAnalysisManager::GetCommonFileName();
	fileName += ":ElectroID_";
	fileName += uniqueID;

	//Create containers for input/output -> to run locally
	AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("ccontainer0_%s",uniqueID.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,fileName);
	//Connect input/output
	mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task, 1, coutput);
	
	return task;
}

