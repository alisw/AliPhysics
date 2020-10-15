AliAnalysisTaskNuclei* AddHighMultNucleiTask(TString name = "name", ULong64_t triggerMask = AliVEvent::kHighMultV0)
{
    // get the manager via the static access member. since it's static, you don't need
    // to create an instance of the class here to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    // get the input event handler, again via a static method. 
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":AntiHe3HM";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskNuclei* task = new AliAnalysisTaskNuclei(name.Data());   
    if(!task) return 0x0;
	//Add task settings here
	task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
	task->SetFilterBit(256);
	task->SetLowPCut(0.1);
	task->SetHighPCut(1e30);
	task->SetEtaCut(0.8);
	task->SetMinNITSCl(1);
	task->SetMaxDCAxyPreCut(1.5);
	task->SetMaxDCAxyFinal(1.5);
	task->SetMaxDCAz(1.5);
	//set PID cuts #### Legacy code not used in analysis anymore ###
	task->SetMaxTPCnSigma(3.0);
	task->SetUseTOFPidCut(kFALSE);//kTRUE or kFALSE for use of TOF
	task->SetMaxTOFnSigma(3.0);
	// momentum p from which a hit/cut in TOF is required
	task->SetMomForTOFanaProt(0.7);
	task->SetMomForTOFanaDeut(1.4);
	task->SetAnalyseAllParticles(kTRUE);
	
	AliESDtrackCuts* trackCuts = new AliESDtrackCuts();
	trackCuts->SetPtRange(0.15, 1e30);
	trackCuts->SetEtaRange(-0.8, 0.8);
	/*
	//TPC
	trackCuts->SetRequireTPCRefit(kTRUE);
	trackCuts->SetMinNCrossedRowsTPC(110);
	trackCuts->SetMaxChi2PerClusterTPC(3.0);
	trackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
	trackCuts->SetMaxFractionSharedTPCClusters(0.1);
	trackCuts->SetMaxChi2TPCConstrainedGlobal(25.);
	//ITS
	trackCuts->SetRequireITSRefit(kTRUE);
	trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
	trackCuts->SetMinNClustersITS(3);//This seems very high?

	//primary selection
	trackCuts->SetAcceptKinkDaughters(kFALSE);
	trackCuts->SetDCAToVertex2D(kFALSE);
	trackCuts->SetRequireSigmaToVertex(kFALSE);
	trackCuts->SetMaxDCAToVertexZ(10.0);
	trackCuts->SetMaxDCAToVertexXY(10.0);
	*/

	task->SetTrackCuts(trackCuts);

    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("MyOutputContainerHighMult", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
