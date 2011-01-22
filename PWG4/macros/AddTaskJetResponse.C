AliAnalysisTaskJetResponse* AddTaskJetResponse(Char_t* jf = "FASTKT", Float_t radius = 0.4, UInt_t filterMask = 256 , Float_t ptTrackMin = 0.15, Int_t iBack = 1, Int_t eventClassMin = 1, Int_t eventClassMax = 5){

  Printf("adding task jet response\n");

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if(!mgr){
	::Error("AddTaskJetResponse", "No analysis manager to connect ot.");
	return NULL;
    }
    if(!mgr->GetInputEventHandler()){
        ::Error("AddTaskJetResponse", "This task requires an input event handler.");
	return NULL;
    }

    AliAnalysisTaskJetResponse *task = new AliAnalysisTaskJetResponse("JetResponse");

    TString branch1 = "jetsAODextraonly";
    branch1 += Form("_%s", jf);
    branch1 += Form("%02d", (int)((radius+0.01)*10.));
    branch1 += Form("_B%d", iBack);
    branch1 += Form("_Filter%05d", filterMask); 
    branch1 += Form("_Cut%05d", (int)((1000.*kPtTrackMin)));
    Printf("Branch1: %s",branch1.Data());

    TString branch2 = "jetsAODextra";
    branch2 += Form("_%s", jf);
    branch2 += Form("%02d", (int)((radius+0.01)*10.));
    branch2 += Form("_B%d", iBack);
    branch2 += Form("_Filter%05d", filterMask); 
    branch2 += Form("_Cut%05d", (int)((1000.*kPtTrackMin)));
    Printf("Branch2: %s",branch2.Data());
 
    task->SetBranchNames(branch1,branch2);
    //task->SetOfflineTrgMask(AliVEvent::kMB);

    task->SetEvtClassMin(eventClassMin);
    task->SetEvtClassMax(eventClassMax);
    task->SetCentMin(0.);
    task->SetCentMax(100.);


    mgr->AddTask(task);

    AliAnalysisDataContainer *coutputJetResponse = mgr->CreateContainer(
         Form("jetresponse%s%02d_B%d_Filter%05d_Cut%05d", jf, (int)((radius+0.01)*10.), iBack, filterMask, (int)((1000.*kPtTrackMin))), TList::Class(), AliAnalysisManager::kOutputContainer,
         Form("%s:PWG4_JetResponse_%s%02d_B%d_Filter%05d_Cut%05d", AliAnalysisManager::GetCommonFileName(), jf, (int)((radius+0.01)*10.), iBack, filterMask, (int)((1000.*kPtTrackMin))));

    mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
    mgr->ConnectOutput(task, 1, coutputJetResponse);

    return task;
}
