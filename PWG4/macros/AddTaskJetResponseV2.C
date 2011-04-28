AliAnalysisTaskJetResponseV2* AddTaskJetResponseV2(Char_t* type = "clusters", Char_t* jf = "FASTKT", Float_t radius = 0.4, UInt_t filterMask = 256 , Float_t ptTrackMin = 0.15, Int_t iBack = 1, Int_t eventClassMin = 1, Int_t eventClassMax = 4, Int_t reactionPlaneBin = -1){

  Printf("adding task jet response\n");

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if(!mgr){
	::Error("AddTaskJetResponseV2", "No analysis manager to connect to.");
	return NULL;
    }
    if(!mgr->GetInputEventHandler()){
        ::Error("AddTaskJetResponseV2", "This task requires an input event handler.");
	return NULL;
    }

    TString suffix = "";
    suffix += Form("_%s", jf);
    suffix += Form("%02d", (int)((radius+0.01)*10.));
    suffix += Form("_B0");                                // no background subtraction for extra-only
    suffix += Form("_Filter%05d", filterMask);
    suffix += Form("_Cut%05d", (int)((1000.*ptTrackMin)));
    if(type=="clusters") suffix += Form("_Skip00");

    TString suffix2 = "";
    suffix2 += Form("_%s", jf);
    suffix2 += Form("%02d", (int)((radius+0.01)*10.));
    suffix2 += Form("_B%d", iBack);
    suffix2 += Form("_Filter%05d", filterMask);
    suffix2 += Form("_Cut%05d", (int)((1000.*ptTrackMin)));
    if(type=="clusters") suffix2 += Form("_Skip00");

    AliAnalysisTaskJetResponseV2 *task = new AliAnalysisTaskJetResponseV2(Form("JetResponseV2%s", suffix2.Data()));

    TString branch1 = Form("%sAODextraonly%s",type, suffix.Data());
    Printf("Branch1: %s",branch1.Data());

    TString branch2 = Form("%sAODextra%s",type, suffix2.Data());
    Printf("Branch2: %s",branch2.Data());
 
    task->SetBranchNames(branch1,branch2);
    task->SetOfflineTrgMask(AliVEvent::kMB);

    task->SetEvtClassMin(eventClassMin);
    task->SetEvtClassMax(eventClassMax);
    task->SetCentMin(0.);
    task->SetCentMax(100.);
    //task->SetJetDeltaEta(0.2);
    //task->SetJetDeltaPhi(0.2);

    task->SetReactionPlaneBin(reactionPlaneBin);
   
    task->SetJetPtMin(0.);   // min jet pt is implicit a cut on delta pT!!
    //task->SetJetEtaMin(-.4); // reduced tracking acceptance in
    //task->SetJetEtaMax(.4);  // LHC10h pass1

    //task->SetNMatchJets(1); // leading jets only


    mgr->AddTask(task);


    if(reactionPlaneBin>=0){
        suffix2 += Form("_rp%d", reactionPlaneBin);
    }

    AliAnalysisDataContainer *coutputJetResponseV2 = mgr->CreateContainer(
         Form("jetresponseV2_%s%s", type,suffix2.Data()), TList::Class(), AliAnalysisManager::kOutputContainer,
         Form("%s:PWG4_JetResponseV2_%s%s", AliAnalysisManager::GetCommonFileName(), type, suffix2.Data()));

    mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
    mgr->ConnectOutput(task, 1, coutputJetResponseV2);

    return task;
}
