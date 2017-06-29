AliAnalysisTaskCmeEse* AddTaskCmeEse(UInt_t filterbit = 768,
                                 Double_t etaCut = 0.8,
                                 Double_t vtxCut = 10.,
                                 Double_t minPt = 0.2,
                                 Double_t maxPt = 5.0,
                                 Int_t noclus = 70,
                                 Bool_t IsLHC10h = kTRUE,
                                 Bool_t IsPileUp = kTRUE,
                                 Bool_t IsRecEff = kFALSE,
                                 Short_t nCentFl = 0,
                                 Bool_t IsQAV0 = kFALSE,
                                 Bool_t IsQATrk = kFALSE,
                                 TString uniqueID = "def")
{

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTaskCmeEse", "No analysis manager to connect to.");
        return NULL;
    }

    
    if (!mgr->GetInputEventHandler()) {
        Error("AddTaskCmeEse", "This task requires an input event handler");
        return NULL;
    }


    AliAnalysisTaskCmeEse* taskD = new AliAnalysisTaskCmeEse(Form("taskCMEESE_%s", uniqueID.Data()));
    
    taskD->SetDebugLevel(3);
    taskD->SetFilterbit(filterbit);
    taskD->SetNoClus(noclus);
    taskD->SetEtaCut(etaCut);
    taskD->SetVtxCut(vtxCut);
    taskD->SetMinPt(minPt);
    taskD->SetMaxPt(maxPt);
    taskD->SetFlagLHC10h(IsLHC10h);
    taskD->SetFlagPileUp(IsPileUp);
    taskD->SetCentFlag(nCentFl);
    taskD->SetFlagRecEff(IsRecEff);
    taskD->SetFlagQAV0(IsQAV0);
    taskD->SetFlagQATrk(IsQATrk);
    
    mgr->AddTask(taskD);

    AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer* cout = mgr->CreateContainer(Form("output_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("cme_psi2V0A.root:%s", uniqueID.Data()));
    mgr->ConnectInput (taskD, 0, cinput);
    mgr->ConnectOutput(taskD, 1, cout);
  
    // Return task pointer at the end
    return taskD;
}
