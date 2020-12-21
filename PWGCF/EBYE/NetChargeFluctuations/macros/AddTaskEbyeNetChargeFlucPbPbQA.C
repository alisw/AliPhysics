
AliAnalysisTaskEbyeChargeFlucPbPbQA* AddTaskEbyeNetChargeFlucPbPbQA(
                                                                TString suffixName = "",
                                                                Bool_t MCthere = kFALSE,
                                                                Int_t zvtxcut = 10,
                                                                Int_t trackBit = 768,
                                                                Int_t MaxTPCclus = 100.,
                                                                Int_t TPCNclus = 80.,
                                                                Double_t DCAxyCut =2.4,
                                                                Double_t DCAzCut =3.2,
                                                                Bool_t MCFill = kFALSE

                                                                )


{ // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    
    // now we create an instance for task
    AliAnalysisTaskEbyeChargeFlucPbPbQA* task = new AliAnalysisTaskEbyeChargeFlucPbPbQA("");
    if(!task) return 0x0;
    
    // add your task to the manager
    mgr->AddTask(task);
    task->SetAnalysisType(MCthere);
    task->Setzvtxcut(zvtxcut);
    task->SettrackBit(trackBit);
    task->SetMaxTPCCluster(MaxTPCclus);
    task->SetNTPCCluster(TPCNclus);
    task->SetDCACut(DCAxyCut,DCAzCut);
    task->SetMCTreeFill(MCFill);
    
    // Create containers for input/output
    TString finDirname      =   "_INT7";
    
    finDirname              +=   suffixName.Data();
    
    TString fileName = mgr->GetCommonFileName(); // new add
    fileName += ":EbyEtask";
    fileName += finDirname.Data();
    printf("container name: %s\n", fileName.Data());
    
    
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("EbyE_%s", suffixName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
//    mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("fTree_%s", suffixName.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("fTreedata_%s", suffixName.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,3,mgr->CreateContainer(Form("fTreeMCrec_%s", suffixName.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,4,mgr->CreateContainer(Form("fTreeMCgen_%s", suffixName.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
    return task;
    
}
