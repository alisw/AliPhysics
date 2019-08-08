

//_______________________________________________________________________________________________________________________________________________________
AliAnalysisTask *AddTaskReducedTreeHyperTritonBindingEnergy (Double_t CentralityMin = 0.0, Double_t CentralityMax = 10.0)  {
    
    
    //Get Analysis Manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        printf("ERROR: Analysis Manager Not Found!!!\n");
        return NULL;
    }
    
    //Retrieve Input Event Handler
    if (!mgr->GetInputEventHandler()) {
        printf("ERROR: Input Event Handler Not Found!!!\n");
        return NULL;
    }
    
    //Get File Name
    TString Filename = mgr -> GetCommonFileName();
    
    //Analysis Task
    AliAnalysisTaskReducedTreeHypertritonBindingEnergy *task = new AliAnalysisTaskReducedTreeHypertritonBindingEnergy  ("taskHyperTritonTree");
    if (CentralityMin==0.0  && CentralityMax==90.0) task -> SelectCollisionCandidates (AliVEvent::kINT7);
    if (CentralityMin==0.0  && CentralityMax==10.0) task -> SelectCollisionCandidates (AliVEvent::kCentral);
    if (CentralityMin==30.0 && CentralityMax==50.0) task -> SelectCollisionCandidates (AliVEvent::kSemiCentral);
    task -> AliAnalysisTaskReducedTreeHypertritonBindingEnergy::SetCentrality (CentralityMin,CentralityMax);
    mgr -> AddTask(task);
    AliAnalysisDataContainer *outputData = mgr -> CreateContainer("InputData",TList::Class(),AliAnalysisManager::kOutputContainer,Filename.Data());
    AliAnalysisDataContainer *outputQA   = mgr -> CreateContainer("InputQA",TList::Class(),AliAnalysisManager::kOutputContainer,Filename.Data());
    mgr -> ConnectInput  (task,0,mgr->GetCommonInputContainer());
    mgr -> ConnectOutput (task,1,outputData);
    mgr -> ConnectOutput (task,2,outputQA);
    
    return task;
}
//_______________________________________________________________________________________________________________________________________________________


