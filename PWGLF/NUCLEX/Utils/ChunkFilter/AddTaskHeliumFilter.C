#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
R__ADD_INCLUDE_PATH($ALICE_ROOT)
#endif


//_________________________________________________________________________________________________________________________________________________
AliAnalysisTaskHeliumFilter *AddTaskHeliumFilter(){


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
    
    //Output File Name
    TString outputFileName = mgr->GetCommonFileName();
    
    //Analysis Task
    AliAnalysisTaskHeliumFilter *task_HeliumFilter = new AliAnalysisTaskHeliumFilter ("task_HeliumFilter");
    mgr -> AddTask(task_HeliumFilter);
    AliAnalysisDataContainer *list = mgr -> CreateContainer("List",TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
    mgr -> ConnectInput (task_HeliumFilter,0,mgr -> GetCommonInputContainer());
    mgr -> ConnectOutput(task_HeliumFilter,1,list);
       
    return task;
}
//_________________________________________________________________________________________________________________________________________________
