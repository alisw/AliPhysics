AliTaskConfigOCDB* AddTaskConfigOCDB(const char *path, Int_t run=0) 
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
       ::Error("AddTaskConfigOCDB", "No analysis manager to connect to.");
       return NULL;
    }   
    TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    if (inputDataType != "ESD") {
       ::Error("AddTaskCDBconnect", "Can only run with ESD input handler");
       return NULL;
    }   
    
    AliTaskConfigOCDB *task= new AliTaskConfigOCDB("CDBconnect", path, run);
    mgr->AddTask(task);
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();    
    mgr->ConnectInput(task,  0, cinput1);
    return task;
}   
