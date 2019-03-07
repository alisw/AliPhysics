///////////////////////////////////////////////////////////////////
//                                                               //            
// AddMyTask                                                     //
// Author: Redmer A. Bertens, Utrecht University, 2012           //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskDiHadCorrelHighPt* AddTaskDiHadCorrelHighPt(TString taskName = "name", Bool_t analysisMC = kFALSE, TString container_name_extension = "",TString fileName_extension = "")
{
    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function
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
    fileName += ":AliAnalysisTaskDiHadCorrelHighPt";      // create a subfolder in the file
    fileName += fileName_extension.Data();
    // now we create an instance of your task
    AliAnalysisTaskDiHadCorrelHighPt* task = new AliAnalysisTaskDiHadCorrelHighPt(taskName.Data(),analysisMC);
    if(!task) return 0x0;
    task->SetPtTrigMin(3);
    task->SetPtAsocMin(1);
    task->SetOStatus(1);
    task->SetCutsCrosscheck(kFALSE);
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    TString container_name = "MyOutputContainer";
    container_name += container_name_extension.Data();
    mgr->ConnectOutput(task,1,mgr->CreateContainer(container_name.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
