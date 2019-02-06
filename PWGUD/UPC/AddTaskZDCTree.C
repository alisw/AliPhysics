//////////////////////////////////////////////////////////////////////////////
//                                                               			//            
// AddTaskZDCTree                                      						//
// Created by  Michal Broz Michal.Broz@cern.ch								//
// Modified by Uliana Dmitrieva uliana.dmitrieva@cern.ch on 11/01/2018      //
//////////////////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskZDCTree* AddTaskZDCTree(TString name = "name")
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":ZDCTreeTask";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskZDCTree* task = new AliAnalysisTaskZDCTree(name.Data());   
    if(!task) return 0x0;
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("fOutput", TList::Class(),AliAnalysisManager::kOutputContainer,fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
