//////////////////////////////////////////////////////////////////////
//                                                                  //            
//  AddMyTaskCorPIDTOFQA                                            //
//  Author: Brennan Schaefer, Oak Ridge National Laboratory, 2016   //
//                                                                  //
//////////////////////////////////////////////////////////////////////
#include <stdlib.h> 

class AliAnalysisDataContainer;

AliAnalysisTaskCorPIDTOFQA* AddMyTaskCorPIDTOFQA(TString name = "name")
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
    // fileName += ":MyTask";      // create a subfolder in the file
    // now we create an instance of your task
    // AliAnalysisTaskCorPIDTOFQA* task = new BSchaefer_devel::AliAnalysisTaskCorPIDTOFQA(name.Data());
    AliAnalysisTaskCorPIDTOFQA* task = new AliAnalysisTaskCorPIDTOFQA(name.Data());   
    if(!task) return 0x0;
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    // if(run_mode == 0)    mgr->ConnectInput(task    ,0,mgr->GetCommonInputContainer());
    // if(run_mode == 1)    mgr->ConnectInput(task_PID,0,mgr->GetCommonInputContainer());
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output

    int run_mode = 0;
    run_mode = atoi(name);
    // cout<<endl<<endl<<endl<<run_mode<<endl<<endl<<endl;

    if(run_mode == 0){   mgr->ConnectOutput(task,1,mgr->CreateContainer("cOutput",     TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));    }
    if(run_mode == 1){   mgr->ConnectOutput(task,1,mgr->CreateContainer("cOutput_PID", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));    }
    if(run_mode == 2){   mgr->ConnectOutput(task,1,mgr->CreateContainer("cOutput_TPC", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));    }
    if(run_mode == 3){   mgr->ConnectOutput(task,1,mgr->CreateContainer("cOutput_nCl", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));    }
    if(run_mode == 4){   mgr->ConnectOutput(task,1,mgr->CreateContainer("cOutput_DCA", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));    }
    if(run_mode == 5){   mgr->ConnectOutput(task,1,mgr->CreateContainer("cOutput_sid", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));    }
    if(run_mode == 6){   mgr->ConnectOutput(task,1,mgr->CreateContainer("cOutput_glo", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));    }
    // mgr->ConnectOutput(task,1,mgr->CreateContainer("cOutput",     TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
