                                                                                                                    
class AliAnalysisDataContainer;
#if !defined (__CINT__) || defined (__CLING__)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#endif
#include "AliAnalysisTaskResonanceDP.h"
AliAnalysisTaskResonanceDP* AddMyTaskDPak(
    const char* runtype="local", 
    const char* taskname="PD13TeV",
    int triggerClass = 1
    )
{
    // get the manager via the static access member. since it's static, you don't need                                                                                                   
    // an instance of the class to call the function                                                                                                                                     
    AliAnalysisManager *mgr = (AliAnalysisManager*) AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskPhiPbPb5TeVRun2", "No analysis manager to connect to ==>");
        return 0x0;
    }

    // get the input event handler, again via a static method.                                                                                                                           
    // this handler is part of the managing system and feeds events to your task                                                                                                                                                                      
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

    // by default, a file is open for writing. here, we get the filename     
    const char* ananame  ="myana";                                                                                                                                                                                                                                                    
    AliAnalysisTaskResonanceDP* task = new AliAnalysisTaskResonanceDP(taskname);
    //if (triggerClass ==0 ) task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
    //if (triggerClass ==1 ) task->SelectCollisionCandidates(AliVEvent::kINT7);
    if (triggerClass ==0 ) task->setTriggerType(AliVEvent::kINT7);
    if (triggerClass ==1 ) task->setTriggerType(AliVEvent::kHighMultV0);

    if(!task) return 0x0;
    // add your task to the manager                                                                                                                                                      
    mgr->AddTask(task);
    
    TString outfilename = "AnalysisResults.root";
    
    // create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("nuclei", TTree::Class(), AliAnalysisManager::kOutputContainer, outfilename.Data());

    // connect input/output
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput1);

return task;
}
