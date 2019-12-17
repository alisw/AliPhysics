///////////////////////////////////////////////////////////////////
//                                                                 //
//            AddTaskMpiUE Macro to run on grids                   //
//  Author: Prabi, Feng Fang, and Paolo Bartalini, CCNU 2016       //
//                                                                 //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskMpiUE* AddTaskMpiUE(Bool_t  useMC  = kFALSE)
{
    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function

const char *name="Results.root";

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    // get the input event handler this handler is part of the managing system and feeds events to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

    // now you create an instance of your task
    AliAnalysisTaskMpiUE* taskUE = new AliAnalysisTaskMpiUE("taskUE");
    if(!taskUE) return 0x0;

    // add your task to the manager
    mgr->AddTask(taskUE);
    taskUE->SetUseMC(useMC);

    //Printf("Set OutputFileName : \n %s\n", name.Data());

    // by default, a file is open for writing. here,  get the desired filename
//    if (mgr->GetMCtruthEventHandler()) name += ":MC";

    // your task needs input: here you connect the manager to your task
    mgr->ConnectInput(taskUE,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(taskUE,1,mgr->CreateContainer("UE", TList::Class(), AliAnalysisManager::kOutputContainer, name));
    // this macro returns a pointer to the task, this will be convenient when you run the analysis in an analysis train on grids
    return taskUE;
}
