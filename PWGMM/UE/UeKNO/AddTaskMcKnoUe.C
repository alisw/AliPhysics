///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskMcKnoUe Macro to run on grids               //
//                                                               //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskMcKnoUe* AddTaskMcKnoUe(const Char_t* taskname="McKnoUe", Bool_t  useMC  = kFALSE, Bool_t performMCclosuretest = kFALSE, Bool_t isPythia=kFALSE,Bool_t IsppData=kFALSE,Bool_t IspPbData=kTRUE, Double_t minpT=0.5,Double_t PtLmin = 5.0, Double_t PtLmax = 40.0)
{
    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function


    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    // get the input event handler this handler is part of the managing system and feeds events to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

    // now you create an instance of your task
    AliAnalysisTaskMcKnoUe* taskUE = new AliAnalysisTaskMcKnoUe("taskKno");
    if(!taskUE) return 0x0;
    taskUE->SetUseMC(useMC);
    taskUE->SetMCclosureTest(performMCclosuretest);
    taskUE->SetParametrizationEfficiency(isPythia);
    taskUE->SetParametrizationEfficiencyppdata(IsppData);
    taskUE->SetParametrizationEfficiencypPbdata(IspPbData);
    // add your task to the manager
    taskUE->SetPtMin(minpT);
    taskUE->SetLeadingPtMin(PtLmin);
    taskUE->SetLeadingPtMax(PtLmax);
    mgr->AddTask(taskUE);


    mgr->ConnectInput(taskUE,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskUE,1,mgr->CreateContainer(Form("cList%s_%1.2f",taskname,minpT), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname)));

    return taskUE;
}
