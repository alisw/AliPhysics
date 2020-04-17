///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskMcKno Macro to run on grids               //
//                                                               //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskMcKno* AddTaskMcKno(const Char_t* taskname="McKno", Bool_t  useMC  = kTRUE, Bool_t performMCclosuretest = kFALSE, Bool_t IspPb = kFALSE, Double_t minpT=0.5, Double_t PtLmin = 5.0, Double_t PtLmax = 40.0, Double_t V0Mmin = 0.0, Double_t V0Mmax = 100.0, Bool_t IsTPConly = kTRUE)
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
    AliAnalysisTaskMcKno* taskKno = new AliAnalysisTaskMcKno("taskKno");
    if(!taskKno) return 0x0;
    taskKno->SetUseMC(useMC);
    taskKno->SetMCclosureTest(performMCclosuretest);
    // add your task to the manager
    taskKno->SetPtMin(minpT);
    taskKno->SetIspPb(IspPb);
    taskKno->SetLeadingPtMin(PtLmin);
    taskKno->SetLeadingPtMax(PtLmax);
    taskKno->SetV0Mmin(V0Mmin);
    taskKno->SetV0Mmax(V0Mmax);
    taskKno->SetNchTScut(IsTPConly);
    mgr->AddTask(taskKno);

    mgr->ConnectInput(taskKno,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskKno,1,mgr->CreateContainer(Form("cList%s_%1.2f_%1.1f_%1.1f_%1.1f_%1.1f",taskname,minpT,PtLmin,PtLmax,V0Mmin,V0Mmax), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname)));

    return taskKno;
}
