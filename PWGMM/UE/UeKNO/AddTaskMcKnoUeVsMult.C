///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskMcKnoUe Macro to run on grids               //
//                                                               //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskMcKnoUeVsMult* AddTaskMcKnoUeVsMult(const Char_t* taskname="UeVsMult", Bool_t  useMC  = kFALSE, Bool_t performMCclosuretest = kFALSE, Bool_t isPythia=kFALSE,Bool_t IsppData=kFALSE,Bool_t IspPbData=kFALSE, Double_t minpT=0.5,Double_t PtLmin = 0.5, Double_t PtLmax = 40.0,Double_t V0Amin = 0.0, Double_t V0Amax = 100.0,Double_t EtaCut = 0.8)
{
    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function


    AliAnalysisTaskMcKnoUeVsMult* taskUE = new AliAnalysisTaskMcKnoUeVsMult("taskUeVsMult");
    if(!taskUE) return 0x0;
    
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    // get the input event handler this handler is part of the managing system and feeds events to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

    // now you create an instance of your task
    
    taskUE->SetUseMC(useMC);
    taskUE->SetMCclosureTest(performMCclosuretest);
    taskUE->SetParametrizationEfficiency(isPythia);
    taskUE->SetParametrizationEfficiencyppdata(IsppData);
    taskUE->SetParametrizationEfficiencypPbdata(IspPbData);
    // add your task to the manager
    taskUE->SetPtMin(minpT);
    taskUE->SetLeadingPtMin(PtLmin);
    taskUE->SetLeadingPtMax(PtLmax);
    taskUE->SetV0Amin(V0Amin);
    taskUE->SetV0Amax(V0Amax);
    taskUE->SetEtaCut(EtaCut);
    mgr->AddTask(taskUE);


    mgr->ConnectInput(taskUE,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskUE,1,mgr->CreateContainer(Form("cList%s_%1.2f_%1.1fto%1.1f_%1.1f",taskname,minpT,V0Amin,V0Amax,EtaCut), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname)));

    return taskUE;
}
