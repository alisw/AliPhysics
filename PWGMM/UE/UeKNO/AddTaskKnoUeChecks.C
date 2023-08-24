///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskKnoUeChecks Macro to run on grids               //
//                                                               //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
using namespace std;

AliAnalysisTaskKnoUeChecks* AddTaskKnoUeChecks(
    const Char_t *taskname = "McKnoUe", Bool_t fUseMC = kFALSE, Bool_t fIsMCclosureTest = kFALSE, Double_t minPt = 0.5, Double_t PtLmin = 5.0, Double_t PtLmax = 40.0)
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
    AliAnalysisTaskKnoUeChecks* taskUE = new AliAnalysisTaskKnoUeChecks("taskKno");
    if(!taskUE) return 0x0;
    taskUE->SetUseMC(fUseMC);
    taskUE->SetMCclosureTest(fIsMCclosureTest);
    taskUE->SetPtMin(minPt);         //0.5  GeV/c
    taskUE->SetLeadingPtMin(PtLmin); //5.0  GeV/c
    taskUE->SetLeadingPtMax(PtLmax); //40.0 GeV/c
    

    // add your task to the manager
    mgr->AddTask(taskUE);

    mgr->ConnectInput(taskUE,0,mgr->GetCommonInputContainer());
    if (fUseMC) {
        mgr->ConnectOutput(taskUE,1,mgr->CreateContainer(Form("cList%s_%1.2f",taskname,minPt), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s_MC", AliAnalysisManager::GetCommonFileName(),taskname)));
    }
    else {
        mgr->ConnectOutput(taskUE,1,mgr->CreateContainer(Form("cList%s_%1.2f",taskname,minPt), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname)));
    }

    return taskUE;
}


