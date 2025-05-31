#include "AliAnalysisManager.h"
#include "AliAnalysisTaskCVEUtil.h"
#include "TError.h"
#include <TString.h>

AliAnalysisTaskCVEUtil* AddTaskCVEUtil(
    bool isMC = false,
    TString period = "LHC18q",
    TString uniqueID = "default"
)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Fatal("AddTaskCVEUtil.C", "No analysis manager to connect to.");
        return nullptr;
    }
    if (!mgr->GetInputEventHandler()) {
        Fatal("AddTaskCVEUtil.C", "No input event handler.");
        return nullptr;
    }

    AliAnalysisTaskCVEUtil* task = new AliAnalysisTaskCVEUtil("TaskCVEUtil");
    if(!task) {
      Fatal("AddTaskCVEUtil.C", "Failed to create task.");
      return nullptr;
    }
    task -> SetMC(isMC);
    task -> SetPeriod(period);

    AliAnalysisDataContainer* inputContainer = mgr->GetCommonInputContainer();
    if (!inputContainer) {
        Fatal("AddTaskCVEUtil.C", "No input container.");
        return nullptr;
    }
    mgr->ConnectInput(task,0,inputContainer);
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    AliAnalysisDataContainer* outputContainer = mgr->CreateContainer(Form("ResultsList_%s", uniqueID.Data()), TList::Class(),
                                                    AliAnalysisManager::kOutputContainer,
                                                    Form("%s:%s", outputFileName.Data(), uniqueID.Data()));
    if (!outputContainer) {
        Fatal("AddTaskCVEUtil.C", "Failed to create output container.");
        return nullptr;
    }
    mgr->ConnectOutput(task,1,outputContainer);
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
