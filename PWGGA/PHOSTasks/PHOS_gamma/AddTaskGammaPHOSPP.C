#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskGammaPHOSPP.h"
#include "AliCaloPhoton.h"
#include <TString.h>
#include <TList.h>
#include <THashList.h>
#include <AliTender.h>
#endif

AliAnalysisTaskGammaPHOSPP* AddTaskGammaPHOSPP( Bool_t isMC = kFALSE,
                                                TString name = "GammaAnalysis")
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
    //fileName += ":MyTask";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskGammaPHOSPP* task = new AliAnalysisTaskGammaPHOSPP(name.Data());   
    if(!task) return 0x0;

    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task, 0 , mgr->GetCommonInputContainer());
    // same for the output

  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("Data", THashList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:GammaPHOSPP",AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput2 =
    mgr->CreateContainer("Data2", THashList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:GammaPHOSPP",AliAnalysisManager::GetCommonFileName()));

    mgr->ConnectOutput(task, 1, coutput1);
    mgr->ConnectOutput(task, 2, coutput2);
    return task;
}
