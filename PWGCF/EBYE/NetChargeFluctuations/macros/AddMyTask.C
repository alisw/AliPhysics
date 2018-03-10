
#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEbyeCharge.h"
#include "AliMultSelection.h"
#include <TString.h>
#include <TList.h>
#endif

//AliAnalysisTaskEbyeCharge* AddMyTask(Int_t MCthere, TString name = "name")
AliAnalysisTaskEbyeCharge* AddMyTask(Int_t MCthere)
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
    
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    if(type.Contains("ESD"))
    {
        ::Error("AddTaskMyTask", "This task requires to run on AOD");
        return NULL;
    }
    
    // by default, a file is open for writing. here, we get the filename
    //   TString fileName = AliAnalysisManager::GetCommonFileName();
    //  fileName += ":MyTask";      // create a subfolder in the file
    
    // in your macro:
    gROOT->ProcessLine(".L $ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask= AddTaskPhysicsSelection(kFALSE,0);  //kFALSE for data, kTRUE for MC
    
    // now we create an instance of your task
    AliAnalysisTaskEbyeCharge* task = new AliAnalysisTaskEbyeCharge("analysismytask");
    if(!task) return 0x0;
    // task->SelectCollisionCandidates(AliVEvent::kINT7);
    
    //  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    //  AliMultSelectionTask * task = AddTaskMultSelection();
    
    if(MCthere){
        task->SetAnalysisType(MCthere);
        /*  const Int_t tmpCentbins  = 10;
         Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80};
         task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);*/
    }
    
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("fOutputList", TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer("fTree", TTree::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(task,3,mgr->CreateContainer("fTreeMCrec", TTree::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(task,4,mgr->CreateContainer("fTreeMCgen", TTree::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName()));
    
    
    
    return task;
    
}
