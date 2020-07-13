#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskNdimFitsTPCsplines.h"
#include <TString.h>
#include <TList.h>
#include <TTree.h>
#endif


//_____________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskNdimFitsTPCsplines *addTask_NdimFitsTPCsplines()  {
    
    //Get Analysis Manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return 0x0;
   
    //Get Input Event Handler
    if (!mgr->GetInputEventHandler()) return 0x0;
     
    //File Name
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":NdimFitsTPCsplines";

    //Analysis Task
    AliAnalysisTaskNdimFitsTPCsplines *task = new AliAnalysisTaskNdimFitsTPCsplines ("task_nDimTPCsplines");
    task -> SelectCollisionCandidates (AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral);
    mgr -> AddTask(task);
    mgr -> ConnectInput (task,0,mgr->GetCommonInputContainer());
    mgr -> ConnectOutput(task,1,mgr->CreateContainer("InputAN", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr -> ConnectOutput(task,2,mgr->CreateContainer("InputQA", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    return task;
}
//_____________________________________________________________________________________________________________________________________________________________

