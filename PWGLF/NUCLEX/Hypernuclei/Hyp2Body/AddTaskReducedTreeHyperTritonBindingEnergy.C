#if !defined (__CINT__) || defined (__CLING__)
#include <Rtypes.h>
#include <TString.h>
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskReducedTreeHypertritonBindingEnergy.h"
#include "AliPID.h"
#endif

//_______________________________________________________________________________________________________________________________________________________
AliAnalysisTask *AddTaskReducedTreeHyperTritonBindingEnergy (const char *centralityTrigger = "MinBias", Double_t CentralityMin = 0.0, Double_t CentralityMax = 90.0)  {
    
    
    //Get Analysis Manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        printf("ERROR: Analysis Manager Not Found!!!\n");
        return NULL;
    }
    
    //Retrieve Input Event Handler
    if (!mgr->GetInputEventHandler()) {
        printf("ERROR: Input Event Handler Not Found!!!\n");
        return NULL;
    }
    
    
    //Input List Names
    const char *nameListData = Form ("InputHypertritonTreeBindingEnergyData_%s",centralityTrigger);
    const char *nameListQA   = Form ("InputHypertritonTreeBindingEnergyQA_%s",centralityTrigger);
    const char *taskHypName  = Form ("TaskHyperTritonTree_%s",centralityTrigger);
    
    
    //Analysis Task
    AliAnalysisTaskReducedTreeHypertritonBindingEnergy *task = new AliAnalysisTaskReducedTreeHypertritonBindingEnergy  (taskHypName);
    if (CentralityMin==0.0  && CentralityMax==90.0) task -> SelectCollisionCandidates (AliVEvent::kINT7);
    if (CentralityMin==0.0  && CentralityMax==10.0) task -> SelectCollisionCandidates (AliVEvent::kCentral);
    if (CentralityMin==30.0 && CentralityMax==50.0) task -> SelectCollisionCandidates (AliVEvent::kSemiCentral);
    task -> AliAnalysisTaskReducedTreeHypertritonBindingEnergy::SetCentrality (CentralityMin,CentralityMax);
    mgr -> AddTask(task);
    AliAnalysisDataContainer *outputData = mgr -> CreateContainer(nameListData,TList::Class(),AliAnalysisManager::kOutputContainer,AliAnalysisManager::GetCommonFileName());
    AliAnalysisDataContainer *outputQA   = mgr -> CreateContainer(nameListQA,TList::Class(),AliAnalysisManager::kOutputContainer,AliAnalysisManager::GetCommonFileName());
    mgr -> ConnectInput  (task,0,mgr->GetCommonInputContainer());
    mgr -> ConnectOutput (task,1,outputData);
    mgr -> ConnectOutput (task,2,outputQA);
    
    return task;
}
//_______________________________________________________________________________________________________________________________________________________
