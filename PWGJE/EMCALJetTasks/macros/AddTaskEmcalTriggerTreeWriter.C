#if !defined (__CINT__) || defined (__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliESDtrackCuts.h"
#include <TList.h>
#include <TString.h>
#endif

AliAnalysisTask* AddTaskEmcalTriggerTreeWriter(){
        AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
         
        if (!mgr) {
             ::Error("AddTaskPtEMCalTrigger", "No analysis manager to connect to.");
             return NULL;
        }
         
        if (!mgr->GetInputEventHandler()) {
             ::Error("AddTaskPtEMCalTrigger", "This task requires an input event handler");
             return NULL;
        }
        
        AliAnalysisTaskEmcalTriggerTreeWriter *treewriter = new AliAnalysisTaskEmcalTriggerTreeWriter("TriggerTreewriterTask");
        //pttriggertask->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kEMC7);                          // Select both INT7 or EMC7 triggered events
        treewriter->SelectCollisionCandidates(AliVEvent::kAny);
        mgr->AddTask(treewriter);

        AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutputTree = mgr->CreateContainer("EMCalTriggerTree", TTree::Class(),    AliAnalysisManager::kOutputContainer, "EMCalTriggerTree.root");
   
        //Connect input/output
        mgr->ConnectInput(treewriter, 0, cinput);
        mgr->ConnectOutput(treewriter, 1, coutputTree);
           
        return treewriter;
}
