#if !defined (__CINT__) || defined (__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskPtEMCalTrigger.h"
#include "AliESDtrackCuts.h"
#include <TList.h>
#endif

AliAnalysisTask* AddTaskPtEMCalTrigger(){
        AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
         
        if (!mgr) {
             ::Error("AddTaskPtEMCalTrigger", "No analysis manager to connect to.");
             return NULL;
        }
         
        if (!mgr->GetInputEventHandler()) {
             ::Error("AddTaskPtEMCalTrigger", "This task requires an input event handler");
             return NULL;
        }
        
        EMCalTriggerPtAnalysis::AliAnalysisTaskPtEMCalTrigger *pttriggertask = new EMCalTriggerPtAnalysis::AliAnalysisTaskPtEMCalTrigger("ptemcaltriggertask");
        pttriggertask->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kEMC7);                          // Select both INT7 or EMC7 triggered events
        mgr->AddTask(pttriggertask);

        // Create charged hadrons pPb standard track cuts
        AliESDtrackCuts *standardTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
        standardTrackCuts->SetName("Standard Track cuts");
        standardTrackCuts->SetMinNCrossedRowsTPC(120);
        standardTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
        pttriggertask->AddTrackCuts(standardTrackCuts);

        // Create hybrid track cuts as used in the jet analysis
        AliESDtrackCuts* hybridTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
        hybridTrackCuts->SetName("Global Hybrid tracks, loose DCA");
        hybridTrackCuts->SetMaxDCAToVertexXY(2.4);
        hybridTrackCuts->SetMaxDCAToVertexZ(3.2);
        hybridTrackCuts->SetDCAToVertex2D(kTRUE);
        hybridTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
        hybridTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
        pttriggertask->AddTrackCuts(hybridTrackCuts);

        AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput = mgr->CreateContainer("results", TList::Class(),    AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
   
        //Connect input/output
        mgr->ConnectInput(pttriggertask, 0, cinput);
        mgr->ConnectOutput(pttriggertask, 1, coutput);
           
        return pttriggertask;
}
