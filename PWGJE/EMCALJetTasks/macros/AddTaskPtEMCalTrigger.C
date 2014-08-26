#if !defined (__CINT__) || defined (__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskPtEMCalTrigger.h"
#include "AliESDtrackCuts.h"
#include <TList.h>
#include <TString.h>
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
        //pttriggertask->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kEMC7);                          // Select both INT7 or EMC7 triggered events
        pttriggertask->SelectCollisionCandidates(AliVEvent::kAny);
        mgr->AddTask(pttriggertask);
        pttriggertask->SetPtRange(2., 100.);

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

        TString containerName = mgr->GetCommonFileName();
        containerName += ":PtEMCalTriggerTask";
        printf("container name: %s\n", containerName.Data());

        AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput = mgr->CreateContainer("results", TList::Class(),    AliAnalysisManager::kOutputContainer, containerName.Data());
   
        //Connect input/output
        mgr->ConnectInput(pttriggertask, 0, cinput);
        mgr->ConnectOutput(pttriggertask, 1, coutput);
           
        return pttriggertask;
}
