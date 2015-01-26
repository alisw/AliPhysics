AliAnalysisTaskEMCalHFEpA *AddTaskEMCalHFEpACorrelation(

                        Bool_t  isMC                    = kFALSE, 
                        Int_t   triggerIndex    = 0, 
                        Int_t   configIndex     = 0, 
                        Int_t   centralityIndex = 0, 
                        Bool_t  isAOD                   = kFALSE,
                        Bool_t isEMCal                  = kFALSE,
                        char * period                   = "b",
                        Int_t EMCalThreshould   = 0
                )
{
        AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
        
        if (!mgr) {
        ::Error("AddTaskEMCalHFEpA", "No analysis manager to connect to.");
        return NULL;
        }
        
        if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskEMCalHFEpA", "This task requires an input event handler");
        return NULL;
        }
        
        //_______________________
        //Config Task
        //gROOT->LoadMacro("ConfigEMCalHFEpACorrelation.C");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/ConfigEMCalHFEpACorrelation.C");
        AliAnalysisTaskEMCalHFEpA *task = ConfigEMCalHFEpACorrelation(isMC,triggerIndex,configIndex,centralityIndex,isAOD,isEMCal,EMCalThreshould);
        
        //_______________________
        //Trigger
        if(!isMC && (period=="d" || period=="e" || period=="f"))
        {
                if(triggerIndex == 0) task->SelectCollisionCandidates(AliVEvent::kINT7);
                if(triggerIndex == 1) task->SelectCollisionCandidates(AliVEvent::kEMC7);
                if(triggerIndex == 2) task->SelectCollisionCandidates(AliVEvent::kEMCEGA);
        }
        else if(!isMC)
        {
                if(triggerIndex == 0) task->SelectCollisionCandidates(AliVEvent::kINT7);
                if(triggerIndex == 1) task->SelectCollisionCandidates(AliVEvent::kEMC7);
                if(triggerIndex == 2) task->SelectCollisionCandidates(AliVEvent::kEMCEGA);
                //if(triggerIndex == 3) task->SelectCollisionCandidates(AliVEvent::kEMC8);
                //if(triggerIndex == 4) task->SelectCollisionCandidates(AliVEvent::kEMCEJE); //Jet Trigger
        }
        
        mgr->AddTask(task);
        
        TString containerName = mgr->GetCommonFileName();
		containerName += ":HFE_EMCal_pPb_elienos";
		containerName += Form("_%d_%d_%d_%d",triggerIndex,configIndex,centralityIndex,EMCalThreshould);
        
        //Create containers for input/output
        AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
        AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("chist_eh_%d_%d_%d_%d",triggerIndex,configIndex,centralityIndex,EMCalThreshould), TList::Class(),    AliAnalysisManager::kOutputContainer, containerName.Data());

        //Connect input/output
        mgr->ConnectInput(task, 0, cinput);
        mgr->ConnectOutput(task, 1, coutput);
        
        return task;
}
