AliAnalysisTaskEMCalHFEpA *AddTaskEMCalHFEpA(

                        Bool_t  isMC            = kFALSE, 
                        Int_t   triggerIndex    = 0, 
                        Int_t   configIndex     = 0, 
                        Int_t   centralityIndex = 0, 
                        Bool_t  isAOD           = kFALSE,
                        Bool_t isEMCal          = kFALSE,
						Bool_t isTrigger 		= kFALSE,
                        char * period           = "b",
                        Int_t EMCalThreshould   = 0,
						Bool_t isTender = kFALSE,
						Int_t   centralityEstimator = 0,
						Bool_t isCentralitySys 		= kFALSE,
						Bool_t isTOFdet 		= kFALSE,
						TString finname = ""
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
        //gROOT->LoadMacro("ConfigEMCalHFEpA.C");
        TString configFile="";
        if( finname.EqualTo("") )
			configFile="$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/ConfigEMCalHFEpA.C";
        else  if(!gSystem->Exec(Form("alien_cp %s  .",finname.Data())))
			configFile=Form("%s/ConfigEMCalHFEpA.C",gSystem->pwd()); 
        else {
			printf("ERROR: couldn't copy file %s/ from grid \n", finname.Data() );
			configFile="$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/ConfigEMCalHFEpA.C";
		} 
        gROOT->LoadMacro(Form("%s",configFile.Data()));
        AliAnalysisTaskEMCalHFEpA *task = ConfigEMCalHFEpA(isMC,triggerIndex,configIndex,centralityIndex,isAOD,isEMCal,isTrigger, EMCalThreshould, isTender, period, centralityEstimator, isCentralitySys, isTOFdet);
        
        //_______________________
        //Trigger
        if(!isMC)
        {
                if(triggerIndex == 0) task->SelectCollisionCandidates(AliVEvent::kINT7);
                if(triggerIndex == 1) task->SelectCollisionCandidates(AliVEvent::kEMC7);
                if(triggerIndex == 2) task->SelectCollisionCandidates(AliVEvent::kEMCEGA);
                
                //if(triggerIndex == 3) task->SelectCollisionCandidates(AliVEvent::kEMC8);
                //if(triggerIndex == 4) task->SelectCollisionCandidates(AliVEvent::kEMCEJE); //Jet Trigger
        }
        
        mgr->AddTask(task);
        
        TString containerName = mgr->GetCommonFileName();
        containerName += ":HFE_EMCal_pPb_cris";
        containerName += Form("_%d_%d_%d_%d",triggerIndex,configIndex,centralityIndex,EMCalThreshould);
        

        
        //Create containers for input/output
        AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	
	
        AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("chist_RpPb_%d_%d_%d_%d",triggerIndex,configIndex,centralityIndex, EMCalThreshould), TList::Class(),    AliAnalysisManager::kOutputContainer, containerName.Data());
	
		
		
        //Connect input/output
        mgr->ConnectInput(task, 0, cinput);
        mgr->ConnectOutput(task, 1, coutput);
        
        return task;
}
