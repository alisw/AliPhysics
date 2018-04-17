AliAnalysisTask *AddTask_slehner_ElectronEfficiency(Bool_t hasITS = kTRUE,
                                                     Int_t trackCut=1,
                                                     Int_t evCut=1,
                                                     TString directoryBaseName = "slehner",
                                                     Char_t* outputFileName="LMEE_output.root",
                                                     Bool_t deactivateTree=kFALSE, // enabling this has priority over 'writeTree'! (needed for LEGO trains)
                                                     TString resolutionfile = "" //Resolution_pp_16l_eta.root
                                                     )
{

	//get the current analysis manager
	AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
	if(!mgr){
		
		Error("AddTask_slehner_ElectronEfficiency", "No analysis manager found.");
		return 0;
	}

	//Base Directory for GRID / LEGO Train  
	TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";

        TString configLMEECutLib("LMEECutLib_slehner.C");
        TString configLMEECutLibPath(configBasePath+configLMEECutLib);
        
	std::cout << "Configpath:  " << configFilePath << std::endl;

        //LOAD CUTLIB
        if(gSystem->Exec(Form("ls %s", configLMEECutLibPath.Data()))==0){

		std::cout << "loading LMEECutLib: " << configLMEECutLibPath.Data() << std::endl;
		gROOT->LoadMacro(configLMEECutLibPath.Data());
	} 
	else{
		std::cout << "LMEECutLib not found: " << configLMEECutLibPath.Data() << std::endl;
		return 0; // if return is not called, the job will fail instead of running wihout this task... (good for local tests, bad for train)
	}
        
	//Do we have an MC handler?
	Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);
	std::cout << "hasMC = " << hasMC << std::endl;


      LMEECutLib* cutlib = new LMEECutLib();

      AliAnalysisTaskMLTreeMaker *task = new AliAnalysisTaskMLTreeMaker(taskname);

      task->SelectCollisionCandidates(AliVEvent::kINT7);
      task->SetTrackCuts(cutlib->GetTrackSelectionAna(trackCut));
      task->SetEventCuts(cutlib->GetEventCuts(evCut));
        
  mgr->AddTask(taskESD);

  TString outputFileName = AliAnalysisManager::GetCommonFileName();

 
  AliAnalysisDataContainer *coutESD = mgr->CreateContainer("output", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(taskESD, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskESD, 1, coutESD);
 
 return taskESD;

}
