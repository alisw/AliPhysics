AliAnalysisTask *AddTask_slehner_ElectronEfficiency(Bool_t hasITS = kTRUE,
                                                     Double_t CentMin = -2,
                                                     Double_t CentMax = 102,
                                                     TString directoryBaseName = "slehner",
                                                     Bool_t getFromAlien=kFALSE,
                                                     TString cFileName = "Config_slehner_ElectronEfficiency.C",
                                                     Char_t* outputFileName="LMEE_Eff_output.root",
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


	TString configFilePath(configBasePath+cFileName);
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
        
        //LOAD CONFIG        
	if(gSystem->Exec(Form("ls %s", configFilePath.Data()))==0){

		std::cout << "loading config: " << configFilePath.Data() << std::endl;
		gROOT->LoadMacro(configFilePath.Data());
	} 
	else{
		std::cout << "config not found: " << configFilePath.Data() << std::endl;
		return 0; // if return is not called, the job will fail instead of running wihout this task... (good for local tests, bad for train)
	}
        

	//Do we have an MC handler?
	Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);
	std::cout << "hasMC = " << hasMC << std::endl;

	// Electron efficiency task
	AliAnalysisTaskElectronEfficiency *task = new AliAnalysisTaskElectronEfficiency("slehner_ElectronEfficiency");
	std::cout << "task created: " << task->GetName() << std::endl;

//	
  mgr->AddTask(task);

  //
  // Create containers for input/output
  //
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(TString::Format("%s_ElectronEfficiency", directoryBaseName.Data()), TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(TString::Format("%s_supportHistos", directoryBaseName.Data()), TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(TString::Format("%s_EffTree", directoryBaseName.Data()), TTree::Class(),
                                                            AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(TString::Format("%s_stats", directoryBaseName.Data()), TH1D::Class(),
                                                            AliAnalysisManager::kOutputContainer,outputFileName);                                                          

  //connect input/output
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
  mgr->ConnectOutput(task,4,coutput4);

  return task;

}//AddTask
