AliAnalysisTaskSimpleTreeMaker *AddTaskSimpleTreeMaker(TString taskName    = "MLtree",
                                                       Bool_t hasSDD       = kTRUE,
                                                       Bool_t useITScorr   = kFALSE,
                                                       Bool_t useTPCcorr   = kFALSE,
                                                       Bool_t useTOFcorr   = kFALSE,
                                                       Bool_t isMC         = kFALSE,
                                                       Bool_t runOnGrid    = kTRUE,
                                                       Bool_t getFromAlien = kFALSE,
						                                           Bool_t ExtraDCA = kFALSE)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskSimpleTreeMaker",  "No analysis manager to connect to.");
        return NULL;
    }

    TString configBasePath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/"); //AliPhysics
    TString configLMEECutLib("LMEECutLib_acapon.C");

    //Load updated macros from private ALIEN path
    TString myCutLib ="alien_cp alien:///alice/cern.ch/user/a/acapon/PWGDQ/dielectron/macrosLMEE/LMEECutLib_acapon.C ."; 
    if(getFromAlien && (!gSystem->Exec(myCutLib))){

        std::cout << "Copy config from Alien" << std::endl;
        configBasePath=Form("%s/",gSystem->pwd());
    }

    TString configLMEECutLibPath(configBasePath+configLMEECutLib);

    //load dielectron configuration files
    if(!gROOT->GetListOfGlobalFunctions()->FindObject(configLMEECutLib.Data())){
        gROOT->LoadMacro(configLMEECutLibPath.Data());
    }

    // Check the analysis type using the event handlers connected to the analysis manager.
    //===========================================================================
    if(!mgr->GetInputEventHandler()){
			::Error("AddTaskSimpleTreeMaker",  "This task requires an input event handler");
			return NULL;
    }

    TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

    if(analysisType != "ESD" && analysisType != "AOD"){
			::Error("AddTaskSimpleTreeMaker",  "analysis type NOT AOD or ESD --> makes no sense!");
			return NULL;
    }

		LMEECutLib* cutLib = new LMEECutLib(hasSDD);
		AliAnalysisTaskSimpleTreeMaker *task = new AliAnalysisTaskSimpleTreeMaker(taskName, ExtraDCA);
		task->analyseMC(isMC);
		task->GRIDanalysis(runOnGrid);
    // ==========================================================================
    // user customization part

		// TPC in MC is already calibrated
    if(useTPCcorr && !isMC){
			TH3D meanTPC = cutLib->SetEtaCorrectionTPCTTree(AliDielectronVarManager::kP,
                                              AliDielectronVarManager::kEta,
                                              AliDielectronVarManager::kRefMultTPConly, 1);

			TH3D widthTPC = cutLib->SetEtaCorrectionTPCTTree(AliDielectronVarManager::kP,
                                               AliDielectronVarManager::kEta,
                                               AliDielectronVarManager::kRefMultTPConly, 2);
			task->SetUseTPCcorr(kTRUE);
			task->SetCorrWidthMeanTPC((TH3D*)widthTPC.Clone(),(TH3D*)meanTPC.Clone());

		}

    if(useITScorr){
			TH3D meanITS = cutLib->SetEtaCorrectionITSTTree(AliDielectronVarManager::kP,
                                              AliDielectronVarManager::kEta,
                                              AliDielectronVarManager::kRefMultTPConly, 1, isMC);

			TH3D widthITS = cutLib->SetEtaCorrectionITSTTree(AliDielectronVarManager::kP,
                                               AliDielectronVarManager::kEta,
                                               AliDielectronVarManager::kRefMultTPConly, 2, isMC);
			task->SetUseITScorr(kTRUE);
			task->SetCorrWidthMeanITS((TH3D*)widthITS.Clone(),(TH3D*)meanITS.Clone());

		}
		if(useTOFcorr){
			TH3D meanTOF = cutLib->SetEtaCorrectionTOFTTree(AliDielectronVarManager::kP,
                                              AliDielectronVarManager::kEta,
                                              AliDielectronVarManager::kRefMultTPConly, 1, isMC);

			TH3D widthTOF = cutLib->SetEtaCorrectionTOFTTree(AliDielectronVarManager::kP,
                                               AliDielectronVarManager::kEta,
                                               AliDielectronVarManager::kRefMultTPConly, 2, isMC);
			task->SetUseTOFcorr(kTRUE);
			task->SetCorrWidthMeanTOF((TH3D*)widthTOF.Clone(),(TH3D*)meanTOF.Clone());

		}
    //Add event filter
		task->SelectCollisionCandidates(AliVEvent::kINT7);
    task->SetupEventCuts(cutLib->GetEventCuts(kFALSE, kFALSE));
		task->SetupTrackCuts(cutLib->GetTrackCuts(LMEECutLib::kTTreeCuts, LMEECutLib::kTTreeCuts));

    // ==========================================================================
    mgr->AddTask(task);

    // Create ONLY the output containers for the data produced by the task.
    // Get and connect other common input/output containers via the manager as below
    //===================================================mcHandler->========

    TString outputFileName = AliAnalysisManager::GetCommonFileName();  

    AliAnalysisDataContainer *coutTree = mgr->CreateContainer("Tree", TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
    AliAnalysisDataContainer *coutHisto1 = mgr->CreateContainer("Histo", TH1F::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
    AliAnalysisDataContainer *coutHisto2 = mgr->CreateContainer("Arm. Plot", TH2F::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutTree);
    mgr->ConnectOutput(task, 2, coutHisto1);
    mgr->ConnectOutput(task, 3, coutHisto2);

    return task;
}
