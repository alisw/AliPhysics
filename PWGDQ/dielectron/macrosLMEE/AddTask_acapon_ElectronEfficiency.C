AliAnalysisTask *AddTask_acapon_ElectronEfficiency_v2(TString configFile="Config_acapon_ElectronEfficiency_v2.C",
                                                     Double_t centMin=0, Double_t centMax=90,
                                                     Bool_t getFromAlien=kFALSE,
                                                     Bool_t deactivateTree=kFALSE, // enabling this has priority over 'writeTree' in config file! (enable for LEGO trains)
                                                     Char_t* outputFileName="acapon_ElectronEfficiency.root",
                                                     Bool_t forcePhysSelAndTrigMask=kTRUE, // possibility to activate UsePhysicsSelection and SetTriggerMask for MC (may be needed for new MC productions according to Mahmut) as well as for AOD data.
                                                     Int_t triggerNames=(AliVEvent::kMB),
                                                     Int_t collCands=AliVEvent::kMB,
                                                     Int_t cutlibPreloaded=0
                                                     )
{
  // get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_acapon_ElectronEfficiency_v2", "No analysis manager found.");
    return 0;
  }

	// Base Directory for GRID / LEGO Train  
	TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
	if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/a/acapon/PWGDQ/dielectron/macrosLMEE/%s .",cFileName.Data()))) ){
		
		configBasePath=Form("%s/",gSystem->pwd());
	}

	TString configFilePath(configBasePath+cFileName);
	std::cout << "Configpath:  " << configFilePath << std::endl;

	if(gSystem->Exec(Form("ls %s", configFilePath.Data()))==0){

		std::cout << "loading config: " << configFilePath.Data() << std::endl;
		gROOT->LoadMacro(configFilePath.Data());
	} 
	else{
		std::cout << "config not found: " << configFilePath.Data() << std::endl;
		return 0; // if return is not called, the job will fail instead of running wihout this task... (good for local tests, bad for train)
	}

  std::cout << "computing binning..." << std::endl;
  Double_t EtaBins[nBinsEta+1];
  for(Int_t i = 0; i <= nBinsEta; i++){
		EtaBins[i] = EtaMin + i*(EtaMax-EtaMin)/nBinsEta;
	}
  Double_t PhiBins[nBinsPhi+1];
  for(Int_t i = 0; i <= nBinsPhi; i++){
		PhiBins[i] = PhiMin + i*(PhiMax-PhiMin)/nBinsPhi;
	}

  const Int_t nBinsPt =  ( sizeof(PtBins) / sizeof(PtBins[0]) )-1;

  Bool_t bESDANA=kFALSE; //Autodetect via InputHandler
  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
    ::Info("AddTaskElectronEfficiency","running on AODs.");
  }
  else if (mgr->GetInputEventHandler()->IsA()==AliESDInputHandler::Class()){
    ::Info("AddTaskElectronEfficiency","switching on ESD specific code, make sure ESD cuts are used.");
    bESDANA=kTRUE;
  }
  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);
  std::cout << "hasMC = " << hasMC << std::endl;

  // Electron efficiency task
  AliAnalysisTaskElectronEfficiency *task = new AliAnalysisTaskElectronEfficiency("acapon_ElectronEfficiency_v2");
  std::cout << "task created: " << task->GetName() << std::endl;
  task->SetUseMultSelection(kTRUE);

  // Event related
  // Note: event cuts are identical for all analysis 'cutInstance's that run together!
  if(!hasMC) task->UsePhysicsSelection();
  if(!hasMC) task->SetTriggerMask(triggerNames);
	// Next should be kFALSE by default! (may be needed for new MC productions according to Mahmut)
  if(forcePhysSelAndTrigMask){ 
    task->UsePhysicsSelection();
    task->SetTriggerMask(triggerNames);
  }
  task->SelectCollisionCandidates(collCands); // didnt check its meaning and effect...!
  //---
  task->SetEventFilter(SetupEventCuts()); //returns eventCuts from Config.
  task->SetCentralityRange(centMin, centMax);
  task->SetNminEleInEventForRej(NminEleInEventForRej);
  //track related
  task->SetEtaRangeGEN(EtaMinGEN, EtaMaxGEN);
  task->SetPtRangeGEN(PtMinGEN, PtMaxGEN);
  task->SetCalcEfficiencyPoslabel(CalcEfficiencyPoslabel);

  // Resolution calculation
  task->SetCalcResolution(CalcResolution);
  if(CalcResolution || CalcEfficiencyRec){
		task->SetResolutionCuts(SetupTrackCutsAndSettings("TTree_Cuts"));
	}
  if(CalcResolution){
    task->SetDeltaMomBinning(NbinsDeltaMom,DeltaMomMin,DeltaMomMax);
    task->SetRelMomBinning(NbinsRelMom,RelMomMin,RelMomMax);
    task->SetDeltaEtaBinning(NbinsDeltaEta,DeltaEtaMin,DeltaEtaMax);
    task->SetDeltaThetaBinning(NbinsDeltaTheta,DeltaThetaMin,DeltaThetaMax);
    task->SetDeltaPhiBinning(NbinsDeltaPhi,DeltaPhiMin,DeltaPhiMax);
  }
  // resolution usage
  task->SetCalcEfficiencyRec(CalcEfficiencyRec);
  if(resolutionfile.Contains("CENTRALITY")){
		resolutionfile.ReplaceAll("CENTRALITY",Form("%02.0f%02.0f",centMin,centMax));
	}
  if(CalcEfficiencyRec && !resolutionfile.IsNull()){
		if( (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/a/acapon/ResolutionFiles/%s .",resolutionfile.Data()))) ){
			TFile *fRes = TFile::Open(Form("%s/%s",gSystem->pwd(),resolutionfile.Data()),"READ");
			if(!fRes){
				Error("---- Resolution file was specified but not opened. Exiting. ------");
				return 0;
			}
			if(!((TObjArray*)fRes->Get("DeltaPtResArr"))){
				std::cout << "Pt res file not found" << std::endl;
			}
	  	task->SetResolutionP ((TObjArray*) fRes->Get("DeltaPtResArr"),kFALSE);
	  	if(bUseEtaResolution){
				if(!((TObjArray*)fRes->Get("EtaResArr"))){
					std::cout << "Eta res file not found" << std::endl;
				}
				task->SetResolutionEta ((TObjArray*) fRes->Get("EtaResArr"));
	  	}else{
				task->SetResolutionTheta((TObjArray*) fRes->Get("ThetaResArr"));
			}
			if(!((TObjArray*)fRes->Get("PhiEleResArr")) || !((TObjArray*)fRes->Get("PhiPosResArr"))){
				std::cout << "Phi res file not found" << std::endl;
			}
	  	task->SetResolutionPhi( (TObjArray*) fRes->Get("PhiEleResArr"), (TObjArray*) fRes->Get("PhiPosResArr"));
		}
		else{
			std::cout << "Resolution file was not copied. Exiting." << std::endl;
			return 0;
		}
	}

  // pair efficiency
  task->SetDoPairing(doPairing);
  if(doPairing){
    task->SetPtCut(PtMinCut, PtMaxCut);
    task->SetEtaCut(EtaMinCut, EtaMaxCut);
    task->SetKineTrackCuts(SetupTrackCutsAndSettings(0));
    //task->SetPairCuts(SetupTrackCutsAndSettings(101));
    // SetupTrackCutsAndSettings(101); // this fills the pair cuts into rejCutMee,rejCutTheta,rejCutPhiV
    task->SetPairCutMee(rejCutMee);
    task->SetPairCutTheta(rejCutTheta);
    task->SetPairCutPhiV(rejCutPhiV);
    Double_t MeeBins[nBinsMee+1];
    for(Int_t i=0;i<=nBinsMee;i++){
			MeeBins[i] = MeeMin + i*(MeeMax-MeeMin)/nBinsMee;
		}
    Double_t PteeBins[nBinsPtee+1];
    for(Int_t i=0;i<=nBinsPtee;i++){
			PteeBins[i] = PteeMin + i*(PteeMax-PteeMin)/nBinsPtee;
		}
    task->SetBins(nBinsPt,PtBins,nBinsEta,EtaBins,nBinsPhi,PhiBins,nBinsMee,MeeBins,nBinsPtee,PteeBins);
  }else{
    task->SetBins(nBinsPt,PtBins,nBinsEta,EtaBins,nBinsPhi,PhiBins);
  }
  //output related
  task->SetRunBins(sRuns);
  if(deactivateTree){
		task->SetWriteTree(kFALSE);
	}else{
		task->SetWriteTree(writeTree);
	}
  task->SetSupportedCutInstance(supportedCutInstance);
  task->CreateHistoGen();

  // Post PID correction
  /* SetupITSSigmaEleCorrection(task); */
  /* SetupTPCSigmaEleCorrection(task); */

  // Monte Carlo Signals
  SetupMCSignals(task);

  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    TString cutDefinition(arrNames->At(i)->GetName());
    std::cout << "##############################" << std::endl;
    std::cout << "NEW CUT SETTING WITH CUT NAME: " << cutDefinition << std::endl;
    Bool_t IsTOFreqCut = false;
    AliAnalysisFilter *trackCuts = SetupTrackCutsAndSettings(cutDefinition); // main function in config file
    if(!trackCuts){
			std::cout << "WARNING: no TrackCuts given - skipping this Cutset ('"<<arrNames->At(i)->GetName()<<"')!" << std::endl;
			continue;
		}
    
    // fill std vectors with all information which is individual per track setting:
    task->AttachTrackCuts(trackCuts);
    task->AttachRejCutMee(rejCutMee);
    task->AttachRejCutTheta(rejCutTheta);
    task->AttachRejCutPhiV(rejCutPhiV);

    task->CreateHistograms(names,i);
    std::cout << "Cutsetting added\n##############################" << std::endl;
  }

  mgr->AddTask(task);

  //
  // Create containers for input/output
  //
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("acapon_ElectronEfficiency_v2_%d", cutlibPreloaded), TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("acapon_supportHistos_%d", cutlibPreloaded), TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(Form("acapon_EffTree_%d", cutlibPreloaded), TTree::Class(),
                                                            AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(Form("acapon_stats_%d", cutlibPreloaded), TH1D::Class(),
                                                            AliAnalysisManager::kOutputContainer,outputFileName);

  //connect input/output
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
  mgr->ConnectOutput(task,4,coutput4);

  printf("ADDTASK RUN SUCCESSFULL\n");

  return task;

}//AddTask
