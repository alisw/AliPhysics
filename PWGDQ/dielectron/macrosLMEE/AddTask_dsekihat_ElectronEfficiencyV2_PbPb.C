AliAnalysisTaskElectronEfficiencyV2* AddTask_dsekihat_ElectronEfficiencyV2_PbPb(
    Bool_t getFromAlien = kFALSE,
    TString configFile="Config_dsekihat_ElectronEfficiencyV2_PbPb.C",
    TString libFile="LMEECutLib_dsekihat.C",
    UInt_t trigger = AliVEvent::kINT7,
    const Int_t CenMin =  0,
    const Int_t CenMax = 10,
    const Bool_t applyPairCut = kTRUE,
    const TString pileupcut = "_woPU",//can be "_woPU", "_onlyPU",""
    const TString generators = "pizero_0;eta_1;etaprime_2;rho_3;omega_4;phi_5;jpsi_6;EvtGenDecay;Pythia CC_0;Pythia BB_0;Pythia B_0;",
    const TString calibFileName = "",
    const std::string resolutionFilename ="",
    const std::string cocktailFilename   ="",
    const std::string centralityFilename ="",
		const TString outname = "LMEE.root"
    ){

  // Configuring Analysis Manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

	TString suffix = "";
	if(generators.Contains("Pythia CC") && (generators.Contains("Pythia BB") || generators.Contains("Pythia B"))) suffix = "_CC_BB";
	else if(generators.Contains("Pythia CC")) suffix = "_CC";
	else if(generators.Contains("Pythia BB") || generators.Contains("Pythia B")) suffix = "_BB";
	else if(generators.Contains("pizero")) suffix = "_LF";
	else suffix = "";

  // Creating an instance of the task
  AliAnalysisTaskElectronEfficiencyV2* task = new AliAnalysisTaskElectronEfficiencyV2(Form("TaskElectronEfficiencyV2%s_Cen%d_%d_kINT7%s",suffix.Data(),CenMin,CenMax,pileupcut.Data()));
  gROOT->GetListOfSpecials()->Add(task);//this is only for ProcessLine(AddMCSignal);

  //TString configBasePath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/");
	TString configBasePath("./");
	if(getFromAlien
			&& (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/macrosLMEE/%s .",configFile.Data())))
			&& (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/macrosLMEE/%s .",libFile.Data())))
		){
		configBasePath=Form("%s/",gSystem->pwd());
	}
	TString configFilePath(configBasePath + configFile);

  //load cut library first
  TString libFilePath(configBasePath + libFile);
  std::cout << "Configpath:  " << configFilePath << std::endl;
  std::cout << "Libpath:  " << libFilePath << std::endl;

  gROOT->LoadMacro(libFilePath.Data());//library first
  gROOT->LoadMacro(configFilePath.Data());

  const Int_t nTC  = Int_t(gROOT->ProcessLine("GetNTC()") );
  const Int_t nPID = Int_t(gROOT->ProcessLine("GetNPID()"));
  const Int_t nPF  = Int_t(gROOT->ProcessLine("GetNPF()") );

	for (Int_t itc=0; itc<nTC; ++itc){
		for (Int_t ipid=0; ipid<nPID; ++ipid){
			for (Int_t ipf=0; ipf<nPF; ++ipf){
				AliAnalysisFilter *filter = reinterpret_cast<AliAnalysisFilter*>(gROOT->ProcessLine(Form("Config_dsekihat_ElectronEfficiencyV2_PbPb(%d,%d,%d,%d)",itc,ipid,ipf,applyPairCut)));
				task->AddTrackCuts(filter);
			}
		}
	}

  //It is important to apply pileup cuts at "task" level in M.C. for efficiency. (i.e. Both denominator and nominator has to be taken into account in exactly the same events.)

  // Event selection. Is the same for all the different cutsettings
  task->SetEnablePhysicsSelection(kTRUE);//always ON in Run2 analyses for both data and MC.
  task->SetTriggerMask(trigger);
  task->SetEventFilter(reinterpret_cast<AliDielectronEventCuts*>(gROOT->ProcessLine(Form("LMEECutLib::SetupEventCuts(%f,%f,%d,\"%s\")",(Float_t)CenMin,(Float_t)CenMax,kTRUE,"V0M"))));//kTRUE is for Run2
  if(pileupcut == ""){
    printf("analyze all events in M.C.\n");
  }
  else if(pileupcut.Contains("woPU")){
    printf("analyze only in clean events in M.C.\n");
    TF1 *f1pu = reinterpret_cast<TF1*>(gROOT->ProcessLine("LMEECutLib::SetupPileupCuts()"));
    dynamic_cast<AliDielectronEventCuts*>(task->GetEventFilter())->SetMinCorrCutFunction(f1pu, AliDielectronVarManager::kNTPCclsEvent, AliDielectronVarManager::kNSDDSSDclsEvent);
  }
  else if(pileupcut.Contains("onlyPU")){
    printf("analyze only in pileup events in M.C.\n");
    TF1 *f1pu = reinterpret_cast<TF1*>(gROOT->ProcessLine("LMEECutLib::SetupPileupCuts()"));
    dynamic_cast<AliDielectronEventCuts*>(task->GetEventFilter())->SetMaxCorrCutFunction(f1pu, AliDielectronVarManager::kNTPCclsEvent, AliDielectronVarManager::kNSDDSSDclsEvent);
  }
  else{
    printf("Nothing with pileup cut in M.C.\n");
    printf("analyze all events in M.C.\n");
  }
  dynamic_cast<AliDielectronEventCuts*>(task->GetEventFilter())->Print();

  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  // Do not set here your analysis pt-cuts
  task->SetMinPtGen(0.1);
  task->SetMaxPtGen(1e+10);
  task->SetMinEtaGen(-1.5);
  task->SetMaxEtaGen(+1.5);

  // Set minimum and maximum values for pairing
  task->SetKinematicCuts(0.2, 10, -0.8, +0.8);

  // Set Binning
  task->SetPtBinsLinear   (0, 10, 100);
  task->SetEtaBinsLinear  (-1, +1, 20);
  task->SetPhiBinsLinear  (0, TMath::TwoPi(), 36);
  task->SetThetaBinsLinear(0, TMath::TwoPi(), 60);

  const Int_t Nmee = 150;
  Double_t mee[Nmee] = {};
  for(Int_t i=0  ;i<110 ;i++) mee[i] = 0.01 * (i-  0) +  0.0;//from 0 to 1.1 GeV/c2, every 0.01 GeV/c2
  for(Int_t i=110;i<Nmee;i++) mee[i] = 0.1  * (i-110) +  1.1;//from 1.1 to 5 GeV/c2, evety 0.1 GeV/c2
	std::vector<double> v_mee(mee,std::end(mee));

  const Int_t NpTee = 146;
  Double_t pTee[NpTee] = {};
  for(Int_t i=0  ;i<50    ;i++) pTee[i] = 0.01 * (i-  0) +  0.0;//from 0 to 0.49 GeV/c, every 0.01 GeV/c
  for(Int_t i=50 ;i<NpTee ;i++) pTee[i] = 0.1  * (i- 50) +  0.5;//from 0.5 to 10 GeV/c, evety 0.1 GeV/c
  std::vector<double> v_pTee(pTee,std::end(pTee));

  task->SetMassBins(v_mee);
  task->SetPairPtBins(v_pTee);
  task->SetPhiVBinsLinear(0, TMath::Pi(), 100);
  task->SetFillPhiV(kFALSE);

  task->SetSmearGenerated(kFALSE);
  task->SetResolutionDeltaPtBinsLinear( -10,  +10, 2000);
  task->SetResolutionRelPtBinsLinear  (   0,    2,  400);
  task->SetResolutionEtaBinsLinear    (-0.4, +0.4,  200);
  task->SetResolutionPhiBinsLinear    (-0.4, +0.4,  200);
  task->SetResolutionThetaBinsLinear  (-0.4, +0.4,  200);

  // Set MCSignal and Cutsetting to fill the support histograms
  task->SetSupportHistoMCSignalAndCutsetting(0,0);//fill support histograms for first MCsignal and first cutsetting

  // Pairing related config
  task->SetDoPairing(kTRUE);
  task->SetULSandLS(kTRUE);
  task->SetDeactivateLS(kFALSE);

  //TString generators = "pizero_0;eta_1;etaprime_2;rho_3;omega_4;phi_5;jpsi_6;";
  //TString generators = "pizero_0;eta_1;etaprime_2;rho_3;omega_4;phi_5;jpsi_6;Pythia CC_0;Pythia BB_0;Pythia B_0;";
  //TString generators = "Pythia CC_0;Pythia BB_0;Pythia B_0;";
  //TString generators = "Pythia CC_0;";
  //TString generators = "Pythia BB_0;Pythia B_0;";
  //TString generators = "Pythia BB_0;";
  //TString generators = "Pythia B_0";

  cout<<"Efficiency based on MC generators: " << generators <<endl;
  TString generatorsPair=generators;
  task->SetGeneratorMCSignalName(generatorsPair);
  task->SetGeneratorULSSignalName(generators);
  //task->SetLHC19f2MC(isLHC19f2);

	//if(resolutionFilename != "") gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/resolution/%s .",resolutionFilename.c_str()));//this is to avoid unnecessary call of alien_cp in task.
	//if(cocktailFilename   != "") gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/cocktail/%s ."    ,cocktailFilename.c_str()));//this is to avoid unnecessary call of alien_cp in task.
	//if(centralityFilename != "") gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/centrality/%s .",centralityFilename.c_str()));//this is to avoid unnecessary call of alien_cp in task.
	//gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/resolution/%s .",resolutionFilename.c_str()));//this is to avoid unnecessary call of alien_cp in task.
	//gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/cocktail/%s ."    ,cocktailFilename.c_str()));//this is to avoid unnecessary call of alien_cp in task.
	//gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/centrality/%s .",centralityFilename.c_str()));//this is to avoid unnecessary call of alien_cp in task.

  // Resolution File, If resoFilename = "" no correction is applied
  task->SetResolutionFile(resolutionFilename);
  task->SetResolutionFileFromAlien("/alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/resolution/" + resolutionFilename);
  task->SetCocktailWeighting(cocktailFilename);
  task->SetCocktailWeightingFromAlien("/alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/cocktail/" + cocktailFilename);
  task->SetCentralityFile(centralityFilename);

  // Add MCSignals. Can be set to see differences of:
  // e.g. secondaries and primaries. or primaries from charm and resonances
  gROOT->ProcessLine(Form("AddSingleLegMCSignal(%s)",task->GetName()));//not task itself, task name
  gROOT->ProcessLine(Form("AddPairMCSignal(%s)"     ,task->GetName()));//not task itself, task name

  //set PID map for ITS TOF in MC.
  TFile *rootfile = 0x0;
  if(calibFileName != "") rootfile = TFile::Open(calibFileName,"READ");
  if(rootfile && rootfile->IsOpen()){
    TH3D *h3mean_ITS  = (TH3D*)rootfile->Get("h3mean_ITS");
    TH3D *h3width_ITS = (TH3D*)rootfile->Get("h3width_ITS");
    h3mean_ITS ->SetDirectory(0);
    h3width_ITS->SetDirectory(0);
    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kITS, h3mean_ITS, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta);
    task->SetWidthCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kITS, h3width_ITS, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta );

    TH3D *h3mean_TOF = (TH3D*)rootfile->Get("h3mean_TOF");
    TH3D *h3width_TOF = (TH3D*)rootfile->Get("h3width_TOF");
    h3mean_TOF ->SetDirectory(0);
    h3width_TOF->SetDirectory(0);
    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, h3mean_TOF, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta);
    task->SetWidthCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, h3width_TOF, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta);

    rootfile->Close();
  }
	TString outlistname = Form("Efficiency_dsekihat%s_Cen%d_%d_kINT7%s",suffix.Data(),CenMin,CenMax,pileupcut.Data());

  //const TString fileName = AliAnalysisManager::GetCommonFileName();
  const TString fileName = outname;
	//const TString dirname = Form("PWGDQ_LMEE_ElectronEfficiencyV2_Cen%d_%d_kINT7",CenMin,CenMax);
  mgr->AddTask(task);
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(outlistname, TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data())));
  return task;
}

