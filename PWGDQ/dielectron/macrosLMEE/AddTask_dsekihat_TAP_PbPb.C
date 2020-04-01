AliAnalysisTaskTagAndProbe* AddTask_dsekihat_TAP_PbPb(
    Bool_t getFromAlien = kFALSE,
    TString cFileName = "Config_dsekihat_TAP_PbPb.C",
    TString lFileName = "LMEECutLib_dsekihat.C",
		TString calibFileName = "",
    UInt_t trigger = AliVEvent::kINT7|AliVEvent::kCentral|AliVEvent::kSemiCentral,
    const Int_t CenMin =  0,
    const Int_t CenMax = 10,
    const Int_t Nmix   = 10,
		const TString type = "PID",//what kind of efficiency do you want to measure by TAP?//PID, Tracking, etc
		const TString cutname = "DefaultTrackCut_Ns01_TPChadrejORTOFrec",//same as real anslysis used in AliAnalysisTaskMultDielectron
		const TString outname = "LMEE.root"
    )
{
  //Add a task AliAnalysisTaskTagAndProbe to the analysis train
  //Author: Daiki Sekihata

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask_dsekihat_PbPb", "No analysis manager to connect to");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask_dsekihat_PbPb", "This task requires an input event handler");
    return NULL;
  }

	//load config files and cut library

  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  //TString configBasePath= "./";
  if(getFromAlien
    && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/macrosLMEE/%s .",cFileName.Data())))
    && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/macrosLMEE/%s .",lFileName.Data())))
    ){
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath + cFileName);
  TString libFilePath(configBasePath + lFileName);
  std::cout << "Configpath:  " << configFilePath << std::endl;
  std::cout << "Libpath:  "    << libFilePath << std::endl;

  //add dielectron analysis with different cuts to the task
  gROOT->LoadMacro(libFilePath.Data());//library first
  gROOT->LoadMacro(configFilePath.Data());

	TString triggername = "NULL";
	if(trigger == (UInt_t)AliVEvent::kINT7)             triggername = "kINT7";
	else if(trigger == (UInt_t)AliVEvent::kCentral)     triggername = "kCentral";
	else if(trigger == (UInt_t)AliVEvent::kSemiCentral) triggername = "kSemiCentral";
	else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral)) triggername = "kCombinedCentralityTriggers";
	else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kCentral))                           triggername = "kCombinedCentral";
	else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kSemiCentral))                       triggername = "kCombinedSemiCentral";

  TString taskname = Form("TAP_%sEff_%s_dsekihat_Cen%d_%d_%s",type.Data(),cutname.Data(),CenMin,CenMax,triggername.Data());
  AliAnalysisTaskTagAndProbe* task = new AliAnalysisTaskTagAndProbe(taskname);
  gROOT->GetListOfSpecials()->Add(task);//this is only for ProcessLine(Config_dsekihat(taskname));

  task->SelectCollisionCandidates(trigger);
	task->SetTriggerMask(trigger);
  task->SetCentralityMin(CenMin);
  task->SetCentralityMax(CenMax);
  task->SetDepthNMixed(Nmix);
  //centrality setting
  task->SetCentralityEstimator("V0M");

	//task->SetPhivCutRange(0.05,2.0);

  task->GetEventFilter()->AddCuts(reinterpret_cast<AliDielectronEventCuts*>(gROOT->ProcessLine(Form("LMEECutLib::SetupEventCuts(%f,%f,%d,\"%s\")",(Float_t)CenMin,(Float_t)CenMax,kTRUE,"V0M"))));//kTRUE is for Run2

	gROOT->ProcessLine(Form("Config_dsekihat_TAP_PbPb(%s,\"%s\",\"%s\")",task->GetName(),type.Data(),cutname.Data()));

	//PID calibration
	THnF *hs_mean_TPC_El  = 0x0;
	THnF *hs_width_TPC_El = 0x0;
	THnF *hs_mean_ITS_El  = 0x0;
	THnF *hs_width_ITS_El = 0x0;
	THnF *hs_mean_TOF_El  = 0x0;
	THnF *hs_width_TOF_El = 0x0;
	THnF *hs_mean_TPC_Pi  = 0x0;
	THnF *hs_width_TPC_Pi = 0x0;

	TFile *rootfile = 0x0;
	if(calibFileName != ""){
		printf("reading : %s for PID calibration\n",calibFileName.Data());
		rootfile = TFile::Open(calibFileName,"READ");
		hs_mean_TPC_El  = (THnF*)rootfile->Get("hs_mean_TPC_El");
		hs_width_TPC_El = (THnF*)rootfile->Get("hs_width_TPC_El");
		hs_mean_ITS_El  = (THnF*)rootfile->Get("hs_mean_ITS_El");
		hs_width_ITS_El = (THnF*)rootfile->Get("hs_width_ITS_El");
		hs_mean_TOF_El  = (THnF*)rootfile->Get("hs_mean_TOF_El");
		hs_width_TOF_El = (THnF*)rootfile->Get("hs_width_TOF_El");
		hs_mean_TPC_Pi  = (THnF*)rootfile->Get("hs_mean_TPC_Pi");
		hs_width_TPC_Pi = (THnF*)rootfile->Get("hs_width_TPC_Pi");
	}

	if(rootfile && rootfile->IsOpen()){
		task->SetPIDCaibinPU(kTRUE);
		//for electron
		task->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kElectron,hs_mean_TPC_El ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
		task->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kElectron,hs_width_TPC_El,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
		task->SetCentroidCorrFunctionPU(AliDielectronPID::kITS,AliPID::kElectron,hs_mean_ITS_El ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
		task->   SetWidthCorrFunctionPU(AliDielectronPID::kITS,AliPID::kElectron,hs_width_ITS_El,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
		task->SetCentroidCorrFunctionPU(AliDielectronPID::kTOF,AliPID::kElectron,hs_mean_TOF_El ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
		task->   SetWidthCorrFunctionPU(AliDielectronPID::kTOF,AliPID::kElectron,hs_width_TOF_El,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
		//for pion
		task->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kPion,hs_mean_TPC_Pi ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
		task->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kPion,hs_width_TPC_Pi,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
	}


  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );

  //TString outputFile = AliAnalysisManager::GetCommonFileName();
  TString outputFile = outname;
  TString prefix = Form("hist_%s",taskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("%s",prefix.Data()), THashList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",outputFile.Data()));
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}

