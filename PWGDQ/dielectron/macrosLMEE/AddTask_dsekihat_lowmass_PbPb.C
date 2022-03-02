AliAnalysisTask *AddTask_dsekihat_lowmass_PbPb(
		Bool_t getFromAlien = kFALSE,
		TString cFileName = "Config_dsekihat_lowmass_PbPb.C",
		TString lFileName = "LMEECutLib_dsekihat.C",
		TString calibFileName = "",
		UInt_t trigger = AliVEvent::kINT7|AliVEvent::kCentral|AliVEvent::kSemiCentral,
		const Int_t CenMin =  0,
		const Int_t CenMax = 10,
		const Int_t Nmix   = 10,
		const TString outname = "LMEE.root",
		const Bool_t isMC = kFALSE,
		const TString suffix=""
		)
{
	//get the current analysis manager
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTask_dsekihat_lowmass_PbPb", "No analysis manager found.");
		return 0;
	}

	//Base Directory for GRID / LEGO Train
	//TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
	TString configBasePath= "./";
  if (!gSystem->AccessPathName(cFileName) && !gSystem->AccessPathName(lFileName)  ) {
    printf("Configfile already present\n");
    configBasePath=Form("%s/",gSystem->pwd());
  }
	else if(getFromAlien
			&& (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/macrosLMEE/%s file:./",cFileName.Data())))
			&& (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/macrosLMEE/%s file:./",lFileName.Data())))
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

  //if (!gROOT->GetListOfClasses()->FindObject("LMEECutLib")) {
  //  printf("Load library now\n");
  //  gROOT->LoadMacro(libFilePath.Data());
  //  //gROOT->AddClass(LMEECutLib::Class());
  //}

  //if (!gROOT->GetListOfGlobalFunctions()->FindObject("Config_dsekihat_lowmass_PbPb")) {
  //  printf("Load macro now\n");
  //  gROOT->LoadMacro(configFilePath.Data());
  //}


	TString triggername = "NULL";
	if(trigger == (UInt_t)AliVEvent::kINT7)             triggername = "kINT7";
	else if(trigger == (UInt_t)AliVEvent::kCentral)     triggername = "kCentral";
	else if(trigger == (UInt_t)AliVEvent::kSemiCentral) triggername = "kSemiCentral";
	else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral)) triggername = "kCombinedCentralityTriggers";
	else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kCentral))                           triggername = "kCombinedCentral";
	else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kSemiCentral))                       triggername = "kCombinedSemiCentral";

	//create task and add it to the manager (MB)
	AliAnalysisTaskMultiDielectron *task = new AliAnalysisTaskMultiDielectron(Form("MultiDielectron_Cen%d_%d_%s%s",CenMin,CenMax,triggername.Data(),suffix.Data()));
	task->UsePhysicsSelection(kTRUE);
	task->SetTriggerMask(trigger);


	//const Int_t nDie = (Int_t)gROOT->ProcessLine("GetN()");
	const Int_t nEC = Int_t(gROOT->ProcessLine("GetNEC()") );//event cuts
	const Int_t nTC = Int_t(gROOT->ProcessLine("GetNTC()") );//track cuts
	const Int_t nPC = Int_t(gROOT->ProcessLine("GetNPID()"));//pid cuts

  printf("isMC : %d\n",isMC);
  TFile *rootfile = 0x0;
  //for data
  THnF *hs_mean_TPC_El  = 0x0;
  THnF *hs_width_TPC_El = 0x0;
  THnF *hs_mean_ITS_El  = 0x0;
  THnF *hs_width_ITS_El = 0x0;
  THnF *hs_mean_TOF_El  = 0x0;
  THnF *hs_width_TOF_El = 0x0;
  THnF *hs_mean_TPC_Pi  = 0x0;
  THnF *hs_width_TPC_Pi = 0x0;
  THnF *hs_mean_TPC_Ka  = 0x0;
  THnF *hs_width_TPC_Ka = 0x0;
  THnF *hs_mean_TPC_Pr  = 0x0;
  THnF *hs_width_TPC_Pr = 0x0;

  //for MC
  TH3D *h3mean_ITS  = 0x0;
  TH3D *h3width_ITS = 0x0;
  TH3D *h3mean_TOF  = 0x0;
  TH3D *h3width_TOF = 0x0;

  if(calibFileName != "") rootfile = TFile::Open(calibFileName,"READ");
  if(rootfile && rootfile->IsOpen()){
    printf("reading : %s for PID calibration\n",calibFileName.Data());
    if(!isMC){//for data
      printf("PID map in data\n");
      hs_mean_TPC_El  = (THnF*)rootfile->Get("hs_mean_TPC_El");
      hs_width_TPC_El = (THnF*)rootfile->Get("hs_width_TPC_El");
      hs_mean_ITS_El  = (THnF*)rootfile->Get("hs_mean_ITS_El");
      hs_width_ITS_El = (THnF*)rootfile->Get("hs_width_ITS_El");
      hs_mean_TOF_El  = (THnF*)rootfile->Get("hs_mean_TOF_El");
      hs_width_TOF_El = (THnF*)rootfile->Get("hs_width_TOF_El");
      hs_mean_TPC_Pi  = (THnF*)rootfile->Get("hs_mean_TPC_Pi");
      hs_width_TPC_Pi = (THnF*)rootfile->Get("hs_width_TPC_Pi");
      hs_mean_TPC_Ka  = (THnF*)rootfile->Get("hs_mean_TPC_Ka");
      hs_width_TPC_Ka = (THnF*)rootfile->Get("hs_width_TPC_Ka");
      hs_mean_TPC_Pr  = (THnF*)rootfile->Get("hs_mean_TPC_Pr");
      hs_width_TPC_Pr = (THnF*)rootfile->Get("hs_width_TPC_Pr");
    }
    else{//for MC
      printf("PID map in MC\n");
      h3mean_ITS  = (TH3D*)rootfile->Get("h3mean_ITS");
      h3width_ITS = (TH3D*)rootfile->Get("h3width_ITS");
      h3mean_TOF  = (TH3D*)rootfile->Get("h3mean_TOF");
      h3width_TOF = (TH3D*)rootfile->Get("h3width_TOF");
    }
  }
  for (Int_t iec=0; iec<nEC; ++iec){
    for (Int_t itc=0; itc<nTC; ++itc){
      for (Int_t ipc=0; ipc<nPC; ++ipc){
        AliDielectron *die = reinterpret_cast<AliDielectron*>(gROOT->ProcessLine(Form("Config_dsekihat_lowmass_PbPb(%d,%d,%d,%d)",iec,itc,ipc,isMC)));
        if(!die) continue;
        TString name = die->GetName();

        if(rootfile && rootfile->IsOpen()){
          if(!isMC){//for data
            die->SetPIDCaibinPU(kTRUE);
            //for electron
            die->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kElectron,hs_mean_TPC_El ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
            die->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kElectron,hs_width_TPC_El,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
            die->SetCentroidCorrFunctionPU(AliDielectronPID::kITS,AliPID::kElectron,hs_mean_ITS_El ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
            die->   SetWidthCorrFunctionPU(AliDielectronPID::kITS,AliPID::kElectron,hs_width_ITS_El,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
            die->SetCentroidCorrFunctionPU(AliDielectronPID::kTOF,AliPID::kElectron,hs_mean_TOF_El ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
            die->   SetWidthCorrFunctionPU(AliDielectronPID::kTOF,AliPID::kElectron,hs_width_TOF_El,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
            //for pion
            die->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kPion,hs_mean_TPC_Pi     ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
            die->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kPion,hs_width_TPC_Pi    ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
            //for kaon
            die->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kKaon,hs_mean_TPC_Ka     ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
            die->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kKaon,hs_width_TPC_Ka    ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
            //for proton
            die->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kProton,hs_mean_TPC_Pr   ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
            die->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kProton,hs_width_TPC_Pr  ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);

          }
          else{//for MC, TPC is perfectly calibrated.
            h3mean_ITS->SetDirectory(0);
            h3width_ITS->SetDirectory(0);
            die->SetCentroidCorrFunctionITS(h3mean_ITS, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta);//for ITS
            die->SetWidthCorrFunctionITS(h3width_ITS, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta);//for ITS

            h3mean_TOF->SetDirectory(0);
            h3width_TOF->SetDirectory(0);
            die->SetCentroidCorrFunctionTOF(h3mean_TOF, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta);//for TOF
            die->SetWidthCorrFunctionTOF(h3width_TOF, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta);//for TOF
          }
        }

        //if(name.Contains("PIDCalib",TString::kIgnoreCase) || name.Contains("noPID",TString::kIgnoreCase)){
        //  printf("No event mixing handler for post PID calibration/noPID\n");
        //}
        //else{
        //  printf("Add event mixing handler\n");
        //  AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;
        //  mix->SetMixType(AliDielectronMixingHandler::kAll);
        //  mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -8., -6., -4., -2. , 0., 2., 4., 6., 8. , 10.");
        //  mix->AddVariable(AliDielectronVarManager::kCentralityNew,"0,5,10,30,50,70,90,101");
        //  mix->AddVariable(AliDielectronVarManager::kNacc,"0,10000");
        //  mix->SetDepth(Nmix);
        //  if(!isMC) die->SetMixingHandler(mix);
        //}

        //mixing handler is added in Config.C
        AliDielectronMixingHandler *mix = dynamic_cast<AliDielectronMixingHandler*>(die->GetMixingHandler());
        if(mix) mix->SetDepth(Nmix);
        else printf("No event mixing handler is found\n");

        task->AddDielectron(die);
      }//pid loop
    }//track loop
  }//event loop

		//if(rootfile && rootfile->IsOpen()) rootfile->Close();

		task->SetEventFilter(reinterpret_cast<AliDielectronEventCuts*>(gROOT->ProcessLine(Form("LMEECutLib::SetupEventCuts(%f,%f,%d,\"%s\")",(Float_t)CenMin,(Float_t)CenMax,kTRUE,"V0M"))));//kTRUE is for Run2

		mgr->AddTask(task);

		//const TString outputFileName = AliAnalysisManager::GetCommonFileName();
		const TString outputFileName = outname;
		//const TString dirname = Form("PWGDQ_LMEE_Cen%d_%d_%s",CenMin,CenMax,triggername.Data());

		//create output container
		AliAnalysisDataContainer *coutput1 =
			mgr->CreateContainer(Form("Tree_dsekihat_Cen%d_%d_%s%s",CenMin,CenMax,triggername.Data(),suffix.Data()),
					TTree::Class(),
					AliAnalysisManager::kExchangeContainer,
					Form("%s",outputFileName.Data()));

		AliAnalysisDataContainer *cOutputHist1 =
			mgr->CreateContainer(Form("Histos_dsekihat_Cen%d_%d_%s%s",CenMin,CenMax,triggername.Data(),suffix.Data()),
					TList::Class(),
					AliAnalysisManager::kOutputContainer,
					Form("%s",outputFileName.Data()));

		AliAnalysisDataContainer *cOutputHist2 =
			mgr->CreateContainer(Form("CF_dsekihat_Cen%d_%d_%s%s",CenMin,CenMax,triggername.Data(),suffix.Data()),
					TList::Class(),
					AliAnalysisManager::kOutputContainer,
					Form("%s",outputFileName.Data()));

		AliAnalysisDataContainer *cOutputHist3 =
			mgr->CreateContainer(Form("EventStat_dsekihat_Cen%d_%d_%s%s",CenMin,CenMax,triggername.Data(),suffix.Data()),
					TH1D::Class(),
					AliAnalysisManager::kOutputContainer,
					Form("%s",outputFileName.Data()));

		mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
		mgr->ConnectOutput(task, 0, coutput1 );
		mgr->ConnectOutput(task, 1, cOutputHist1);
		mgr->ConnectOutput(task, 2, cOutputHist2);
		mgr->ConnectOutput(task, 3, cOutputHist3);

		return task;
	}

