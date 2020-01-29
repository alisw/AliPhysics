AliAnalysisTask *AddTask_dsekihat_lowmass_PbPb(
    Bool_t getFromAlien = kFALSE,
    TString cFileName = "Config_dsekihat_lowmass_PbPb.C",
    TString lFileName = "LMEECutLib_dsekihat.C",
		TString calibFileName = "",
    UInt_t trigger = AliVEvent::kINT7|AliVEvent::kCentral|AliVEvent::kSemiCentral,
    const Int_t CenMin =  0,
    const Int_t CenMax = 10,
    const Int_t Nmix   = 10,
    const Float_t PtMin = 0.2,
    const Float_t PtMax = 10,
    const Float_t EtaMin = -0.8,
    const Float_t EtaMax = +0.8,
		const TString outname = "LMEE.root"
    )
{
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_dsekihat_lowmass_PbPb", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train
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

  TString triggername = "NULL";
  if(trigger == (UInt_t)AliVEvent::kINT7)             triggername = "kINT7";
  else if(trigger == (UInt_t)AliVEvent::kCentral)     triggername = "kCentral";
  else if(trigger == (UInt_t)AliVEvent::kSemiCentral) triggername = "kSemiCentral";
  else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral)) triggername = "kCombinedCentralityTriggers";
  else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kCentral))                           triggername = "kCombinedCentral";
  else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kSemiCentral))                       triggername = "kCombinedSemiCentral";

  //create task and add it to the manager (MB)
  AliAnalysisTaskMultiDielectron *task = new AliAnalysisTaskMultiDielectron(Form("MultiDielectron_Cen%d_%d_%s",CenMin,CenMax,triggername.Data()));
  task->UsePhysicsSelection(kTRUE);
  task->SetTriggerMask(trigger);

  //add dielectron analysis with different cuts to the task
  gROOT->LoadMacro(libFilePath.Data());//library first
  gROOT->LoadMacro(configFilePath.Data());

  //const Int_t nDie = (Int_t)gROOT->ProcessLine("GetN()");
  const Int_t nTC  = Int_t(gROOT->ProcessLine("GetNTC()") );
  const Int_t nPID = Int_t(gROOT->ProcessLine("GetNPID()"));
  const Int_t nPF  = Int_t(gROOT->ProcessLine("GetNPF()") );

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

  for (Int_t itc=0; itc<nTC; ++itc){
    for (Int_t ipid=0; ipid<nPID; ++ipid){
      for (Int_t ipf=0; ipf<nPF; ++ipf){
        AliDielectron *diel = reinterpret_cast<AliDielectron*>(gROOT->ProcessLine(Form("Config_dsekihat_lowmass_PbPb(%d,%d,%d,%f,%f,%f,%f)",itc,ipid,ipf,PtMin,PtMax,EtaMin,EtaMax)));
        if(!diel) continue;

        //AliDielectronVarCuts*  centCuts = new AliDielectronVarCuts("centCuts",Form("kPbPb%02d%02d",CenMin,CenMax));
        //centCuts->AddCut(AliDielectronVarManager::kCentralityNew, (Float_t)CenMin, (Float_t)CenMax);
        //diel->GetEventFilter().AddCuts(centCuts);

				//if(calibFileName!=""){
				if(rootfile && rootfile->IsOpen()){
					diel->SetPIDCaibinPU(kTRUE);
					//for electron
					diel->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kElectron,hs_mean_TPC_El ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
					diel->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kElectron,hs_width_TPC_El,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
					diel->SetCentroidCorrFunctionPU(AliDielectronPID::kITS,AliPID::kElectron,hs_mean_ITS_El ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
					diel->   SetWidthCorrFunctionPU(AliDielectronPID::kITS,AliPID::kElectron,hs_width_ITS_El,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
					diel->SetCentroidCorrFunctionPU(AliDielectronPID::kTOF,AliPID::kElectron,hs_mean_TOF_El ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
					diel->   SetWidthCorrFunctionPU(AliDielectronPID::kTOF,AliPID::kElectron,hs_width_TOF_El,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
					//for pion
					diel->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kPion,hs_mean_TPC_Pi ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
					diel->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kPion,hs_width_TPC_Pi ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta,AliDielectronVarManager::kNclsITS1,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM);
				}

        TString name = diel->GetName();
        if(name.Contains("PIDCalib",TString::kIgnoreCase) || name.Contains("noPID",TString::kIgnoreCase)){
          printf("No event mixing handler for post PID calibration/noPID\n");
        }
        else{
          printf("Add event mixing handler\n");
          AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;
          mix->SetMixType(AliDielectronMixingHandler::kAll);
          mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -8., -6., -4., -2. , 0., 2., 4., 6., 8. , 10.");
          mix->AddVariable(AliDielectronVarManager::kCentralityNew,"0,5,10,30,50,70,90,101");
          mix->AddVariable(AliDielectronVarManager::kNacc,"0,10000");
          mix->SetDepth(Nmix);
          diel->SetMixingHandler(mix);
        }
        task->AddDielectron(diel);
      }//pre-filter loop
    }//PID loop
  }//track cut loop

	//if(rootfile && rootfile->IsOpen()) rootfile->Close();

  task->SetEventFilter(reinterpret_cast<AliDielectronEventCuts*>(gROOT->ProcessLine(Form("LMEECutLib::SetupEventCuts(%f,%f,%d,\"%s\")",(Float_t)CenMin,(Float_t)CenMax,kTRUE,"V0M"))));//kTRUE is for Run2

  mgr->AddTask(task);

  //const TString outputFileName = AliAnalysisManager::GetCommonFileName();
  const TString outputFileName = outname;
  //const TString dirname = Form("PWGDQ_LMEE_Cen%d_%d_%s",CenMin,CenMax,triggername.Data());

  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(Form("Tree_dsekihat_Cen%d_%d_%s",CenMin,CenMax,triggername.Data()),
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         Form("%s",outputFileName.Data()));

  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer(Form("Histos_dsekihat_Cen%d_%d_%s",CenMin,CenMax,triggername.Data()),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         Form("%s",outputFileName.Data()));

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer(Form("CF_dsekihat_Cen%d_%d_%s",CenMin,CenMax,triggername.Data()),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         Form("%s",outputFileName.Data()));

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer(Form("EventStat_dsekihat_Cen%d_%d_%s",CenMin,CenMax,triggername.Data()),
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

