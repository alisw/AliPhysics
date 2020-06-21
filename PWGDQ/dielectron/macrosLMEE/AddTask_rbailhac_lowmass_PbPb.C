AliAnalysisTask *AddTask_rbailhac_lowmass_PbPb(Bool_t getFromAlien=kFALSE,
					       TString cFileName = "Config_rbailhac_lowmass_PbPb.C",
					       TString calibFileName = "",
					       UInt_t trigger = AliVEvent::kINT7|AliVEvent::kCentral|AliVEvent::kSemiCentral,
					       const Int_t CenMin =  0,
					       const Int_t CenMax = 10,
					       const Bool_t isMC = kFALSE,
					       const Bool_t isMix = kTRUE,
					       const Int_t Nmix   = 10,
					       Char_t* outputFileName="LMEE.root",
					       Bool_t rejpileup = kTRUE
					       )
{

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_rbailhac_test", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/r/rbailhac/PWGDQ/dielectron/macrosLMEE/%s .",cFileName.Data()))) ){
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);

  std::cout << "Configpath:  " << configFilePath << std::endl;

  if (!gROOT->GetListOfGlobalFunctions()->FindObject("Config_rbailhac_lowmass_PbPb")) {
    printf("Load macro now\n");
    gROOT->LoadMacro(configFilePath.Data());
  }

  // trigger
  TString triggername = "NULL";
  if(trigger == (UInt_t)AliVEvent::kINT7)             triggername = "kINT7";
  else if(trigger == (UInt_t)AliVEvent::kCentral)     triggername = "kCentral";
  else if(trigger == (UInt_t)AliVEvent::kSemiCentral) triggername = "kSemiCentral";
  else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral)) triggername = "kCombinedCentralityTriggers";
  else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kCentral))                           triggername = "kCombinedCentral";
  else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kSemiCentral))                       triggername = "kCombinedSemiCentral";
  
  //create task and add it to the manager (MB)
  TString appendix;
  appendix += TString::Format("Cen%d_%d_%s_rejpileup%d",CenMin,CenMax,triggername.Data(),rejpileup);
  printf("appendix %s\n", appendix.Data());
  AliAnalysisTaskMultiDielectron *task = new AliAnalysisTaskMultiDielectron(Form("MultiDielectron_%s",appendix.Data()));
  task->UsePhysicsSelection();
  task->SetTriggerMask(trigger);

  // PID postcalibration
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

  // PID post-calibration
  Bool_t pidcalib = kFALSE;
  gSystem->Exec(Form("alien_cp %s pidpostcalibration.root",calibFileName.Data()));
  TFile fpid("pidpostcalibration.root");
  if (fpid.IsOpen()){
    pidcalib = kTRUE;
    printf("reading : %s for PID calibration\n",calibFileName.Data());
    hs_mean_TPC_El  = (THnF*)fpid.Get("hs_mean_TPC_El");
    hs_width_TPC_El = (THnF*)fpid.Get("hs_width_TPC_El");
    hs_mean_ITS_El  = (THnF*)fpid.Get("hs_mean_ITS_El");
    hs_width_ITS_El = (THnF*)fpid.Get("hs_width_ITS_El");
    hs_mean_TOF_El  = (THnF*)fpid.Get("hs_mean_TOF_El");
    hs_width_TOF_El = (THnF*)fpid.Get("hs_width_TOF_El");
    hs_mean_TPC_Pi  = (THnF*)fpid.Get("hs_mean_TPC_Pi");
    hs_width_TPC_Pi = (THnF*)fpid.Get("hs_width_TPC_Pi");
    hs_mean_TPC_Ka  = (THnF*)fpid.Get("hs_mean_TPC_Ka");
    hs_width_TPC_Ka = (THnF*)fpid.Get("hs_width_TPC_Ka");
    hs_mean_TPC_Pr  = (THnF*)fpid.Get("hs_mean_TPC_Pr");
    hs_width_TPC_Pr = (THnF*)fpid.Get("hs_width_TPC_Pr");
  }

  // Number of cuts
  const Int_t nDie = (Int_t)gROOT->ProcessLine("GetN()");
  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *diel = reinterpret_cast<AliDielectron*>(gROOT->ProcessLine(Form("Config_rbailhac_lowmass_PbPb(%d,%d)",i,isMC)));
    if(!diel) continue;

    if(pidcalib){
      diel->SetPIDCaibinPU(kTRUE);
      //for electron
      diel->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kElectron,hs_mean_TPC_El ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      diel->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kElectron,hs_width_TPC_El,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      diel->SetCentroidCorrFunctionPU(AliDielectronPID::kITS,AliPID::kElectron,hs_mean_ITS_El ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      diel->   SetWidthCorrFunctionPU(AliDielectronPID::kITS,AliPID::kElectron,hs_width_ITS_El,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      diel->SetCentroidCorrFunctionPU(AliDielectronPID::kTOF,AliPID::kElectron,hs_mean_TOF_El ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      diel->   SetWidthCorrFunctionPU(AliDielectronPID::kTOF,AliPID::kElectron,hs_width_TOF_El,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      //for pion
      diel->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kPion,hs_mean_TPC_Pi     ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      diel->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kPion,hs_width_TPC_Pi    ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      //for kaon
      diel->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kKaon,hs_mean_TPC_Ka     ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      diel->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kKaon,hs_width_TPC_Ka    ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      //for proton
      diel->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kProton,hs_mean_TPC_Pr   ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      diel->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kProton,hs_width_TPC_Pr  ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      
    }

    if(isMix) {
      printf("Add event mixing handler\n");
      AliDielectronMixingHandler *mix = new AliDielectronMixingHandler;
      mix->SetMixType(AliDielectronMixingHandler::kAll);
      mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -8., -6., -4., -2. , 0., 2., 4., 6., 8. , 10.");
      mix->AddVariable(AliDielectronVarManager::kCentralityNew,"0,5,10,30,50,70,90,101");
      mix->AddVariable(AliDielectronVarManager::kNacc,"0,10000");
      mix->SetDepth(Nmix);
      if(!isMC) diel->SetMixingHandler(mix);
    }
    
    task->AddDielectron(diel);
    
  }//loop

  //Add event filter
  task->SetEventFilter((reinterpret_cast<AliDielectronEventCuts*>(gROOT->ProcessLine(Form("GetEventCuts(%f,%f,%d,\"%s\")",(Float_t)CenMin,(Float_t)CenMax,rejpileup,"V0M"))));

  mgr->AddTask(task);

 

  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(Form("tree_lowmass_%s", appendix.Data()),
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         outputFileName);

  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer(Form("Histos_diel_lowmass_%s", appendix.Data()),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);
  
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer(Form("CF_diel_lowmass_%s", appendix.Data()),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer(Form("rbailhac_lowmass_EventStat_%s", appendix.Data()),
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);



    return task;
    
}
