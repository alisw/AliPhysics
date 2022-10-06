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
    const TString outname = "LMEE.root",
    const Bool_t isMC = kFALSE
    )
{
  //Add a task AliAnalysisTaskTagAndProbe to the analysis train
  //Author: Daiki Sekihata

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask_dsekihat_TAP_PbPb", "No analysis manager to connect to");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask_dsekihat_TAP_PbPb", "This task requires an input event handler");
    return NULL;
  }

  //load config files and cut library

  //TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  TString configBasePath= "./";
  if(getFromAlien
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

  TString triggername = "NULL";
  if(trigger == (UInt_t)AliVEvent::kINT7)             triggername = "kINT7";
  else if(trigger == (UInt_t)AliVEvent::kCentral)     triggername = "kCentral";
  else if(trigger == (UInt_t)AliVEvent::kSemiCentral) triggername = "kSemiCentral";
  else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral)) triggername = "kCombinedCentralityTriggers";
  else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kCentral))                           triggername = "kCombinedCentral";
  else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kSemiCentral))                       triggername = "kCombinedSemiCentral";

  TString taskname = Form("TAP_%sEff_Cen%d_%d_%s",type.Data(),CenMin,CenMax,triggername.Data());
  AliAnalysisTaskTagAndProbe* task = new AliAnalysisTaskTagAndProbe(taskname);
  gROOT->GetListOfSpecials()->Add(task);//this is only for ProcessLine(Config_dsekihat(taskname));

  task->SelectCollisionCandidates(trigger);
  task->SetTriggerMask(trigger);
  task->SetCentralityMin(CenMin);
  task->SetCentralityMax(CenMax);
  task->SetDepthNMixed(Nmix);
  //centrality setting
  task->SetCentralityEstimator("V0M");

  task->SetMC(isMC);
  //task->SetPhivCutRange(0.05,2.0);

  task->GetEventFilter()->AddCuts(reinterpret_cast<AliDielectronEventCuts*>(gROOT->ProcessLine(Form("LMEECutLib::SetupEventCuts(%f,%f,%d,\"%s\")",(Float_t)CenMin,(Float_t)CenMax,kTRUE,"V0M"))));//kTRUE is for Run2

  gROOT->ProcessLine(Form("Config_dsekihat_TAP_PbPb(%s,\"%s\",%d)",task->GetName(),type.Data(),isMC));

  //PID calibration
  TFile *rootfile = 0x0;

  if(calibFileName.Contains("MC")){//for MC
    printf("reading : %s for PID calibration\n",calibFileName.Data());
    rootfile = TFile::Open(calibFileName,"READ");
    task->SetPIDCalibinPU(kFALSE);
    if(rootfile && rootfile->IsOpen()){
      TH3D *h3mean_ITS  = (TH3D*)rootfile->Get("h3mean_ITS");
      TH3D *h3width_ITS = (TH3D*)rootfile->Get("h3width_ITS");
      h3mean_ITS ->SetDirectory(0);
      h3width_ITS->SetDirectory(0);
      task->SetCentroidCorrFunctionITS(h3mean_ITS,  AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta);
      task->SetWidthCorrFunctionITS(   h3width_ITS, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta);

      TH3D *h3mean_TOF = (TH3D*)rootfile->Get("h3mean_TOF");
      TH3D *h3width_TOF = (TH3D*)rootfile->Get("h3width_TOF");
      h3mean_TOF ->SetDirectory(0);
      h3width_TOF->SetDirectory(0);
      task->SetCentroidCorrFunctionTOF(h3mean_TOF, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta);
      task->SetWidthCorrFunctionTOF(  h3width_TOF, AliDielectronVarManager::kNSDDSSDclsEvent, AliDielectronVarManager::kPIn,  AliDielectronVarManager::kEta);

      rootfile->Close();
    }
  }
  else if(calibFileName != ""){//for data in PU
    printf("reading : %s for PID calibration\n",calibFileName.Data());
    rootfile = TFile::Open(calibFileName,"READ");
    THnF *hs_mean_TPC_El  = (THnF*)rootfile->Get("hs_mean_TPC_El");
    THnF *hs_width_TPC_El = (THnF*)rootfile->Get("hs_width_TPC_El");
    THnF *hs_mean_ITS_El  = (THnF*)rootfile->Get("hs_mean_ITS_El");
    THnF *hs_width_ITS_El = (THnF*)rootfile->Get("hs_width_ITS_El");
    THnF *hs_mean_TOF_El  = (THnF*)rootfile->Get("hs_mean_TOF_El");
    THnF *hs_width_TOF_El = (THnF*)rootfile->Get("hs_width_TOF_El");
    THnF *hs_mean_TPC_Pi  = (THnF*)rootfile->Get("hs_mean_TPC_Pi");
    THnF *hs_width_TPC_Pi = (THnF*)rootfile->Get("hs_width_TPC_Pi");
    THnF *hs_mean_TPC_Ka  = (THnF*)rootfile->Get("hs_mean_TPC_Ka");
    THnF *hs_width_TPC_Ka = (THnF*)rootfile->Get("hs_width_TPC_Ka");
    THnF *hs_mean_TPC_Pr  = (THnF*)rootfile->Get("hs_mean_TPC_Pr");
    THnF *hs_width_TPC_Pr = (THnF*)rootfile->Get("hs_width_TPC_Pr");

    if(rootfile && rootfile->IsOpen()){
      task->SetPIDCalibinPU(kTRUE);
      //for electron
      task->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kElectron,hs_mean_TPC_El ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      task->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kElectron,hs_width_TPC_El,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      task->SetCentroidCorrFunctionPU(AliDielectronPID::kITS,AliPID::kElectron,hs_mean_ITS_El ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      task->   SetWidthCorrFunctionPU(AliDielectronPID::kITS,AliPID::kElectron,hs_width_ITS_El,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      task->SetCentroidCorrFunctionPU(AliDielectronPID::kTOF,AliPID::kElectron,hs_mean_TOF_El ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      task->   SetWidthCorrFunctionPU(AliDielectronPID::kTOF,AliPID::kElectron,hs_width_TOF_El,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      //for pion
      task->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kPion,hs_mean_TPC_Pi     ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      task->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kPion,hs_width_TPC_Pi    ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      //for kaon
      task->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kKaon,hs_mean_TPC_Ka     ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      task->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kKaon,hs_width_TPC_Ka    ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      //for proton
      task->SetCentroidCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kProton,hs_mean_TPC_Pr   ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);
      task->   SetWidthCorrFunctionPU(AliDielectronPID::kTPC,AliPID::kProton,hs_width_TPC_Pr  ,AliDielectronVarManager::kNSDDSSDclsEvent,AliDielectronVarManager::kTPCpileupZ,AliDielectronVarManager::kTPCpileupM,AliDielectronVarManager::kPIn,AliDielectronVarManager::kEta);

    }
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

