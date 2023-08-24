AliAnalysisTask *AddTask_hmurakam_pp(Bool_t getFromAlien    = kFALSE,
				     TString year           = "16",
				     Bool_t hasSpline       = kFALSE,
				     ULong64_t triggerMask  = AliVEvent::kINT7,
				     Float_t cent_min       = 0.0,
				     Float_t cent_max       = 1.0, //data: "16" -> 0.05, "17","18" -> 0.1
				     TString cFileName      = "Config_hmurakam_pp.C",
				     TString outputFileName = "LMEE.root",
				     Int_t  version          = 0,
				     Bool_t isMC            = kFALSE,
				     Bool_t isMix           = kFALSE,
				     Bool_t rejectPileup    = kTRUE,
				     Bool_t RandomRej       = kFALSE
				     )
{

  //=== get the current analysis manager ===========================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_hmurakam_minbiaspp", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  //TString configBasePath= "./";//for local
  if (!gSystem->AccessPathName(cFileName)) {
    printf("Configfile already present\n");
    configBasePath=Form("%s/",gSystem->pwd());
  }
  else if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/h/hmurakam/PWGDQ/dielectron/macrosLMEE/%s file:./",cFileName.Data()))) ){
    printf("Copy configfile from alien\n");
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);
  std::cout << "Configpath: " << configFilePath << std::endl;
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(cFileName.Data())) {
    printf("Load macro now\n");
    gROOT->LoadMacro(configFilePath.Data());
  }

  //=== Create the main dielectron task =============================
  TString appendix="";
  if(triggerMask == (UInt_t)AliVEvent::kINT7)             appendix = "";
  else if(triggerMask == (UInt_t)AliVEvent::kHighMultV0)  appendix = "_hm";
  printf("appendix %s\n", appendix.Data());
  AliAnalysisTaskMultiDielectron *task = new AliAnalysisTaskMultiDielectron(Form("MultiDielectron_%s",appendix.Data()));
  task->UsePhysicsSelection();
  task->SetTriggerMask(triggerMask);
  //task->SetTriggerOnV0AND(kTRUE); // only for cross-check

  // SPD pile-up rejection in mult. bins
  if (rejectPileup) {
    task->SetRejectPileup(kTRUE);
    task->SetPileupRejTool(AliDielectronEventCuts::kSPDInMultBins);//necessary ? todo check
  }

  // randomize daughters
  task->SetRandomizeDaughters(RandomRej);//default on ? todo:check

  //=== Add event filter ============================================
  task->SetEventFilter((reinterpret_cast<AliDielectronEventCuts*>(gROOT->ProcessLine(Form("GetEventCuts(%f,%f)",(Float_t)cent_min,(Float_t)cent_max)))));

  // PID postcalibration
  TGrid::Connect("alien://");

  TString calibFileNameTOF = "alien:///alice/cern.ch/user/h/hmurakam/PWGDQ/dielectron/calibLMEE/calMaps_TOF.root";
  if(isMC)calibFileNameTOF = "alien:///alice/cern.ch/user/h/hmurakam/PWGDQ/dielectron/calibLMEE/calMaps_TOF_mc.root";//switch data/mc
  TString calibFileNameTPC = "alien:///alice/cern.ch/user/h/hmurakam/PWGDQ/dielectron/calibLMEE/calMaps_TPC.root";

  TH2F* histMean2DTPC  = 0x0;
  TH2F* histWidth2DTPC = 0x0;
  TFile *rootfileTPC   = 0x0;
  if(!isMC){//Data only
    // TPC 2016,17,18[w/ w/o spl]
    gSystem->Exec(Form("alien_cp %s file:./",calibFileNameTPC.Data()));
    printf("reading : %s for TPC PID calibration\n",calibFileNameTPC.Data());
    rootfileTPC = TFile::Open(calibFileNameTPC,"READ");
    if(year =="18" && !hasSpline){
      histMean2DTPC  = (TH2F*)rootfileTPC->Get(Form("m%sns",year.Data()));
      histWidth2DTPC = (TH2F*)rootfileTPC->Get(Form("w%sns",year.Data()));
      printf("2018 no spline\n");
    }
    histMean2DTPC  = (TH2F*)rootfileTPC->Get(Form("m%s",year.Data()));
    histWidth2DTPC = (TH2F*)rootfileTPC->Get(Form("w%s",year.Data()));
    printf("%s and %s\n",Form("m%s",year.Data()),Form("w%s",year.Data()));
    for (Int_t i = 0; i <= histMean2DTPC->GetNbinsX()+1; i++){
      for (Int_t k = 0; k <= histMean2DTPC->GetNbinsY()+1; k++){
	if ( (i == 0) || (k == 0) || (i > histMean2DTPC->GetNbinsX()) || (k > histMean2DTPC->GetNbinsY())) { // under/overflows
	  histMean2DTPC->SetBinContent(i, k, 0.0 );
	  histWidth2DTPC->SetBinContent(i, k, 1.0 );
	}
      }
    }
  }

  TH2F* histMean2DTOF  = 0x0;
  TH2F* histWidth2DTOF = 0x0;
  TFile *rootfileTOF   = 0x0;
  // TOF
  gSystem->Exec(Form("alien_cp %s file:./",calibFileNameTOF.Data()));
  printf("reading : %s for TOF PID calibration\n",calibFileNameTOF.Data());
  rootfileTOF = TFile::Open(calibFileNameTOF,"READ");
  histMean2DTOF  = (TH2F*)rootfileTOF->Get(Form("m%s",year.Data()));
  histWidth2DTOF = (TH2F*)rootfileTOF->Get(Form("w%s",year.Data()));
  printf("%s and %s\n",Form("m%s",year.Data()),Form("w%s",year.Data()));
  for (Int_t i = 0; i <= histMean2DTOF->GetNbinsX()+1; i++){
    for (Int_t k = 0; k <= histMean2DTOF->GetNbinsY()+1; k++){
      if ( (i == 0) || (k == 0) || (i > histMean2DTOF->GetNbinsX()) || (k > histMean2DTOF->GetNbinsY())) { // under/overflows
	histMean2DTOF->SetBinContent(i, k, 0.0 );
	histWidth2DTOF->SetBinContent(i, k, 1.0 );
      }
    }
  }

  // Number of cuts
  const Int_t nDie = (Int_t)gROOT->ProcessLine("GetN()");

  //add dielectron analysis with different cuts to the task
  for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
    AliDielectron *diele = reinterpret_cast<AliDielectron*>(gROOT->ProcessLine(Form("Config_hmurakam_pp(%d,%d,%d)",i,isMC,isMix)));
    if(!diele) continue;
    if(!isMC){//TPC PID Data only
      diele->SetCentroidCorrFunction(histMean2DTPC, AliDielectronVarManager::kP, AliDielectronVarManager::kEta);
      diele->SetWidthCorrFunction(histWidth2DTPC, AliDielectronVarManager::kP, AliDielectronVarManager::kEta);
    }
    //TOF PID
    diele->SetCentroidCorrFunctionTOF(histMean2DTOF, AliDielectronVarManager::kP, AliDielectronVarManager::kEta);
    diele->SetWidthCorrFunctionTOF(histWidth2DTOF, AliDielectronVarManager::kP, AliDielectronVarManager::kEta);

    if(isMix){
      AliDielectronMixingHandler *mix=new AliDielectronMixingHandler;
      printf("Add event mixing handler\n");
      mix->AddVariable(AliDielectronVarManager::kZvPrim,"-10., -7.5, -5., -2.5 , 0., 2.5, 5., 7.5 , 10.");
      mix->SetMixType(AliDielectronMixingHandler::kAll);
      mix->SetDepth(20);
      mix->AddVariable(AliDielectronVarManager::kNacc,"0,10000");
      diele->SetMixingHandler(mix);
    }
    task->AddDielectron(diele);
  }

  mgr->AddTask(task);

  //=== create output containers ===========================
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(Form("tree_lowmass%s",appendix.Data()),
			 TTree::Class(),
			 AliAnalysisManager::kExchangeContainer,
			 outputFileName.Data());

  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer(Form("Output_Histos%s",appendix.Data()),
			 TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 outputFileName.Data());

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer(Form("Output_CF%s",appendix.Data()),
			 TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 outputFileName.Data());

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer(Form("Output_EventStat%s",appendix.Data()),
			 TH1D::Class(),
			 AliAnalysisManager::kOutputContainer,
			 outputFileName.Data());

  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);

  return task;

}

