void AddTask_ClusterQA( TString   V0ReaderEventCutNumber        = "00000003",
                        TString   V0ReaderPhotonCutNumber       = "060000084001001500000000",
                        TString   TaskEventCutnumber            = "00000003",
                        TString   TaskClusterCutnumberEMC       = "1111100017032220000",
                        TString   TaskClusterCutnumberDMC       = "1111100017032220000",
                        TString   TaskMesonCutnumber            = "0163300000000000",
                        Int_t     minNLM                        = 1,
                        Int_t     maxNLM                        = 1,
                        Bool_t    isMC                          = kFALSE,
                        Int_t     IsHeavyIon                    = 0,
                        Bool_t    kHistograms                   = kTRUE,
                        Double_t  kTree                         = 1.0,  // 0. / 0 / kFALSE for no, 1. / 1 / kTRUE for yes,  x > 1.0 will use only 1/x of the event statistics for the tree
                        TString   V0ReaderCutNumberAODBranch    = "0000000060084001001500000",
                        Bool_t    doEtaShiftV0Reader            = kFALSE,
                        Bool_t    enableV0findingEffi           = kFALSE,              // enables V0finding efficiency histograms
                        TString   periodNameV0Reader            = ""
                    ){
  
  Bool_t enableTriggerOverlapRej = kTRUE;
  Float_t   maxFacPtHard                  = 3.;
  TString corrTaskSetting = "";
  Int_t     enableExtMatchAndQA           = 5;
  
  Bool_t doSaveSurroundingCells                   = 1;
  Int_t nSurroundingCellsSaved      = 12;
  Bool_t doSaveClusterCells               = 1;
  Bool_t doSaveEventProp                  = 1;
  
  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMC);
  }

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",V0ReaderEventCutNumber.Data(),V0ReaderPhotonCutNumber.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());

    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
    fV0ReaderV1->SetProduceV0FindingEfficiency(enableV0findingEffi);

    AliConvEventCuts *fEventCuts=NULL;
    if(V0ReaderEventCutNumber!=""){
      fEventCuts= new AliConvEventCuts(V0ReaderEventCutNumber.Data(),V0ReaderEventCutNumber.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
      if(fEventCuts->InitializeCutsFromCutString(V0ReaderEventCutNumber.Data())){
        fV0ReaderV1->SetEventCuts(fEventCuts);
        fEventCuts->SetFillCutHistograms("",kTRUE);
        if (IsHeavyIon==2){
          fEventCuts->SelectCollisionCandidates(AliVEvent::kINT7);
          fEventCuts->DoEtaShift(doEtaShiftV0Reader);
        }
      }
    }

    // Set AnalysisCut Number
    AliConversionPhotonCuts *fCuts=NULL;
    if(V0ReaderPhotonCutNumber!=""){
      fCuts= new AliConversionPhotonCuts(V0ReaderPhotonCutNumber.Data(),V0ReaderPhotonCutNumber.Data());
      fCuts->SetPreSelectionCutFlag(kTRUE);
      fCuts->SetIsHeavyIon(IsHeavyIon);
      fCuts->SetV0ReaderName(V0ReaderName);
      if(fCuts->InitializeCutsFromCutString(V0ReaderPhotonCutNumber.Data())){
        fV0ReaderV1->SetConversionCuts(fCuts);
        fCuts->SetFillCutHistograms("",kTRUE);
      }
    }


    if(inputHandler->IsA()==AliAODInputHandler::Class()){
      // AOD mode
      fV0ReaderV1->AliV0ReaderV1::SetDeltaAODBranchName(Form("GammaConv_%s_gamma",V0ReaderCutNumberAODBranch.Data()));
    }
    fV0ReaderV1->Init();

    AliLog::SetGlobalLogLevel(AliLog::kInfo);

    //connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);

  } else {
    Error("AddTask_V0ReaderV1", "Cannot execute AddTask, V0ReaderV1 already exists.");
  }

  AliConvEventCuts *analysisEventCuts = new AliConvEventCuts();
  analysisEventCuts->SetV0ReaderName(V0ReaderName);
  
  analysisEventCuts->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
  analysisEventCuts->SetMaxFacPtHard(maxFacPtHard);
  analysisEventCuts->SetCorrectionTaskSetting(corrTaskSetting);
  if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts->SetPeriodEnum(periodNameV0Reader);  
  analysisEventCuts->InitializeCutsFromCutString(TaskEventCutnumber.Data());
  analysisEventCuts->SetFillCutHistograms("",kFALSE);

  TString caloCutPosEMC = TaskClusterCutnumberEMC;
  caloCutPosEMC.Resize(1);
  TString TrackMatcherNameEMC = Form("CaloTrackMatcher_%s",caloCutPosEMC.Data());
  if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherNameEMC.Data()) ){
    AliCaloTrackMatcher* fTrackMatcherEMC = new AliCaloTrackMatcher(TrackMatcherNameEMC.Data(),caloCutPosEMC.Atoi());
    fTrackMatcherEMC->SetV0ReaderName(V0ReaderName);
    //fTrackMatcherEMC->SetCorrectionTaskSetting(corrTaskSetting);
    mgr->AddTask(fTrackMatcherEMC);
    mgr->ConnectInput(fTrackMatcherEMC,0,cinput);
  }
  TString caloCutPosDMC = TaskClusterCutnumberDMC;
  caloCutPosDMC.Resize(1);
  TString TrackMatcherNameDMC = Form("CaloTrackMatcher_%s",caloCutPosDMC.Data());
  if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherNameDMC.Data()) ){
    AliCaloTrackMatcher* fTrackMatcherDMC = new AliCaloTrackMatcher(TrackMatcherNameDMC.Data(),caloCutPosDMC.Atoi());
    fTrackMatcherDMC->SetV0ReaderName(V0ReaderName);
    //fTrackMatcherDMC->SetCorrectionTaskSetting(corrTaskSetting);
    mgr->AddTask(fTrackMatcherDMC);
    mgr->ConnectInput(fTrackMatcherDMC,0,cinput);
  }
  
  AliCaloPhotonCuts *analysisClusterCutsEMC = new AliCaloPhotonCuts();
  analysisClusterCutsEMC->SetV0ReaderName(V0ReaderName);
  //analysisClusterCutsEMC->SetCorrectionTaskSetting(corrTaskSetting);
  analysisClusterCutsEMC->SetCaloTrackMatcherName(TrackMatcherNameEMC);
  analysisClusterCutsEMC->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsEMC->InitializeCutsFromCutString(TaskClusterCutnumberEMC.Data());
  analysisClusterCutsEMC->SetFillCutHistograms("");
  
  AliCaloPhotonCuts *analysisClusterCutsDMC = new AliCaloPhotonCuts();
  analysisClusterCutsDMC->SetV0ReaderName(V0ReaderName);
  //analysisClusterCutsDMC->SetCorrectionTaskSetting(corrTaskSetting);
  analysisClusterCutsDMC->SetCaloTrackMatcherName(TrackMatcherNameEMC);
  analysisClusterCutsDMC->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsDMC->InitializeCutsFromCutString(TaskClusterCutnumberDMC.Data());
  analysisClusterCutsDMC->SetFillCutHistograms("");

  AliConversionMesonCuts *analysisMesonCuts    = new AliConversionMesonCuts();
  analysisMesonCuts->SetEnableOpeningAngleCut(kFALSE);
  analysisMesonCuts->SetIsMergedClusterCut(1);
  analysisMesonCuts->InitializeCutsFromCutString(TaskMesonCutnumber.Data());
  analysisMesonCuts->SetFillCutHistograms("");


  AliAnalysisTaskClusterQA *fQA = new AliAnalysisTaskClusterQA(Form("%s_%s_%s_QA",TaskEventCutnumber.Data(),TaskClusterCutnumberEMC.Data(),TaskClusterCutnumberEMC.Data()));
  fQA->SetEventCuts(analysisEventCuts,IsHeavyIon);
  fQA->SetClusterCutsEMC(analysisClusterCutsEMC,IsHeavyIon);
  fQA->SetClusterCutsDMC(analysisClusterCutsDMC,IsHeavyIon);
  fQA->SetMesonCuts(analysisMesonCuts,IsHeavyIon);
  fQA->SetMinMaxNLMCut(minNLM,maxNLM);
  fQA->FillType(kTree,kHistograms);
  fQA->SetIsMC(isMC);
  if(isMC)
    fQA->SetSaveMCInformation(kTRUE);
  fQA->SetSaveSurroundingCells(doSaveSurroundingCells);
  fQA->SetNSurroundingCells(nSurroundingCellsSaved);
  fQA->SetSaveClusterCells(doSaveClusterCells);
  fQA->SetSaveEventProperties(doSaveEventProp);
  fQA->SetV0ReaderName(V0ReaderName);
  mgr->AddTask(fQA);

  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("GammaCaloQA_%s_%s_%s", TaskEventCutnumber.Data(), TaskClusterCutnumberEMC.Data(),TaskClusterCutnumberEMC.Data()), TList::Class(), AliAnalysisManager::kOutputContainer,"ClusterTree.root");
  mgr->ConnectOutput(fQA,  1, coutput);
  mgr->ConnectInput(fQA,   0, cinput);

  return;
}

