
AliAnalysisTask* AddTask_acapon(TString outputFileName = "AnalysisResult.root",
                                TString names          = "kCutSet1",
                                Bool_t SDDstatus       = kTRUE,  // SDD used during data taking?
                                Bool_t hasMC           = kFALSE,
                                Int_t wagonNum         = 0,      // Needed when running multiple wagons
                                Int_t pairing          = 0,      // 0=No pairing, 1=Pairing, 2=Pairing+pairCuts
                                Bool_t doEventMixing   = kTRUE,  // Event mixing for R factor
                                Int_t rejectPileUp     = 0,      // 0=None, 1=SPD, 2=SPDinMultBins, 4=MultiVertexer
                                Bool_t trackVarPlots   = kTRUE,  // Simple track QA plots
                                Int_t whichDetPlots    = 0,      // 0=None,1=ITS,2=TPC,4=TOF,7=All3
                                // Use PID post calibration maps (for electrons)
                                Int_t usePIDcorrMaps   = 4,      // 0=None, 1=ITS,2=TPC,3=TOF,4=All Three
                                // Option to use AliEventCuts class for additional event cuts
                                Int_t whichAliEvtCuts  = 0,      // 0=None, 1=Use, 2=Also use correlation cuts
                                Bool_t plots3D         = kFALSE,
                                Bool_t v0plots         = kFALSE,  // Plots for PID calibration
                                Bool_t getFromAlien    = kFALSE) // Pull config+CutLib from alien directory
{

    TObjArray* arrNames = names.Tokenize(";");
    Int_t nDie = arrNames->GetEntries();
    Printf("Number of implemented cuts: %i", nDie);

    // Set appropriate flags
    // To which detectors ePID correction maps will be applied
    Bool_t useITScorr = ((usePIDcorrMaps == 1) || (usePIDcorrMaps == 4)) ? kTRUE : kFALSE;
    Bool_t useTPCcorr = ((usePIDcorrMaps == 2) || (usePIDcorrMaps == 4)) ? kTRUE : kFALSE;
    Bool_t useTOFcorr = ((usePIDcorrMaps == 3) || (usePIDcorrMaps == 4)) ? kTRUE : kFALSE;
    // Pairing flags
    Bool_t doPairing     = ((pairing == 1) || (pairing == 2)) ? kTRUE : kFALSE;
    Bool_t applyPairCuts = (pairing == 2)                     ? kTRUE : kFALSE;


    // Print out analysis control settings
    // Many options not yet implemented for "noCutLib", so no need
    // to print them out.
    std::cout << "Output file name: " << outputFileName << std::endl;
    std::cout << "Cut set         : " << names          << std::endl;
    std::cout << "Use SDD         : " << SDDstatus      << std::endl;
    std::cout << "Monte Carlo     : " << hasMC          << std::endl;
    std::cout << "Wagon number    : " << wagonNum       << std::endl;
    std::cout << "Pairing         : " << doPairing      << std::endl;
    std::cout << "Pair cuts       : " << applyPairCuts  << std::endl;
    std::cout << "Event mixing    : " << doEventMixing  << std::endl;
    std::cout << "rejPileUp       : " << rejectPileUp   << std::endl;
    std::cout << "Track plots     : " << trackVarPlots  << std::endl;
    std::cout << "Which det plots : " << whichDetPlots  << std::endl;
    std::cout << "Use ITScorr     : " << useITScorr     << std::endl;
    std::cout << "Use TPCcorr     : " << useTPCcorr     << std::endl;
    std::cout << "Use TOFcorr     : " << useTOFcorr     << std::endl;
    std::cout << "Use AliEventCuts: " << whichAliEvtCuts<< std::endl;
    std::cout << "3D plots        : " << plots3D        << std::endl;
    std::cout << "v0 plots        : " << v0plots        << std::endl;

    // Get the current analysis manager
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if(!mgr){
      ::Error("AddTask_acapon", "No analysis manager found.");
      return 0;
    }

    TString configBasePath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/");
    TString configLMEECutLib("LMEECutLib_acapon.C");
    TString configFile = "Config_acapon.C";

    // Determine if ESDs or AODs are being analysed
    // CutLibrary version only works/tested for AODs
    // noCutLib version is mixed depending on cut setting
    if(mgr->GetInputEventHandler()->IsA() == AliAODInputHandler::Class()){
      ::Info("AddTask_acapon", "AOD configuration");
    }
    else if(mgr->GetInputEventHandler()->IsA() == AliESDInputHandler::Class()){
      ::Info("AddTask_acapon","ESD configuration");
      configLMEECutLib = "LMEECutLib_acapon_ESD.C";
    }
    ::Info("AddTask_acapon",Form("Use LMeeCutLib: %s",configLMEECutLib.Data()));


    // Load updated macros from private ALIEN path
    TString myConfig = "alien_cp alien:///alice/cern.ch/user/a/acapon/PWGDQ/dielectron/macrosLMEE/Config_acapon.C .";
    TString myCutLib = "alien_cp alien:///alice/cern.ch/user/a/acapon/PWGDQ/dielectron/macrosLMEE/LMEECutLib_acapon.C .";
    if(getFromAlien && (!gSystem->Exec(myConfig)) && (!gSystem->Exec(myCutLib))){
      std::cout << "Copy config from Alien" << std::endl;
      configBasePath = TString::Format("%s/",gSystem->pwd());
    }

    TString configFilePath(configBasePath+configFile);
    TString configLMEECutLibPath(configBasePath+configLMEECutLib);

    // Load dielectron configuration files
    if(!gROOT->GetListOfGlobalFunctions()->FindObject(configLMEECutLib.Data())){
      gROOT->LoadMacro(configLMEECutLibPath.Data());
    }
    if(!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data())){
      gROOT->LoadMacro(configFilePath.Data());
    }


    // Create task and add it to the manager
    AliAnalysisTaskMultiDielectron* task = new AliAnalysisTaskMultiDielectron(::Form("DielectronTask%d", wagonNum));
    if(!task){
      ::Error("AddTask_acapon", "MultiDielectron trask not created");
      return 0x0;
    }

    // Add event filter (the same for all cut sets and configs)
    Int_t triggerNames = (AliVEvent::kINT7);
    task->SetTriggerMask(triggerNames);
    task->UsePhysicsSelection();
    if(rejectPileUp != 0){
      task->SetRejectPileup(kTRUE);
      if(rejectPileUp == 1){
        task->SetPileupRejTool(AliDielectronEventCuts::kSPD);
      }
      else if(rejectPileUp == 2){
        task->SetPileupRejTool(AliDielectronEventCuts::kSPDInMultBins);
      }
      else if(rejectPileUp == 3){
        task->SetPileupRejTool(AliDielectronEventCuts::kMultiVertexer);
      }
    }
    // Set correct flags for AliEventCuts
    Bool_t reqAliEvtCuts     = (whichAliEvtCuts == 0) ? kFALSE : kTRUE;
    Bool_t reqAliEvtCutsCorr = (whichAliEvtCuts == 2) ? kTRUE  : kFALSE;
    // Event cuts are the same for all cut settings regardless
    // of config setup
    LMEECutLib* cutLib = new LMEECutLib(SDDstatus);
    task->SetEventFilter(cutLib->GetEventCuts(reqAliEvtCuts, reqAliEvtCutsCorr));

    // Add the task to the manager
    mgr->AddTask(task);
    // Add dielectron analysis with different cuts to the task
    for(Int_t i = 0; i < nDie; ++i){
      TString dielTaskName(arrNames->At(i)->GetName());
      AliDielectron* diel_low = Config_acapon(dielTaskName, hasMC, SDDstatus,
                                              doPairing, applyPairCuts, doEventMixing,
                                              trackVarPlots, whichDetPlots, v0plots,
                                              useITScorr, useTPCcorr, useTOFcorr,
                                              plots3D);
      if(!diel_low){
        continue;
      }
      task->AddDielectron(diel_low);
      Printf("successfully added AliDielectron: %s\n",diel_low->GetName());
    }// End cut settings initialisation loop


    // Create output container
    AliAnalysisDataContainer* coutput1 =
    mgr->CreateContainer(::Form("acapon_tree_%d",wagonNum),
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         outputFileName.Data());

    AliAnalysisDataContainer* cOutputHist1 =
    mgr->CreateContainer(::Form("acapon_out_%d", wagonNum),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName.Data());

    AliAnalysisDataContainer* cOutputHist2 =
    mgr->CreateContainer(::Form("acapon_CF_%d", wagonNum),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName.Data());

    AliAnalysisDataContainer* cOutputHist3 =
    mgr->CreateContainer(::Form("acapon_EventStat_%d", wagonNum),
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
