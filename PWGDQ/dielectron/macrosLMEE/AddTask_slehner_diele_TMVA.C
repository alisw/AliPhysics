AliAnalysisTask *AddTask_slehner_diele_TMVA(  Double_t centMin=0.,
                                              Double_t centMax=100.,
                                              Bool_t SetTPCCorrection=kFALSE,
                                              Bool_t useAODFilterCuts=kFALSE,
                                              TString TMVAweight = "TMVAClassification_BDTG.weights_094.xml" ){
  
  TString directoryBaseName = "slehner_LMEE_TMVA";
  TString outputFileName = "AnalysisResults.root";
  Bool_t isNano=kFALSE;
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_slehner_diele_TMVA", "No analysis manager found.");
    return 0;
  }

  Bool_t bESDANA=kFALSE; //Autodetect via InputHandler

  TString configBasePath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/");
  TString configFile("Config_slehner_diele_TMVA.C");
  TString configLMEECutLib("LMEECutLib_slehner.C");
  TString configLMEECutLibPath(configBasePath+configLMEECutLib);
  TString configFilePath(configBasePath+configFile);
     
  //load dielectron configuration files
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data()))    gROOT->LoadMacro(configFilePath.Data());
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configLMEECutLibPath.Data()))    gROOT->LoadMacro(configLMEECutLibPath.Data());

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDiEData_slehner_TMVA");

  task->UsePhysicsSelection();

  Int_t triggerNames = AliVEvent::kINT7;//PbPb Min Bias, can be set also from outside
  if(!isNano){
    task->SelectCollisionCandidates(triggerNames);
    task->SetTriggerMask(triggerNames);
    // task->SetRejectPileup(); // to be done differently (too strong cuts at the moment in dielectron framework) 
  }

  LMEECutLib* cutlib = new LMEECutLib();
  Bool_t isRun2 = kTRUE;
  
  task->SetEventFilter(cutlib->GetEventCuts(centMin, centMax)); // All cut sets have same event cuts

  mgr->AddTask(task);
  
  //add dielectron analyses with various cuts to the task
  
  for(Int_t TrCut = 0; TrCut <= 0; ++TrCut){
    for(Int_t PIDCut = 0; PIDCut <= 0; ++PIDCut){
      for(Int_t MVACut = 9; MVACut <= 9; ++MVACut){
        AliDielectron * diel_low = Config_slehner_diele_TMVA(TrCut,PIDCut,MVACut,kFALSE);
        if(!diel_low){
          Printf("=======================================");
          Printf("No AliDielectron object loaded -> EXIT ");
          Printf("=======================================");
          return NULL; 
        }    

        std::cout << "CutTr: "<<TrCut<<" CutPID: "<<PIDCut<<" MVAcut:"<<MVACut<<" being added"<< std::endl;
        diel_low->GetTrackFilter()->AddCuts(SetupTrackCutsAndSettings(TrCut, PIDCut, MVACut, useAODFilterCuts,TMVAweight));   
        task->AddDielectron(diel_low);
        printf("successfully added AliDielectron: %s\n",diel_low->GetName());  
      }
    }
  }

  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(Form("%s_tree",directoryBaseName.Data()),
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         outputFileName.Data());
  
  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer(Form("%s_out",directoryBaseName.Data()),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName.Data());
  
  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer(Form("%s_CF",directoryBaseName.Data()),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName.Data());
  //                         "slehner_diele_PbPb_CF.root");
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer(Form("%s_EventStat",directoryBaseName.Data()),
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
