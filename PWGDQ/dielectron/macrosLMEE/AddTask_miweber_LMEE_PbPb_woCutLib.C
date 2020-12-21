// ROOT6 modifications
#ifdef __CLING__
#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>
#include <AliDielectronVarCuts.h>

// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGDQ/dielectron/macrosLMEE/Config_miweber_LMEE_PbPb_woCutLib.C>

#endif

AliAnalysisTask *AddTask_miweber_LMEE_PbPb_woCutLib(Int_t cutDefinition = 0,
						      TString outputFileName = "AnalysisResult.root",
						      TString directoryBaseName = "miweber_LMEE_PbPb",
						      Bool_t isNano = kFALSE,
						      Bool_t bCutQA = kTRUE,
						      Bool_t useTPCCorr=kFALSE,
						      Bool_t useRotation=kFALSE,
						      Bool_t useMixing=kTRUE,
						      Bool_t noPairing=kFALSE,
						      Bool_t bUsePileUpCutsTPCClusters = kFALSE,
						      Float_t pileUpCutsTPCClustersMin = 0.,
						      Float_t pileUpCutsTPCClustersMax = 0.,
						      
						      Double_t centMin = -1.,
						      Double_t centMax = -1.,
						      Bool_t reqAliEventCuts = kFALSE,
						      Bool_t reqAliEventCutsCorrelated = kFALSE
						    ){


  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_miweber_LMEE_PbPb_woCutLib", "No analysis manager found.");
    return 0;
  }

  Bool_t bESDANA=kFALSE; //Autodetect via InputHandler

  // ROOT6 modifications
#ifndef __CLING__
  TString configBasePath("$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/");
  TString configFile("Config_miweber_LMEE_PbPb_woCutLib.C");
  TString configFilePath(configBasePath+configFile);
   
  //load dielectron configuration files
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFile.Data()))
    gROOT->LoadMacro(configFilePath.Data());

#endif

  //Do we have an MC handler?
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //AOD Usage currently tested with Input handler
  if (mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class()){
    ::Info("AddTask_miweber_LMEE_PbPb_woCutLib", "no dedicated AOD configuration");
  }
  else if (mgr->GetInputEventHandler()->IsA()==AliESDInputHandler::Class()){
    ::Info("AddTask_miweber_LMEE_PbPb_woCutLib","switching on ESD specific code");
    bESDANA=kTRUE;
  }
  
  //create task and add it to the manager
  AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDiEData_miweber_PbPb_woCutLib");

  // for MC no  need for physics selection and for Nano AODs this has been done already
  Int_t triggerNames = AliVEvent::kINT7;//PbPb Min Bias
  if (!hasMC && !isNano){
    task->UsePhysicsSelection();
    task->SetTriggerMask(triggerNames); //PbPb Min Bias, can be set also from outside
    // task->SelectCollisionCandidates(triggerNames); // not needed, since UsePhysicsSelection::UsePhysicsSelection is called; but can be set also from outside
    // task->SetRejectPileup(); // to be done differently (too strong cuts at the moment in dielectron framework) 
  }

  // Note: event cuts are identical for all analysis 'cutDefinition's that run together!
  //Add event filter
  task->SetEventFilter( GetEventCuts(centMin, centMax, reqAliEventCuts, reqAliEventCutsCorrelated) );
  
  //add dielectron analysis with selected cut to the task
  AliDielectron *diel_low = Config_miweber_LMEE_PbPb_woCutLib(cutDefinition,bESDANA,bCutQA,kFALSE,useTPCCorr,useRotation,useMixing,noPairing,hasMC);
  if(diel_low){
    AliDielectronVarCuts *eventplaneCuts = new AliDielectronVarCuts("eventplaneCuts","eventplaneCuts");

    // use event plane cuts only for this cut set
    if(cutDefinition==671){
      eventplaneCuts->AddCut(AliDielectronVarManager::kQnTPCrpH2,-999.,kTRUE); // makes sure that the event has an eventplane
      eventplaneCuts->Print();
      diel_low->GetEventFilter().AddCuts(eventplaneCuts);
    }

    // use pile-up rejection cuts based on TPC clusters/event
    if(bUsePileUpCutsTPCClusters){

      Printf("Using TPC cluster based pile-up cuts with parameter 0 from %f to %f",pileUpCutsTPCClustersMin,pileUpCutsTPCClustersMax);

      TF1* fFitMin = new TF1("fFit","pol6",0,90);
      TF1* fFitMax = new TF1("fFit","pol6",0,90);
      fFitMin->SetParameters(pileUpCutsTPCClustersMin,-95678.946999,2152.010478,-50.119000,0.780528,-0.006150,0.000019);
      fFitMax->SetParameters(pileUpCutsTPCClustersMax,-95678.946999,2152.010478,-50.119000,0.780528,-0.006150,0.000019);
      
      AliDielectronEventCuts *pileUpCuts = new AliDielectronEventCuts("pileUpCuts","pileUpCuts");
      pileUpCuts->SetMinCorrCutFunction(fFitMin, AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kNTPCclsEvent);
      pileUpCuts->SetMaxCorrCutFunction(fFitMax, AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kNTPCclsEvent);

      pileUpCuts->Print();
      diel_low->GetEventFilter().AddCuts(pileUpCuts);
    }
      
    task->AddDielectron(diel_low);
    printf("successfully added AliDielectron: %s\n",diel_low->GetName());
  }
  else{
    Printf("=======================================");
    Printf("No AliDielectron object loaded -> EXIT ");
    Printf("=======================================");
    return NULL;
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
  //                         "miweber_LMEE_PbPb_CF.root");
  
  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer(Form("%s_EventStat",directoryBaseName.Data()),
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName.Data());
  
  
  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  if((AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ZDCEPExchangeContainer"))
    mgr->ConnectInput(task,  1,(AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ZDCEPExchangeContainer"));
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
}
