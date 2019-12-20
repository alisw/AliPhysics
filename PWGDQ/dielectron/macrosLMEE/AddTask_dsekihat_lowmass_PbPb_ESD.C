#ifdef __CLING__
#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>

// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#endif

AliAnalysisTask *AddTask_dsekihat_lowmass_PbPb_ESD(
      Bool_t getFromAlien=kFALSE,
      TString cFileName = "Config_dsekihat_lowmass_PbPb_ESD.C",
      UInt_t trigger = AliVEvent::kINT7|AliVEvent::kCentral|AliVEvent::kSemiCentral,
      const Int_t CenMin =  0,
      const Int_t CenMax = 10,
      const Int_t Nmix   = 10
      )
{

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_dsekihat_lowmass_PbPb_ESD", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
  //TString configBasePath= "./";
  if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/d/dsekihat/PWGDQ/dielectron/macrosLMEE/%s .",cFileName.Data()))) ){
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);

  std::cout << "Configpath:  " << configFilePath << std::endl;

  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);


  TString triggername = "NULL";
  if(trigger == (UInt_t)AliVEvent::kINT7)             triggername = "kINT7";
  else if(trigger == (UInt_t)AliVEvent::kCentral)     triggername = "kCentral";
  else if(trigger == (UInt_t)AliVEvent::kSemiCentral) triggername = "kSemiCentral";
  else if(trigger == (UInt_t)(AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral)) triggername = "kCombinedCentralityTriggers";

  //create task and add it to the manager (MB)
  AliAnalysisTaskMultiDielectron *task = new AliAnalysisTaskMultiDielectron(Form("MultiDielectron_Cen%d_%d_%s",CenMin,CenMax,triggername.Data()));
  task->SetTriggerMask(trigger);

  mgr->AddTask(task);

  //add dielectron analysis with different cuts to the task
  #if defined(__CLING__)
    gROOT->LoadMacro("./Config_dsekihat_lowmass_PbPb_ESD.C");

    const Int_t nDie = (Int_t)gROOT->ProcessLine("GetN()");
    //std::cout << "nDie:  " << nDie << std::endl;

    for (Int_t i=0; i<nDie; ++i){ //nDie as argument of add task
      AliDielectron *diel = reinterpret_cast<AliDielectron*>(gROOT->ProcessLine(Form("Config_dsekihat_lowmass_PbPb_ESD(%d,%d,%d,%d)",i,CenMin,CenMax,Nmix)));
      if(!diel)continue;
      task->AddDielectron(diel);
    }
  #elif defined(__CINT__)
    gROOT->LoadMacro(configFilePath.Data()); //ROOT5 syntax
    const Int_t nDie = GetN();
    for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
      AliDielectron *diel = Config_dsekihat_lowmass_PbPb_ESD(i); // also ROOT5
      if(!diel)continue;
      task->AddDielectron(diel);
    }
  #endif

  // std::cout << "Mixing is set to: " << diel_low->GetMixing() << syd::endl;

  const TString outputFileName = AliAnalysisManager::GetCommonFileName();
  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer(Form("Tree_diel_Cen%d_%d_%s",CenMin,CenMax,triggername.Data()),
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         outputFileName.Data());

  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer(Form("Histos_diel_Cen%d_%d_%s",CenMin,CenMax,triggername.Data()),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName.Data());

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer(Form("CF_diel_Cen%d_%d_%s",CenMin,CenMax,triggername.Data()),
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName.Data());

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer(Form("EventStat_diel_Cen%d_%d_%s",CenMin,CenMax,triggername.Data()),
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

