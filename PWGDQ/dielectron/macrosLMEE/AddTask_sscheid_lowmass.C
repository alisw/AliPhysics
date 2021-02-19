// ROOT6 modifications
#ifdef __CLING__
#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>
#include <AliDielectronVarCuts.h>

// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#endif

AliAnalysisTask *AddTask_sscheid_lowmass(Bool_t getFromAlien=kFALSE,
                                             TString cFileName = "Config_sscheid_lowmass_06.C",
                                             int nDie = 0,
                                             TString outputFileName="LMEE.root",
                                             ULong64_t triggerMask = AliVEvent::kINT7
                                             )
{

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_sscheid_lowmass", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train
  TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";

  if (!gSystem->AccessPathName(cFileName))
  {
    printf("file already present\n");
    configBasePath=Form("%s/",gSystem->pwd());
  }
  else if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/h/hscheid/PWGDQ/dielectron/macrosLMEE/%s .",cFileName.Data()))) )
  {
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);

  std::cout << "Configpath:  " << configFilePath << std::endl;

  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //if (!gROOT->GetListOfGlobalFunctions()->FindObject(cFileName.Data()))

  //create task and add it to the manager (MB)
  AliAnalysisTaskMultiDielectron *task = new AliAnalysisTaskMultiDielectron("MultiDielectron");
  if (!hasMC) task->UsePhysicsSelection();
  task->SetTriggerMask(triggerMask);
  //  taskMB->SetRejectPileup();

  // task->SetRandomizeDaughters(randomizeDau); //default kFALSE

  //Add event filter
  // task->SetEventFilter( GetEventCuts() );

  mgr->AddTask(task);

  //add dielectron analysis with different cuts to the task
  #if defined(__CLING__)
    TMacro conf_die(gSystem->ExpandPathName(configFilePath.Data())); //ROOT6
    for (Int_t i=0; i<nDie; ++i){ //nDie as argument of add task
      AliDielectron *diel_low = reinterpret_cast<AliDielectron *>(conf_die.Exec(Form("%i",i)));
      if(!diel_low)continue;
      task->AddDielectron(diel_low);
    }
  #elif defined(__CINT__)
    gROOT->LoadMacro(configFilePath.Data()); //ROOT5 syntax
    for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
      AliDielectron *diel_low = Config_sscheid_lowmass(i); // also ROOT5
      if(!diel_low)continue;
      task->AddDielectron(diel_low);
    }
  #endif


  // std::cout << "Mixing is set to: " << diel_low->GetMixing() << syd::endl;


  //create output container
  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("tree_lowmass",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         outputFileName.Data());

  AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("Histos_diel_lowmass",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName.Data());

  AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("CF_diel_lowmass",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName.Data());

  AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("sscheid_lowmass_EventStat",
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
