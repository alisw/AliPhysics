// ROOT6 modifications
#ifdef __CLING__
#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>
#include <AliDielectronVarCuts.h>

// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#endif

AliAnalysisTask *AddTask_jjung_lowmass(Bool_t getFromAlien=kFALSE,
                                      TString cFileName = "Config_jjung_lowmass.C",
                                      TString outputFileName="LMEE.root",
                                      Int_t wagonnr=0,
                                      UInt_t triggerMask = AliVEvent::kCentral				      
                                      )
{
  std::cout << "AddTask_jjung_lowmass" << std::endl;
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_jjung_lowmass", "No analysis manager found.");
    return 0;
  }

  //Base Directory for GRID / LEGO Train
  TString configBasePath= "/home/jerome/analysis/localLegotrain/005_tests_PbPbAOD/";
  if (!gSystem->AccessPathName(cFileName))
  {
    printf("file already present\n");
    configBasePath=Form("%s/",gSystem->pwd());
  }
  else if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/j/jjung/%s file:./",cFileName.Data()))) ){
    std::cout << "ALIEN???!  " << std::endl;
    configBasePath=Form("%s/",gSystem->pwd());
  }

  TString configFilePath(configBasePath+cFileName);

  std::cout << "Configpath:  " << configFilePath << std::endl;
  
  //Do we have an MC handler?
  Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  if (!gROOT->GetListOfGlobalFunctions()->FindObject(cFileName.Data()))
  gROOT->LoadMacro(configFilePath.Data());  //old root5

  //create task and add it to the manager (MB)
  AliAnalysisTaskMultiDielectron *task = new AliAnalysisTaskMultiDielectron(Form("MultiDielectron_%d", wagonnr));
  if (!hasMC) task->UsePhysicsSelection();
  task->SetTriggerMask(UInt_t(triggerMask));
  task->SelectCollisionCandidates(UInt_t(triggerMask));

  std::cout << "Task created" << std::endl;

  task->SetRandomizeDaughters(kTRUE); //default kFALSE
  std::cout << "Adding event cuts" << std::endl;

  mgr->AddTask(task);
  std::cout << "Task Added" << std::endl;

  const Int_t nDie = (Int_t)gROOT->ProcessLine("GetN()");
  //add dielectron analysis with different cuts to the task
  #if defined(__CLING__)
    Int_t dot = cFileName.First('.'); 
    Int_t len = cFileName.Length(); 
    if (len > 0)
      cFileName.Remove(dot,len-dot);
    //TMacro conf_die(gSystem->ExpandPathName(configFilePath.Data())); //ROOT6
    for (Int_t i=0; i<nDie; ++i){ //nDie as argument of add task
      //AliDielectron *diel_low = reinterpret_cast<AliDielectron *>(conf_die.Exec(Form("%i",i)));
      AliDielectron *diel_low = reinterpret_cast<AliDielectron*>(gROOT->ProcessLine(Form("%s(%i,%i,%i)",cFileName.Data(),i,kFALSE,wagonnr)));
      if(!diel_low)continue;
      task->AddDielectron(diel_low);
    }
  #elif defined(__CINT__)
    gROOT->LoadMacro(configFilePath.Data()); //ROOT5 syntax
    for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
      AliDielectron *diel_low = Config_jjung_lowmass(i); // also ROOT5
      if(!diel_low)continue;
      task->AddDielectron(diel_low);
    }
  #endif

  std::cout << "Cuts Settings Added" << std::endl;

  //create output container
  AliAnalysisDataContainer *coutput1 = 0x0;
  AliAnalysisDataContainer *cOutputHist1 = 0x0;
  AliAnalysisDataContainer *cOutputHist2 = 0x0;
  AliAnalysisDataContainer *cOutputHist3 = 0x0;


  if(wagonnr == 0){
    coutput1 = mgr->CreateContainer("tree_lowmass", TTree::Class(),
                         AliAnalysisManager::kExchangeContainer, outputFileName.Data());
  
    cOutputHist1 = mgr->CreateContainer("Histos_diel_lowmass", TList::Class(),
                         AliAnalysisManager::kOutputContainer, outputFileName.Data());

    cOutputHist2 = mgr->CreateContainer("CF_diel_lowmass", TList::Class(),
                         AliAnalysisManager::kOutputContainer, outputFileName.Data());

    cOutputHist3 = mgr->CreateContainer("jjung_lowmass_EventStat", TH1D::Class(),
                         AliAnalysisManager::kOutputContainer, outputFileName.Data());
  }

  else{
    coutput1 = mgr->CreateContainer(Form("tree_lowmass_%d", wagonnr), TTree::Class(),
                         AliAnalysisManager::kExchangeContainer, outputFileName.Data());
  
    cOutputHist1 = mgr->CreateContainer(Form("Histos_diel_lowmass_%d", wagonnr), TList::Class(),
                         AliAnalysisManager::kOutputContainer, outputFileName.Data());

    cOutputHist2 = mgr->CreateContainer(Form("CF_diel_lowmass_%d", wagonnr), TList::Class(),
                         AliAnalysisManager::kOutputContainer, outputFileName.Data());

    cOutputHist3 = mgr->CreateContainer(Form("jjung_lowmass_EventStat_%d", wagonnr), TH1D::Class(),
                         AliAnalysisManager::kOutputContainer, outputFileName.Data());

  }

  mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, coutput1 );
  mgr->ConnectOutput(task, 1, cOutputHist1);
  mgr->ConnectOutput(task, 2, cOutputHist2);
  mgr->ConnectOutput(task, 3, cOutputHist3);
  
  return task;
  
}

