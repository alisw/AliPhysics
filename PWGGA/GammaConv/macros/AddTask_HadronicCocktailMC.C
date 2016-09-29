void AddTask_HadronicCocktailMC(Double_t maxy = 0.8, Bool_t runPi0 = kTRUE, Bool_t runLightOutput = kFALSE) {
  
  // ================= Load Librariers =================================
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCDB");
  gSystem->Load("libSTEER");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libTender");
  gSystem->Load("libTenderSupplies");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
  gSystem->Load("libPWGGAGammaConv");

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_HadronicCocktailMC", "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  //            find input container
  AliAnalysisTaskHadronicCocktailMC *task=NULL;
  task= new AliAnalysisTaskHadronicCocktailMC(Form("HadronicCocktailMC_%1.2f",maxy));
  task->SetMaxY(maxy);
  task->SetLightOutput(runLightOutput);
  task->SetAnalyzePi0(runPi0);          // switch to run for pi0/eta
  
  TString analyzedParticle;
  if(runPi0)analyzedParticle  = "pi0";
  else analyzedParticle       = "eta";
  
  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("HadronicCocktailMC_%s_%1.2f",analyzedParticle.Data(),maxy), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:HadronicCocktailMC",AliAnalysisManager::GetCommonFileName()));
  
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  
  return;
  
}
