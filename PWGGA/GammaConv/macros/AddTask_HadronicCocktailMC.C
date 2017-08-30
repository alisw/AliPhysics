void AddTask_HadronicCocktailMC(Int_t particleFlag = 0, Bool_t runLightOutput = kFALSE, TString maxyset = "0.8") {

  Double_t maxy = maxyset.Atof();
  maxy         /= 100;  // needed to enable subwagon feature on grid

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
  task = new AliAnalysisTaskHadronicCocktailMC(Form("HadronicCocktailMC_%1.2f",maxy));
  task->SetMaxY(maxy);
  task->SetLightOutput(runLightOutput);
  task->SetAnalyzedParticle(particleFlag);          // switch to run: 0 - pi0, 1 - eta, 2 - pi+-
  
  TString                   analyzedParticle = "";
  if (particleFlag==0)      analyzedParticle = "pi0";
  else if (particleFlag==1) analyzedParticle = "eta";
  else if (particleFlag==2) analyzedParticle = "pi+-";
  
  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("HadronicCocktailMC_%s_%1.2f",analyzedParticle.Data(),maxy), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:HadronicCocktailMC",AliAnalysisManager::GetCommonFileName()));
  
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  
  return;
  
}
