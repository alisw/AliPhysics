//void AddTask_LMeeCocktailMC(Double_t maxy = 0.8) {
void AddTask_LMeeCocktailMC(Int_t collisionSystem = 200, Float_t maxy = 0.8, Float_t minpt = 0.2, Bool_t WriteTTree = kFALSE) {

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
    Error("AddTask_LMeeCocktailMC", "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  //            find input container
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================

  AliAnalysisTaskLMeeCocktailMC *task=NULL;
  task= new AliAnalysisTaskLMeeCocktailMC(Form("LMeeCocktailMC_%1.2f",maxy));
  task->SetMaxY(maxy);
  task->SetMinPt(minpt);
  task->SetCollisionSystem(collisionSystem);
  task->SetWriteTTree(WriteTTree);
  
  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("LMeeCocktailMC_%1.2f",maxy), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:LMeeCocktailMC",AliAnalysisManager::GetCommonFileName()));
    
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  
  return;
}
// null
