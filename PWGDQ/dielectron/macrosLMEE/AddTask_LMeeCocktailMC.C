void AddTask_LMeeCocktailMC(Int_t CollisionSystem = 200, Float_t MaxEta = 0.8, Float_t MinPt = 0.2, Float_t MaxPt = 8.0, Bool_t WriteTTree = kFALSE, Int_t ResolType = 2 , Int_t ALTweightType = 1, TString resFileName = "") {

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
  task= new AliAnalysisTaskLMeeCocktailMC(Form("LMeeCocktailMC_%1.2f",MaxEta));
  task->SetCollisionSystem(CollisionSystem);
  task->SetMaxEta(MaxEta);
  task->SetMinPt(MinPt);
  task->SetMaxPt(MaxPt);
  task->SetWriteTTree(WriteTTree);
  task->SetResolType(ResolType);
  task->SetALTweight(ALTweightType);
  if(resFileName != ""){
    Printf("Set resolution file name to %s",resFileName.Data());
    task->SetResFileName(resFileName);
  }
  
  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("LMeeCocktailMC_%1.2f",MaxEta), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:LMeeCocktailMC",AliAnalysisManager::GetCommonFileName()));
    
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  
  return;
}
