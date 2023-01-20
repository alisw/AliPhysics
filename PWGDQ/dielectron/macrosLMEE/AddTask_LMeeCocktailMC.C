#include <RVersion.h>
AliAnalysisTaskLMeeCocktailMC *AddTask_LMeeCocktailMC(Int_t CollisionSystem = 200, Float_t MaxEta = 0.8, Float_t MinPt = 0.2, Float_t MaxPt = 8.0, Bool_t WriteTTree = kFALSE, Int_t ResolType = 2 , Bool_t local = kFALSE, Int_t ALTweightType = 1, TString resFileName = "",TString effFileName = "", Int_t version = 0,Float_t NMee = 1200, Float_t MinMee = 0, Float_t MaxMee = 6,Float_t NPtee  = 400, Float_t MinPtee = 0, Float_t MaxPtee = 10, Float_t MinOpAng = 0.0, Bool_t DoRapidityCut = kFALSE) {

// Only load libraries manually for ROOT 5. In ROOT 6, loading libraries
// multiple times leads to segfaults.
#if ROOT_VERSION_CODE < ROOT_VERSION(5, 99, 0)
  // gSystem->Load(...) lines go here
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
#endif

  // ================= Load Librariers =================================
  AliAnalysisTaskLMeeCocktailMC *task=NULL;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_LMeeCocktailMC", "No analysis manager found.");
    return task;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  //            find input container
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================

 
  task= new AliAnalysisTaskLMeeCocktailMC(Form("LMeeCocktailMC_%1.2f_%d",MaxEta,version));
  task->SetCollisionSystem(CollisionSystem);
  task->SetMaxEta(MaxEta);
  task->SetMinPt(MinPt);
  task->SetMaxPt(MaxPt);
  //mee
  task->SetNBinsMee(NMee);
  task->SetMinMee(MinMee);
  task->SetMaxMee(MaxMee);
  //ptee
  task->SetNBinsPtee(NPtee);
  task->SetMinPtee(MinPtee);
  task->SetMaxPtee(MaxPtee);

  task->SetMinOpAng(MinOpAng);

  //Rapidiy cut on mother particle
  task->SetRapidityCut(DoRapidityCut);

  task->SetWriteTTree(WriteTTree);
  task->SetResolType(ResolType);
  task->SetResFileLocal(local);
  task->SetALTweight(ALTweightType);

  // resolution file to be set always
  Printf("Set resolution file name to %s",resFileName.Data());
  task->SetResFileName(resFileName);

  // efficiency file to be set always
  Printf("Set eff file name to %s",effFileName.Data());
  task->SetEffFileName(effFileName);

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("LMeeCocktailMC_%1.2f_%d",MaxEta,version), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:LMeeCocktailMC",AliAnalysisManager::GetCommonFileName()));
    
  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  
  return task;
}
