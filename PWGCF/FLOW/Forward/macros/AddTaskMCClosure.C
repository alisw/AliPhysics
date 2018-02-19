/**
 * @file   FTAddMyTask.C
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 * 
 * @brief  Add Q-cummulant forward task to train 
 * 
 * 
 * @ingroup pwglf_forward_scripts_tasks
 */
/** 
 * @defgroup pwglf_forward_flow Flow 
 *
 * Code to deal with flow 
 *
 * @ingroup pwglf_forward_topical
 */
/** 
 * Add Flow task to train 
 * 
 * @ingroup pwglf_forward_flow
 */
AliAnalysisTaskSE* AddTaskForwardMCClosure()
{
  Bool_t etagap = false;
  Int_t mode = kRECON;
  bool doNUA = false;

  std::cout << "AddTaskForwardMCClosure" << std::endl;

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
    Fatal("","No analysis manager to connect to.");
  
  const char* name = Form("ForwardFlowQC");
  AliForwardMCClosure* task = new AliForwardMCClosure(name);
  
  TString resName = "awesome";

  task->fSettings.doNUA = doNUA;


  if (task->fSettings.doNUA){

    //TString nua_filepath = std::getenv("NUA_FILE");
    //if (!nua_filepath) {
      TString nua_filepath = "/home/thoresen/Documents/PhD/analysis/nua.root";
      std::cerr << "Environment variable 'NUA_FILE' not found (this should be a path to nua.root).\n";
      std::cerr << "   Using default value: '" << nua_filepath << "'\n";
    //}

    TFile *file = new TFile(nua_filepath);

    file->GetObject("nuacentral", task->fSettings.nuacentral);  

    task->fSettings.nuacentral->SetDirectory(0);
    file->GetObject("nuaforward", task->fSettings.nuaforward);   
    task->fSettings.nuaforward->SetDirectory(0);
    file->Close(); 
  }

    task->fSettings.fNDiffEtaBins = 20;

  if (etagap){
    // if etagap otherwise comment out, and it will be standard
    task->fSettings.fFlowFlags = task->fSettings.kEtaGap;
    task->fSettings.fNRefEtaBins = 1;
    task->fSettings.gap = 0.5;
  }
  else {
    task->fSettings.fNRefEtaBins = 1; // eller skal det være et andet antal?
  }

  if (mode == kRECON) {
    AliAnalysisDataContainer *coutput_recon =
    mgr->CreateContainer(resName,
     TList::Class(),
     AliAnalysisManager::kOutputContainer,
     mgr->GetCommonFileName());
    task->fSettings.fDataType = task->fSettings.kRECON;
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutput_recon);
  }
  else if (mode == kTRUTH) {
    AliAnalysisDataContainer *coutput_truth =
    mgr->CreateContainer(resName,
     TList::Class(),
     AliAnalysisManager::kOutputContainer,
     mgr->GetCommonFileName());

    task->fSettings.fDataType = task->fSettings.kMCTRUTH;
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutput_truth);
  }
  else {
    ::Error("AddTaskForwardMCClosure", "Invalid mode specified");
  }


  return task;
}
/*
 * EOF
 *
 */
