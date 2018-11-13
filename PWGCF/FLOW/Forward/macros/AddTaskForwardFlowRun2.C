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
AliAnalysisTaskSE* AddTaskForwardFlowRun2(bool nua_mode, bool doetagap, bool doNUA, Double_t gap, bool mc, bool prim, Int_t tracktype, TString centrality, bool esd)
{
  std::cout << "______________________________________________________________________________" << std::endl;

  std::cout << "AddTaskForwardFlowRun2" << std::endl;

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    Fatal("","No analysis manager to connect to.");

  AliForwardFlowRun2Task* task = new AliForwardFlowRun2Task("hej");
  TString resName = "ForwardFlow";

  if (doetagap){
    // if etagap otherwise comment out, and it will be standard
    task->fSettings.fFlowFlags = task->fSettings.kEtaGap;
    task->fSettings.fNRefEtaBins = 1;
    task->fSettings.gap = gap;
    resName += "_etagap";
    //resName += std::to_string(gap);
  }
  else {
    task->fSettings.fNRefEtaBins = 1; // eller skal det være et andet antal?
    task->fSettings.gap = 0.0;
  }

  task->fSettings.mc = mc;
  std::cout << "Using tracktype = " << tracktype << std::endl;
  if (tracktype == 0){
    task->fSettings.fFlowFlags = task->fSettings.kSPD;
    resName += "_SPD";
  }
  else{
    task->fSettings.fFlowFlags = task->fSettings.kTPC;
    if (tracktype == 768){
      task->fSettings.tracktype = AliForwardFlowRun2Settings::kHybrid;
      resName += "_hybrid";
    }
    else if (tracktype == 128){
      task->fSettings.tracktype = AliForwardFlowRun2Settings::kTPCOnly;
      resName += "_TPConly";
    }
    else if (tracktype == 32){
      task->fSettings.tracktype = AliForwardFlowRun2Settings::kGlobal;
      resName += "_global";
    }
    else if (tracktype == 64){
      task->fSettings.tracktype = AliForwardFlowRun2Settings::kGlobalLoose;
      resName += "_globalLoose";
    }
    else if (tracktype == 96){
      task->fSettings.tracktype = AliForwardFlowRun2Settings::kGlobalComb;
      resName += "_globalComb";
    }
    else{
      std::cout << "INVALID TRACK TYPE FOR TPC" << std::endl;
    }
  }

  task->fSettings.doNUA = doNUA;
  if (task->fSettings.doNUA){
      //TString nua_filepath = std::getenv("NUA_FILE");
      //if (!nua_filepath) {
    TString nua_filepath = "/home/thoresen/Documents/PhD/Analysis/nua.root";
    // std::cerr << "Environment variable 'NUA_FILE' not found (this should be a path to nua.root).\n";
    // std::cerr << "   Using default value: '" << nua_filepath << "'\n";
    //}
    TFile *file;

    if (nua_mode) {
      // INTERPOLATE
      resName += "_NUA_extrapolated";
        file = new TFile("/home/thoresen/Programs/alice/AliPhysics/PWGCF/FLOW/Forward/corrections/LHC15o_pass5_lowIR_selected.root");
    }
    else {
      // FILL
        file = new TFile("/home/thoresen/Programs/alice/AliPhysics/PWGCF/FLOW/Forward/corrections/LHC15o_pass5_lowIR_selected.root");
        resName += "_NUA_filled";
    }

    file->GetObject("nuacentral", task->fSettings.nuacentral);
    task->fSettings.nuacentral->SetDirectory(0);
    file->GetObject("nuaforward", task->fSettings.nuaforward);
    task->fSettings.nuaforward->SetDirectory(0);
    file->Close();
  }

  resName += "_" + centrality;
  if (mc) resName += "_mc";

  task->fSettings.centrality_estimator = centrality; // "V0M";// RefMult08; // "V0M" // "SPDTracklets";
  std::cout << "Container name: " << resName << std::endl;
  //resName = "hej";
  std::cout << "______________________________________________________________________________" << std::endl;

  AliAnalysisDataContainer *coutput_recon =
  mgr->CreateContainer(resName,
   TList::Class(),
   AliAnalysisManager::kOutputContainer,
   mgr->GetCommonFileName());

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput_recon);

  return task;
}
/*
 * EOF
 *
 */
