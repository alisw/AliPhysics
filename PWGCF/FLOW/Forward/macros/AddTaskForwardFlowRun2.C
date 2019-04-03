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
AliAnalysisTaskSE* AddTaskForwardFlowRun2( bool doNUA, bool makeFakeHoles, 
                                          TString nua_file, TString ref_nua_file, 
                                          UShort_t nua_mode, 
                                          bool doetagap, Double_t gap, 
                                          bool mc,  bool esd, bool prim_cen,bool prim_fwd , 
                                          UInt_t tracktype, TString centrality,
                                          Double_t minpt,Double_t maxpt,
                                          TString sec_file_cen,TString sec_file_fwd,
                                          TString suffix)
{
  std::cout << "______________________________________________________________________________" << std::endl;

  std::cout << "AddTaskForwardFlowRun2" << std::endl;

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    Fatal("","No analysis manager to connect to.");

  TString name = suffix;

  AliForwardFlowRun2Task* task = new AliForwardFlowRun2Task(name);
  TString resName = "ForwardFlow";

  if (doetagap){
    // if etagap otherwise comment out, and it will be standard
    task->fSettings.fFlowFlags = task->fSettings.kEtaGap;
    task->fSettings.fNRefEtaBins = 1;
    task->fSettings.gap = gap;
    resName += "_etagap";
//    resName += std::to_string((int)(10*gap));

    //resName += std::to_string(gap);
  }
  else {
    task->fSettings.fNRefEtaBins = 1; // eller skal det være et andet antal?
    task->fSettings.gap = 0.0;
  }

    if (makeFakeHoles){
      task->fSettings.makeFakeHoles = kTRUE;
      resName += "_fakeholes";
    }

  task->fSettings.use_primaries_cen = prim_cen;
  if (mc) resName += (prim_cen ? "_primcen" : "_trcen");

  task->fSettings.use_primaries_fwd = prim_fwd;
  if (mc) resName += (prim_fwd ? "_primfwd" : "_trfwd");

  task->fSettings.mc = mc;
  task->fSettings.esd = esd;
  std::cout << "Using tracktype = " << tracktype << std::endl;
  if (tracktype == 0){
    task->fSettings.fFlowFlags = task->fSettings.kSPD;
    task->fSettings.useSPD = kTRUE;
    resName += "_SPD";
  }
  else{
    task->fSettings.fFlowFlags = task->fSettings.kTPC;
    if (tracktype == 768){
      task->fSettings.tracktype = AliForwardSettings::kHybrid;
      resName += "_hybrid";
    }
    else if (tracktype == 128){
      task->fSettings.tracktype = AliForwardSettings::kTPCOnly;
      resName += "_TPConly";
    }
    else if (tracktype == 32){
      task->fSettings.tracktype = AliForwardSettings::kGlobal;
      resName += "_global";
    }
    else if (tracktype == 64){
      task->fSettings.tracktype = AliForwardSettings::kGlobalLoose;
      resName += "_globalLoose";
    }
    else if (tracktype == 96){
      task->fSettings.tracktype = AliForwardSettings::kGlobalComb;
      resName += "_globalComb";
    }
    else{
      std::cout << "INVALID TRACK TYPE FOR TPC" << std::endl;
    }
  }


  task->fSettings.minpt = minpt;
  task->fSettings.maxpt = maxpt;
  

  task->fSettings.doNUA = doNUA;
  if (task->fSettings.doNUA){
    if (nua_file.Contains("alien:")) TGrid::Connect("alien:");


    TFile *file = TFile::Open(nua_file);//new TFile(nua_file);


    // if (nua_mode == AliForwardSettings::kInterpolate)
    //   resName += "_NUA_interpolate";
    // if (nua_mode == AliForwardSettings::kFill)
    //   resName += "_NUA_fill";
    // if (nua_mode == AliForwardSettings::kNormal)
    //   resName += "_NUA_normal";

    task->fSettings.nua_mode = nua_mode; 

    file->GetObject("nuacentral", task->fSettings.nuacentral);
    task->fSettings.nuacentral->SetDirectory(0);
    file->GetObject("nuaforward", task->fSettings.nuaforward);
    task->fSettings.nuaforward->SetDirectory(0);
    file->Close();

    TFile *file_ref = TFile::Open(ref_nua_file);
    file_ref->GetObject("nuacentral", task->fSettings.nuacentral_ref);
    task->fSettings.nuacentral_ref->SetDirectory(0);
    file_ref->GetObject("nuaforward", task->fSettings.nuaforward_ref);
    task->fSettings.nuaforward_ref->SetDirectory(0);

    file_ref->Close();

  }

  if (sec_file_fwd != ""){
    TFile *file1 = new TFile(sec_file_fwd);

    file1->GetObject("correction", task->fSettings.seccorr_fwd);
    task->fSettings.seccorr_fwd->SetDirectory(0);
    file1->Close();
  }


  if (sec_file_cen != ""){
    TFile *file2 = new TFile(sec_file_cen);

    file2->GetObject("correction", task->fSettings.seccorr_cen);
    task->fSettings.seccorr_cen->SetDirectory(0);
    file2->Close();
  }

  resName += "_" + centrality;
  if (mc) resName += "_mc";

  task->fSettings.centrality_estimator = centrality; // "V0M";// RefMult08; // "V0M" // "SPDTracklets";
   TString combName = suffix;

  std::cout << "Container name: " << combName << std::endl;
  std::cout << "______________________________________________________________________________" << std::endl;

  AliAnalysisDataContainer *coutput_recon =
  mgr->CreateContainer(suffix,//combName,
   TList::Class(),
   AliAnalysisManager::kOutputContainer,
   mgr->GetCommonFileName());


  AliAnalysisDataContainer* valid = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("event_selection_xchange");
  task->ConnectInput(1,valid);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput_recon);

  return task;
}
/*
 * EOF
 *
 */
