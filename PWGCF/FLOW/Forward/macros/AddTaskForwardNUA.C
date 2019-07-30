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

AliAnalysisTaskSE* AddTaskForwardNUA(UShort_t nua_mode, bool makeFakeHoles, bool mc,  bool esd,bool prim_cen,bool prim_fwd , 
                                     Int_t tracktype, TString centrality,Double_t minpt,Double_t maxpt,TString suffix="")
{
  std::cout << "______________________________________________________________________________" << std::endl;

  std::cout << "AddTaskForwardNUA" << std::endl;

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    Fatal("","No analysis manager to connect to.");

  AliForwardNUATask* task = new AliForwardNUATask(suffix);
  TString resName = suffix;

    task->fSettings.use_primaries_cen = prim_cen;
    if (mc) resName += (prim_cen ? "_primcen" : "_trcen");

    task->fSettings.use_primaries_fwd = prim_fwd;
    if (mc) resName += (prim_fwd ? "_primfwd" : "_trfwd");

  task->fSettings.mc = mc;
  task->fSettings.esd = esd;
  std::cout << "Using tracktype = " << tracktype << std::endl;
  if (tracktype == 0){
    resName += "_SPD";
  }
  else{
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

  //TString name; name.Form("%d",type);
   //return name.Data();

//TString combinedName;
//combinedName.Form("%s%s", suffix);
TString ptmaxname;
TString ptminname;
ptminname.Form("%d",(int)(minpt*10));
ptmaxname.Form("%d",(int)(maxpt*10));
    task->fSettings.minpt = minpt;
    resName += "_minpt";
    resName +=  ptminname;//std::to_string

    task->fSettings.maxpt = maxpt;
    resName += "_maxpt";
    resName += ptmaxname;//std::to_string


  resName += "_" + centrality;
  if (mc) resName += "_mc";

  if (makeFakeHoles){
    task->fSettings.makeFakeHoles = kTRUE;
    resName += "_fakeholes";
  }

  task->fSettings.centrality_estimator = centrality; // "V0M";// RefMult08; // "V0M" // "SPDTracklets";

  if (nua_mode == AliForwardSettings::kInterpolate)
    resName += "_NUA_interpolate";
  if (nua_mode == AliForwardSettings::kFill)
    resName += "_NUA_fill";
  if (nua_mode == AliForwardSettings::kNormal)
    resName += "_NUA_normal";

  task->fSettings.nua_mode = nua_mode; // "V0M";// RefMult08; // "V0M" // "SPDTracklets";

  TString combName = resName + '_' + suffix;
  std::cout << "Container name: " << combName << std::endl;

  //resName = "hej";
  AliAnalysisDataContainer* valid = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("event_selection_xchange");
  task->ConnectInput(1,valid);

  AliAnalysisDataContainer *coutput_recon =
  mgr->CreateContainer(suffix, //combName
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
