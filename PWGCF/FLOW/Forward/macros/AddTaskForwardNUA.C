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

AliAnalysisDataContainer* makeWeightContainer(TString nua_file, TString containerName){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer* weights;
  if (nua_file.Contains("alien:")) {
    std::cout << "I-AddTaskForwardSecondaries: Connecting to alien" << std::endl;
    TGrid::Connect("alien:");
  }
  std::cout << nua_file << std::endl;
  TFile* file;
  file = TFile::Open(nua_file.Data(), "READ");

  if(!file) { printf("E-AddTaskForwardSecondaries: Input file with differential weights not found!\n"); return NULL; }

  TList* weights_list = new TList();
  weights_list->SetName("nuaWeights");
  
  TH3F* nuacentral = new TH3F();
  TH3F* nuaforward = new TH3F();

  file->GetObject("nuacentral", nuacentral);
  nuacentral->SetDirectory(0);
  nuacentral->SetNameTitle("nuacentral","nuacentral");

  file->GetObject("nuaforward", nuaforward);
  nuaforward->SetDirectory(0);
  nuaforward->SetNameTitle("nuaforward","nuaforward");
  file->Close();

  weights_list->Add(nuacentral);
  weights_list->Add(nuaforward);

  weights = mgr->CreateContainer(containerName,TList::Class(), AliAnalysisManager::kInputContainer,Form("%s", mgr->GetCommonFileName()));
  weights->SetData(weights_list);
  return weights;
}

void connectContainer(AliAnalysisDataContainer* container,AliForwardNUATask* task){

  task->fSettings.correct_nua_mc = static_cast<TH3F*>( static_cast<TList*>(container->GetData())->FindObject("nuaforward") );
  task->fSettings.correct_nua_mc->SetDirectory(0);
}


AliAnalysisDataContainer* makeEmptyContainer(){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer* weights;

  TList* weights_list = new TList();
  weights_list->SetName("empty");
  
  weights = mgr->CreateContainer("emptyContainer",TList::Class(), AliAnalysisManager::kInputContainer,Form("%s", mgr->GetCommonFileName()));

  return weights;
}


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
  task->fSettings.doNUE = kFALSE;
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
  Bool_t doNUA = kFALSE;
  TString nua_file =  "/home/thoresen/Documents/PhD/plots/nua/Datasets/LHC17i2f/fmdreco_tpcreco_runsz.root";

  if (doNUA){
    std::cout << "doing NUA" << std::endl;
    TObjArray* taskContainers = mgr->GetContainers();
    AliAnalysisDataContainer* weights;
    
    TObjArray *tx = nua_file.Tokenize("/");
    TObjArray *ty = ((TObjString *)(tx->At(tx->GetEntries()-1)))->String().Tokenize(".");
    TString nuaobject =  ((TObjString *)(ty->At(0)))->String();
    std::cout << "I-AddTaskForwardNUA: NUA container " << nuaobject << std::endl;

    weights = (AliAnalysisDataContainer*) taskContainers->FindObject(nuaobject);
    
    if (!weights) {
      std::cout << "I-AddTaskForwardNUA: " << nuaobject << " weights not defined - reading now. " << std::endl;
      weights = makeWeightContainer(nua_file,nuaobject);
    }
    connectContainer( weights, task);
  }
  return task;
}
/*
 * EOF
 *
 */
