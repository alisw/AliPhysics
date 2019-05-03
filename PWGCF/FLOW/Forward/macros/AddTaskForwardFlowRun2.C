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
  if (nua_file.Contains("alien:")) TGrid::Connect("alien:");
  TFile* file;
  file = TFile::Open(nua_file.Data(), "READ");

  if(!file) { printf("E-AddTaskForwardFlowRun2: Input file with differential weights not found!\n"); return NULL; }

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


AliAnalysisDataContainer* makeWeightContainerSec(TString sec_file, TString containerName){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer* weights;
  if (sec_file.Contains("alien:") ) TGrid::Connect("alien:");
  TFile* file;
  if (sec_file.Contains(containerName))
    file = TFile::Open(sec_file.Data(), "READ");

  if(!file) { printf("E-AddTaskForwardFlowRun2: Input file with secondary weights not found!\n"); return NULL; }

  TList* weights_list = new TList();
  weights_list->SetName("sec_corr");
  
  TH3F* sechist = new TH3F();
  TH3F* nuaforward = new TH3F();

  file->GetObject("correction", sechist);
  sechist->SetDirectory(0);
  sechist->SetNameTitle("sec_corr","sec_corr");

  file->Close();

  weights_list->Add(sechist);

  weights = mgr->CreateContainer(containerName,TList::Class(), AliAnalysisManager::kInputContainer,Form("%s", mgr->GetCommonFileName()));
  weights->SetData(weights_list);
  return weights;
}



AliAnalysisDataContainer* makeEmptyContainer(){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer* weights;

  TList* weights_list = new TList();
  weights_list->SetName("empty");
  
  weights = mgr->CreateContainer("emptyContainer",TList::Class(), AliAnalysisManager::kInputContainer,Form("%s", mgr->GetCommonFileName()));

  return weights;
}


void connectContainer(AliAnalysisDataContainer* container,AliForwardFlowRun2Task* task){

  task->fSettings.nuacentral = static_cast<TH3F*>( static_cast<TList*>(container->GetData())->FindObject("nuacentral") );
  task->fSettings.nuaforward = static_cast<TH3F*>( static_cast<TList*>(container->GetData())->FindObject("nuaforward") );
  task->fSettings.nuacentral->SetDirectory(0);
  task->fSettings.nuaforward->SetDirectory(0);
}

void connectSecContainer(AliAnalysisDataContainer* container,AliForwardFlowRun2Task* task){

  task->fSettings.seccorr_fwd = static_cast<TH3F*>( static_cast<TList*>(container->GetData())->FindObject("sec_corr") );
  task->fSettings.seccorr_fwd->SetDirectory(0);
}


AliAnalysisTaskSE* AddTaskForwardFlowRun2( bool doNUA, 
                                          TString nua_file, 
                                          UShort_t nua_mode, 
                                          bool doetagap, Double_t gap, 
                                          bool mc,  bool esd, bool prim_cen,bool prim_fwd , 
                                          UInt_t tracktype, TString centrality,
                                          Double_t minpt,Double_t maxpt,
                                          TString sec_file_fwd,
                                          TString suffix)
{
  std::cout << "______________________________________________________________________________" << std::endl;

  std::cout << "AddTaskForwardFlowRun2" << std::endl;

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    Fatal("","No analysis manager to connect to.");

  AliForwardFlowRun2Task* task = new AliForwardFlowRun2Task(suffix);

  if (doetagap){
    // if etagap otherwise comment out, and it will be standard
    task->fSettings.fNRefEtaBins = 1;
    task->fSettings.gap = gap;
  }
  else {
    task->fSettings.fNRefEtaBins = 1; // eller skal det være et andet antal?
    task->fSettings.gap = 0.0;
  }

  task->fSettings.use_primaries_cen = prim_cen;
  task->fSettings.use_primaries_fwd = prim_fwd;

  task->fSettings.mc = mc;
  task->fSettings.esd = esd;
  std::cout << "Using tracktype = " << tracktype << std::endl;
  if (tracktype == 0){
    task->fSettings.useSPD = kTRUE;
  }
  else{
    if (tracktype == 768){
      task->fSettings.tracktype = AliForwardSettings::kHybrid;
    }
    else if (tracktype == 128){
      task->fSettings.tracktype = AliForwardSettings::kTPCOnly;
    }
    else if (tracktype == 32){
      task->fSettings.tracktype = AliForwardSettings::kGlobal;
    }
    else if (tracktype == 64){
      task->fSettings.tracktype = AliForwardSettings::kGlobalLoose;
    }
    else if (tracktype == 96){
      task->fSettings.tracktype = AliForwardSettings::kGlobalComb;
    }
    else{
      std::cout << "INVALID TRACK TYPE FOR TPC" << std::endl;
    }
  }

  task->fSettings.minpt = minpt;
  task->fSettings.maxpt = maxpt;
  
  task->fSettings.doNUA = doNUA;
  task->fSettings.nua_mode = nua_mode; 
  task->fSettings.centrality_estimator = centrality; // "V0M";// RefMult08; // "V0M" // "SPDTracklets";

  TString combName = suffix;

  std::cout << "Container name: " << combName << std::endl;
  std::cout << "______________________________________________________________________________" << std::endl;

  mgr->AddTask(task);
  AliAnalysisDataContainer *coutput_recon =
  mgr->CreateContainer(suffix,
   TList::Class(),
   AliAnalysisManager::kOutputContainer,
   mgr->GetCommonFileName());
  mgr->ConnectOutput(task, 1, coutput_recon);
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());

  AliAnalysisDataContainer* valid = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("event_selection_xchange");
  mgr->ConnectInput(task,1,valid);


  if (doNUA){
    TObjArray* taskContainers = mgr->GetContainers();
    AliAnalysisDataContainer* weights;
    
    // mc recon.
    if (nua_file.Contains("reco")){
      weights = (AliAnalysisDataContainer*) taskContainers->FindObject("reco");

      if (!weights) weights = makeWeightContainer(nua_file,"reco");   
      if (nua_file.Contains("reco")) connectContainer( weights, task);
    }
    // mc prim.
    if (nua_file.Contains("prim") ){
      weights = (AliAnalysisDataContainer*) taskContainers->FindObject("prim");

      if (!weights) weights = makeWeightContainer(nua_file,"prim");   
      if (nua_file.Contains("prim")) connectContainer(weights, task);
    }
    
    // global tracks NUA
    if  (nua_file.Contains("cluster70")) 
    {
      weights = (AliAnalysisDataContainer*) taskContainers->FindObject("cluster70");
    
      if (!weights){
        std::cout << "I-AddTaskForwardFlowRun2: Globalz tracks weights not defined - reading now. " << std::endl;
        weights = makeWeightContainer(nua_file,"cluster70");
      }  
      if (nua_file.Contains("cluster70")) connectContainer(weights, task);
    }
        // global tracks NUA
    if  (nua_file.Contains("cluster80")) 
    {
      weights = (AliAnalysisDataContainer*) taskContainers->FindObject("cluster80");
    
      if (!weights){
        std::cout << "I-AddTaskForwardFlowRun2: Globalz tracks weights not defined - reading now. " << std::endl;
        weights = makeWeightContainer(nua_file,"cluster80");
      }  
      if (nua_file.Contains("cluster80")) connectContainer(weights, task);
    }
        // global tracks NUA
    if  (nua_file.Contains("cluster90")) 
    {
      weights = (AliAnalysisDataContainer*) taskContainers->FindObject("cluster90");
    
      if (!weights){
        std::cout << "I-AddTaskForwardFlowRun2: Globalz tracks weights not defined - reading now. " << std::endl;
        weights = makeWeightContainer(nua_file,"cluster90");
      }  
      if (nua_file.Contains("cluster90")) connectContainer(weights, task);
    }
        // global tracks NUA
    if  (nua_file.Contains("cluster100")) 
    {
      weights = (AliAnalysisDataContainer*) taskContainers->FindObject("cluster100");
    
      if (!weights){
        std::cout << "I-AddTaskForwardFlowRun2: Globalz tracks weights not defined - reading now. " << std::endl;
        weights = makeWeightContainer(nua_file,"cluster100");
      }  
      connectContainer(weights, task);
    }


    // global tracks NUA
    if  (nua_file.Contains("globaltrackz")) 
    {
      weights = (AliAnalysisDataContainer*) taskContainers->FindObject("globaltrackz");
    
      if (!weights){
        std::cout << "I-AddTaskForwardFlowRun2: Globalz tracks weights not defined - reading now. " << std::endl;
        weights = makeWeightContainer(nua_file,"globaltrackz");
      }  
      if (nua_file.Contains("globaltrackz")) connectContainer(weights, task);
    }

    // global tracks NUA
    if  (nua_file.Contains("globaltracks")) 
    {
      weights = (AliAnalysisDataContainer*) taskContainers->FindObject("globaltracks");
    
      if (!weights){
        std::cout << "I-AddTaskForwardFlowRun2: Global tracks weights not defined - reading now. " << std::endl;
        weights = makeWeightContainer(nua_file,"globaltracks");
      }  
      if (nua_file.Contains("globaltracks")) connectContainer(weights, task);
    }
    // hybrid tracks NUA
    if  (nua_file.Contains("hybridtracks") ) 
    {
      weights = (AliAnalysisDataContainer*) taskContainers->FindObject("hybridtracks");
    
      if (!weights){
        std::cout << "I-AddTaskForwardFlowRun2: Hybrid tracks weights not defined - reading now. " << std::endl;
        weights = makeWeightContainer(nua_file,"hybridtracks");
      }  
      if (nua_file.Contains("hybridtracks")) connectContainer(weights, task);
    }
    //else printf("E-AddTaskForwardFlowRun2: Invalid detector for NUA weights!\n");
  }

  if (sec_file_fwd != ""){
    TObjArray* taskContainers = mgr->GetContainers();
    AliAnalysisDataContainer* sec_weights;
    
    sec_weights = (AliAnalysisDataContainer*) taskContainers->FindObject("ft_fft");

    if (!sec_weights) sec_weights = makeWeightContainerSec(sec_file_fwd,"ft_fft");   
    connectSecContainer(sec_weights, task);
  }

  return task;
}
/*
 * EOF
 *
 */
