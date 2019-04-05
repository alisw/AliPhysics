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

  if (doetagap){
    // if etagap otherwise comment out, and it will be standard
    task->fSettings.fFlowFlags = task->fSettings.kEtaGap;
    task->fSettings.fNRefEtaBins = 1;
    task->fSettings.gap = gap;
  }
  else {
    task->fSettings.fNRefEtaBins = 1; // eller skal det være et andet antal?
    task->fSettings.gap = 0.0;
  }

  if (makeFakeHoles){
    task->fSettings.makeFakeHoles = kTRUE;
  }

  task->fSettings.use_primaries_cen = prim_cen;
  task->fSettings.use_primaries_fwd = prim_fwd;

  task->fSettings.mc = mc;
  task->fSettings.esd = esd;
  std::cout << "Using tracktype = " << tracktype << std::endl;
  if (tracktype == 0){
    task->fSettings.fFlowFlags = task->fSettings.kSPD;
    task->fSettings.useSPD = kTRUE;
  }
  else{
    task->fSettings.fFlowFlags = task->fSettings.kTPC;
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
  mgr->CreateContainer(suffix,//combName,
   TList::Class(),
   AliAnalysisManager::kOutputContainer,
   mgr->GetCommonFileName());
  mgr->ConnectOutput(task, 1, coutput_recon);

  AliAnalysisDataContainer* valid = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("event_selection_xchange");
  task->ConnectInput(1,valid);

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());


  //if (false){
  if (task->fSettings.doNUA){

    TObjArray* taskContainers = mgr->GetContainers();
    AliAnalysisDataContainer* weights;
    if (nua_file.Contains("ITS") || ref_nua_file.Contains("ITS")) 
    {
      weights = (AliAnalysisDataContainer*) taskContainers->FindObject("ITSclusters");
    
      if (!weights){
        std::cout << "I-AddTaskForwardFlowRun2: Weights not defined - reading now. " << std::endl;
        if (nua_file.Contains("alien:") || ref_nua_file.Contains("alien:")) TGrid::Connect("alien:");
        TFile* file;
        if (nua_file.Contains("ITS"))
          file = TFile::Open(nua_file.Data(), "READ");
        else
          file = TFile::Open(ref_nua_file.Data(), "READ");

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

        weights = mgr->CreateContainer("ITSclusters",TList::Class(), AliAnalysisManager::kInputContainer,Form("%s", mgr->GetCommonFileName()));
        weights->SetData(weights_list);
      }
      if (nua_file.Contains("ITS")) task->ConnectInput(2,weights);
      if (ref_nua_file.Contains("ITS")) task->ConnectInput(3,weights);
    }
    else if  (nua_file.Contains("SPD")||ref_nua_file.Contains("SPD")) 
    {
      weights = (AliAnalysisDataContainer*) taskContainers->FindObject("SPDtracklets");
    
      if (!weights){
        std::cout << "I-AddTaskForwardFlowRun2: Weights not defined - reading now. " << std::endl;
        if (nua_file.Contains("alien:") || ref_nua_file.Contains("alien:")) TGrid::Connect("alien:");
        TFile* file;
        if (nua_file.Contains("SPD"))
          file = TFile::Open(nua_file.Data(), "READ");
        else
          file = TFile::Open(ref_nua_file.Data(), "READ");

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

        weights = mgr->CreateContainer("SPDtracklets",TList::Class(), AliAnalysisManager::kInputContainer,Form("%s", mgr->GetCommonFileName()));
        weights->SetData(weights_list);
      }
      if (nua_file.Contains("SPD")) task->ConnectInput(2,weights);
      if (ref_nua_file.Contains("SPD")) task->ConnectInput(3,weights);
    }
    else if  (nua_file.Contains("global")||ref_nua_file.Contains("global")) 
    {
      weights = (AliAnalysisDataContainer*) taskContainers->FindObject("globaltracks");
    
      if (!weights){
        std::cout << "I-AddTaskForwardFlowRun2: Weights not defined - reading now. " << std::endl;
        if (nua_file.Contains("alien:") || ref_nua_file.Contains("alien:")) TGrid::Connect("alien:");
        TFile* file;
        if (nua_file.Contains("global"))
          file = TFile::Open(nua_file.Data(), "READ");
        else
          file = TFile::Open(ref_nua_file.Data(), "READ");

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

        weights = mgr->CreateContainer("globaltracks",TList::Class(), AliAnalysisManager::kInputContainer,Form("%s", mgr->GetCommonFileName()));
        weights->SetData(weights_list);
      }  
      if (nua_file.Contains("global")) task->ConnectInput(2,weights);
      if (ref_nua_file.Contains("global")) task->ConnectInput(3,weights);
    }
    else printf("E-AddTaskForwardFlowRun2: Invalid central detector for NUA weights!\n");
  }
  return task;
}
/*
 * EOF
 *
 */
