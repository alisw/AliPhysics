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

void connectContainer(AliAnalysisDataContainer* container,AliForwardSecondariesTask* task){

  task->fSettings.nuacentral = static_cast<TH3F*>( static_cast<TList*>(container->GetData())->FindObject("nuacentral") );
  task->fSettings.nuaforward = static_cast<TH3F*>( static_cast<TList*>(container->GetData())->FindObject("nuaforward") );
  task->fSettings.nuacentral->SetDirectory(0);
  task->fSettings.nuaforward->SetDirectory(0);
}


AliAnalysisDataContainer* makeEmptyContainer(){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer* weights;

  TList* weights_list = new TList();
  weights_list->SetName("empty");
  
  weights = mgr->CreateContainer("emptyContainer",TList::Class(), AliAnalysisManager::kInputContainer,Form("%s", mgr->GetCommonFileName()));

  return weights;
}

AliAnalysisTaskSE* AddTaskForwardSecondaries(TString taskname)
{
  std::cout << "AddTaskForwardSecondaries" << std::endl;


  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    Fatal("","No analysis manager to connect to.");

  const char* name = Form("ForwardFlowQC");
  AliForwardSecondariesTask* task = new AliForwardSecondariesTask(name);
  TString resName = taskname;

  //TString nua_file =  "/home/thoresen/Documents/PhD/plots/nua/Datasets/LHC17i2f/ESD/maxstrip0/fmdreco_tpcreco_runsz.root";
  TString nua_file = "/home/thoresen/Documents/PhD/plots/nua/Datasets/LHC17i2f/ESD/fmdreco_tpcreco_runsz.root";

  Bool_t doNUA = kFALSE;

  task->fSettings.fileName = resName;
  mgr->AddTask(task);

  AliAnalysisDataContainer *coutput_recon =
  mgr->CreateContainer(resName,
   AliForwardFlowResultStorage::Class(), //TList::Class(),
   AliAnalysisManager::kOutputContainer,
   mgr->GetCommonFileName());
  //task->fSettings.fDataType = task->fSettings.kRECON;
  mgr->ConnectOutput(task, 1, coutput_recon);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  AliAnalysisDataContainer* valid = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("event_selection_xchange");
  task->ConnectInput(1,valid);


  if (doNUA){
    std::cout << "doing NUA" << std::endl;
    TObjArray* taskContainers = mgr->GetContainers();
    AliAnalysisDataContainer* weights;
    
    TObjArray *tx = nua_file.Tokenize("/");
    TObjArray *ty = ((TObjString *)(tx->At(tx->GetEntries()-1)))->String().Tokenize(".");
    TString nuaobject =  ((TObjString *)(ty->At(0)))->String();
    std::cout << "I-AddTaskForwardSecondaries: NUA container " << nuaobject << std::endl;

    weights = (AliAnalysisDataContainer*) taskContainers->FindObject(nuaobject);
    
    if (!weights) {
      std::cout << "I-AddTaskForwardSecondaries: " << nuaobject << " weights not defined - reading now. " << std::endl;
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
