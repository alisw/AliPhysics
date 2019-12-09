class AliAnalysisDataContainer;
class AliAnalysisTaskESEFlow;

AliAnalysisTaskESEFlow* AddESEFlowTask(TString name = "name",TString dirname ="MyTask", TString sWeightsFile = "", TString sVWeights = "", TString sqSelCuts = "")
{
  // get the manager via the static access member. since it's static, you don't need
  // to create an instance of the class here to call the function
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
      return 0x0;
  }
  // get the input event handler, again via a static method. 
  // this handler is part of the managing system and feeds events
  // to your task
  if (!mgr->GetInputEventHandler()) {
      return 0x0;
  }

  Bool_t bUseOwnWeights = kFALSE;
  Bool_t bUseqSelCuts = kTRUE;


  // by default, a file is open for writing. here, we get the filename
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += Form(":Task%s", dirname.Data());      // create a subfolder in the file
  // now we create an instance of your task
  AliAnalysisTaskESEFlow* task = new AliAnalysisTaskESEFlow(name.Data());   
  if(!task) return 0x0;
  task->SelectCollisionCandidates(AliVEvent::kAnyINT);
  // add your task to the manager
  mgr->AddTask(task);
  // your task needs input: here we connect the manager to your task
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  // same for the output
  mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("MyOutputContainer%s",dirname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("Observables%s",dirname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,3,mgr->CreateContainer(Form("c_n{n} distributions%s",dirname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,4,mgr->CreateContainer(Form("d_n{n} distributions%s",dirname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,5,mgr->CreateContainer(Form("q_n distributions%s",dirname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,6,mgr->CreateContainer(Form("d_n{n} dist after q selection%s",dirname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,7,mgr->CreateContainer(Form("c_n{n} dist after q selection%s",dirname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,8,mgr->CreateContainer(Form("fQAEvents%s",dirname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  // in the end, this macro returns a pointer to your task. this will be convenient later on
  // when you will run your analysis in an analysis train on grid

  // RUN BY RUN
  if(bUseOwnWeights) 
  {
    TObjArray* taskContainers = mgr->GetContainers();
    if(!taskContainers) { printf("E-AddTaskESEFlow: Task containers does not exists!\n"); return NULL; }

    // check if the input weights are already loaded (e.g. in different subwagon)
    AliAnalysisDataContainer* weights = (AliAnalysisDataContainer*) taskContainers->FindObject("inputWeights");
    if(!weights) 
    {  
      // if it does not exists create it

      // in case of non-local run, establish connection to ALiEn for loading the weights
      if(sWeightsFile.Contains("alien://")) { gGrid->Connect("alien://"); }

      TFile* weights_file = TFile::Open(sWeightsFile.Data(),"READ");
      if(!weights_file) { printf("E-AddTaskESEFlow: Input file with weights not found!\n"); return NULL; }

      TList* weights_list = (TList*) weights_file->Get("weights");
      if(!weights_list) { printf("E-AddTaskESEFlow: Input list with weights not found!\n"); weights_file->ls(); return NULL; }

      AliAnalysisDataContainer* cInputWeights = mgr->CreateContainer("inputWeights",TList::Class(), AliAnalysisManager::kInputContainer);
      cInputWeights->SetData(weights_list);
      mgr->ConnectInput(task,1,cInputWeights);
    }
    else 
    {
      // connect existing container
      mgr->ConnectInput(task,1,weights);
    }
  }
  if(!bUseOwnWeights)
  {
    TObjArray* taskContainersVy = mgr->GetContainers();
    if(!taskContainersVy) { printf("E-AddTaskESEFlow: Task containers does not exists!\n"); return NULL; }

    // check if the input weights are already loaded (e.g. in different subwagon)
    AliAnalysisDataContainer* weightsVy = (AliAnalysisDataContainer*) taskContainersVy->FindObject("inputWeights");
    if(!weightsVy) {  
      // if it does not exists create it
      // in case of non-local run, establish connection to ALiEn for loading the weights
      if(sVWeights.Contains("alien://")) { gGrid->Connect("alien://"); }

      TFile* weights_fileVy = TFile::Open(sVWeights.Data(),"READ");
      if(!weights_fileVy) { printf("E-AddTaskESEFlow: Input file with weights not found!\n"); return NULL; }

      TList* weights_listVy = (TList*) weights_fileVy->Get("WeightList");
      if(!weights_listVy) { printf("E-AddTaskESEFlow: Input list with weights not found!\n"); weights_fileVy->ls(); return NULL; }

      AliAnalysisDataContainer* cInputWeightsVy = mgr->CreateContainer("inputWeights",TList::Class(), AliAnalysisManager::kInputContainer);
      cInputWeightsVy->SetData(weights_listVy);
      mgr->ConnectInput(task,1,cInputWeightsVy);
    }
    else {
      // connect existing container
      mgr->ConnectInput(task,1,weightsVy);
    }
  }

  if(bUseqSelCuts){
    TObjArray* taskContainersCuts = mgr->GetContainers();
    if(!taskContainersCuts) { printf("E-AddTaskESEFlow: Task containers does not exists!\n"); return NULL; }

    AliAnalysisDataContainer* qcuts = (AliAnalysisDataContainer*) taskContainersCuts->FindObject("inputCuts");
    if(!qcuts) 
    { 
      if(sWeightsFile.Contains("alien://")) { gGrid->Connect("alien://"); }

      TFile* qcuts_file = TFile::Open(sqSelCuts.Data(),"READ");
      if(!qcuts_file) { printf("E-AddTaskESEFlow: Input file with q-cuts not found!\n"); return NULL; }

      TList* qcuts_list = (TList*) qcuts_file->Get("qCuts");
      if(!qcuts_list) { printf("E-AddTaskESEFlow: Input list with weights not found!\n"); qcuts_file->ls(); return NULL; }

      AliAnalysisDataContainer* cInputCuts = mgr->CreateContainer("inputCuts",TList::Class(), AliAnalysisManager::kInputContainer);
      cInputCuts->SetData(qcuts_list);
      mgr->ConnectInput(task,2,cInputCuts);
    }
    else 
    {
      mgr->ConnectInput(task,2,qcuts);
    }
  }
    
  return task;
}