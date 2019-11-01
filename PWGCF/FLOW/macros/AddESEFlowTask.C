AliAnalysisDataContainer* makeWeightContainer(TString nua_file, TString containerName){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer* weights;
  if (nua_file.Contains("alien:")) TGrid::Connect("alien:");
  TFile* file;
  file = TFile::Open(nua_file.Data(), "READ");

  if(!file) { printf("E-MyAddTask: Input file with differential weights not found!\n"); return NULL; }

  TList* weights_list_newtemp = new TList();
  weights_list_newtemp->SetName("nuaWeights");
 
  TH2F* nuacentral = new TH2F();

  file->GetObject("PhiEtaWeights", nuacentral);
  nuacentral->SetDirectory(0);
  nuacentral->SetNameTitle("nuacentral","nuacentral");

  file->Close();

  weights_list_newtemp->Add(nuacentral);

  weights = mgr->CreateContainer(containerName,TList::Class(), AliAnalysisManager::kInputContainer,Form("%s", mgr->GetCommonFileName()));
  weights->SetData(weights_list_newtemp);
  return weights;
}

void connectContainer(AliAnalysisDataContainer* container,AliAnalysisTaskESEFlow* task)
{
  task->nuacentral = static_cast<TH2F*>( static_cast<TList*>(container->GetData())->FindObject("nuacentral") );
  task->nuacentral->SetDirectory(0);
}

AliAnalysisTaskESEFlow* AddESEFlowTask(TString name = "name",TString dirname ="MyTask",TString nua_file = "", TString sWeightsFile = "", TString sVWeights = "", TString sqSelCuts = "")
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

    Bool_t bUseMyWeights = kFALSE; //dont ever change this to kTRUE

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
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid

  

  if(bUseMyWeights)
  {
    // read new histogram
    TObjArray* taskContainers = mgr->GetContainers();
    AliAnalysisDataContainer* weights_newtemp;
   
    TObjArray *tx = nua_file.Tokenize("/");
    TObjArray *ty = ((TObjString *)(tx->At(tx->GetEntries()-1)))->String().Tokenize(".");
    TString nuaobject =  ((TObjString *)(ty->At(0)))->String();
    std::cout << nuaobject << std::endl;

    weights_newtemp = (AliAnalysisDataContainer*) taskContainers->FindObject(nuaobject);

    if (!weights_newtemp) {
      std::cout << "I-AddTask: " << nuaobject << " weights not defined - reading now. " << std::endl;
      weights_newtemp = makeWeightContainer(nua_file,nuaobject);
    }
    connectContainer( weights_newtemp, task);
  }


  // RUN BY RUN
  if(bUseOwnWeights) 
  {
    TObjArray* taskContainers = mgr->GetContainers();
    if(!taskContainers) { printf("E-AddTaskUniFlow: Task containers does not exists!\n"); return NULL; }

    // check if the input weights are already loaded (e.g. in different subwagon)
    AliAnalysisDataContainer* weights = (AliAnalysisDataContainer*) taskContainers->FindObject("inputWeights");
    if(!weights) 
    {  
      // if it does not exists create it

      // in case of non-local run, establish connection to ALiEn for loading the weights
      if(sWeightsFile.Contains("alien://")) { gGrid->Connect("alien://"); }

      TFile* weights_file = TFile::Open(sWeightsFile.Data(),"READ");
      if(!weights_file) { printf("E-AddTaskUniFlow: Input file with weights not found!\n"); return NULL; }

      TList* weights_list = (TList*) weights_file->Get("weights");
      if(!weights_list) { printf("E-AddTaskUniFlow: Input list with weights not found!\n"); weights_file->ls(); return NULL; }

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
    if(!taskContainersVy) { printf("E-AddTaskUniFlow: Task containers does not exists!\n"); return NULL; }

    // check if the input weights are already loaded (e.g. in different subwagon)
    AliAnalysisDataContainer* weightsVy = (AliAnalysisDataContainer*) taskContainersVy->FindObject("inputWeights");
    if(!weightsVy) {  
      // if it does not exists create it
      // in case of non-local run, establish connection to ALiEn for loading the weights
      if(sVWeights.Contains("alien://")) { gGrid->Connect("alien://"); }

      TFile* weights_fileVy = TFile::Open(sVWeights.Data(),"READ");
      if(!weights_fileVy) { printf("E-AddTaskUniFlow: Input file with weights not found!\n"); return NULL; }

      TList* weights_listVy = (TList*) weights_fileVy->Get("WeightList");
      if(!weights_listVy) { printf("E-AddTaskUniFlow: Input list with weights not found!\n"); weights_fileVy->ls(); return NULL; }

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
    TObjArray* taskContainersqSel = mgr->GetContainers();
    if(!taskContainersqSel) { printf("Task containers for q Selection does not exist!\n"); return NULL; }

    AliAnalysisDataContainer* qSelcuts = (AliAnalysisDataContainer*) taskContainersqSel->FindObject("inputqCuts");

    if(!qSelcuts){
      if(sqSelCuts.Contains("alien://")) { gGrid->Connect("alien://"); }

      TFile* qCuts_file = TFile::Open(sqSelCuts.Data(),"READ");
      if(!qCuts_file) { printf("Input file with q selections cuts not found! \n"); return NULL; }

      TTree* qCuts_Tree = (TTree*) qCuts_file->Get("tree");
      if(!qCuts_Tree) { printf("Input tree with q selection cuts not found! \n"); qCuts_file->ls(); return NULL; }

      AliAnalysisDataContainer* cinputqCuts = mgr->CreateContainer("inputqCuts",TTree::Class(), AliAnalysisManager::kInputContainer);
      cinputqCuts->SetData(qCuts_Tree);
      mgr->ConnectInput(task,2,cinputqCuts);
    }
    else{ 
      mgr->ConnectInput(task,2,qSelcuts);
    }
  }
    
  return task;
}
