AliAnalysisDataContainer* makeWeightContainer(TString nua_file, TString containerName)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer* weights;
  if (nua_file.Contains("alien:")) TGrid::Connect("alien:");
  TFile* file;
  file = TFile::Open(nua_file.Data(), "READ");

  if(!file) { printf("E-MyAddTask: Input file with differential weights not found!\n"); return NULL; }

  TList* weights_list = new TList();
  weights_list->SetName("nuaWeights");
  
  TH2F* nuacentral = new TH2F();

  file->GetObject("fHistPhiEta", nuacentral);
  nuacentral->SetDirectory(0);
  nuacentral->SetNameTitle("nuacentral","nuacentral");

  file->Close();

  weights_list->Add(nuacentral);

  weights = mgr->CreateContainer(containerName,TList::Class(), AliAnalysisManager::kInputContainer,Form("%s", mgr->GetCommonFileName()));
  weights->SetData(weights_list);
  return weights;
}


void connectContainer(AliAnalysisDataContainer* container,AliAnalysisDecorrTask* task)
{

  task->nuacentral = static_cast<TH2F*>( static_cast<TList*>(container->GetData())->FindObject("nuacentral") );
  task->nuacentral->SetDirectory(0);
}

AliAnalysisDecorrTask* AddDecorrTask(TString name = "name", TString dirname = "", TString sWeightsFile = "", TString sVytauWeights = "")
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
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += Form(":%s",dirname.Data());      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisDecorrTask* task = new AliAnalysisDecorrTask(name.Data());   
    if(!task) return 0x0;

    Bool_t useWeights3D = task->GetUseWeights3D();
    
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("FlowList_%s",dirname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid

    if(!useWeights3D) 
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
    else 
    {
      TObjArray* taskContainersVy = mgr->GetContainers();
      if(!taskContainersVy) { printf("E-AddTaskUniFlow: Task containers does not exists!\n"); return NULL; }

      // check if the input weights are already loaded (e.g. in different subwagon)
      AliAnalysisDataContainer* weightsVy = (AliAnalysisDataContainer*) taskContainersVy->FindObject("inputWeights");
      if(!weightsVy) 
      {  
        // in case of non-local run, establish connection to ALiEn for loading the weights
        if(sVytauWeights.Contains("alien://")) { gGrid->Connect("alien://"); }

        TFile* weights_fileVy = TFile::Open(sVytauWeights.Data(),"READ");
        if(!weights_fileVy) { printf("E-AddTaskUniFlow: Input file with weights not found!\n"); return NULL; }

        TList* weights_listVy = (TList*) weights_fileVy->Get("WeightList");
        if(!weights_listVy) { printf("E-AddTaskUniFlow: Input list with weights not found!\n"); weights_fileVy->ls(); return NULL; }

        AliAnalysisDataContainer* cInputWeightsVy = mgr->CreateContainer("inputWeights",TList::Class(), AliAnalysisManager::kInputContainer);
        cInputWeightsVy->SetData(weights_listVy);
        mgr->ConnectInput(task,1,cInputWeightsVy);
      }
      else 
      {
        // connect existing container
        mgr->ConnectInput(task,1,weightsVy);
      }
    }
  return task;
}
