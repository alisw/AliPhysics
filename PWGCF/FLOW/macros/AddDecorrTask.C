AliAnalysisDecorrTask* AddDecorrTask(TString name = "name", Bool_t IsMC = kFALSE, Bool_t use3DWeights = kTRUE, Bool_t useOwnWeights = kFALSE, TString s2DWeightsFile = "", TString s3DWeightsFile = "", const char* suffix = "")
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
    fileName += Form(":%s",suffix);      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisDecorrTask* task = new AliAnalysisDecorrTask(name.Data(),IsMC);   
    if(!task) return 0x0;

    task->SetUseWeights3D(use3DWeights);
    task->SetUseOwnWeights(useOwnWeights);
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("FlowList_%s",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("FlowWeights_%s",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,3,mgr->CreateContainer(Form("QAlist_%s",suffix),TList::Class(),AliAnalysisManager::kOutputContainer,fileName.Data()));

    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid

    if(useOwnWeights)
    {
      if(!use3DWeights) 
      {
          TObjArray* taskContainers = mgr->GetContainers();
          if(!taskContainers) { printf("E-AddTaskUniFlow: Task containers does not exists!\n"); return NULL; }

          // check if the input weights are already loaded (e.g. in different subwagon)
          AliAnalysisDataContainer* weights2D = (AliAnalysisDataContainer*) taskContainers->FindObject("inputWeights");
          if(!weights2D) 
          {  
              // if it does not exists create it
              // in case of non-local run, establish connection to ALiEn for loading the weights
              if(s2DWeightsFile.Contains("alien://")) { gGrid->Connect("alien://"); }

              TFile* weights2D_file = TFile::Open(s2DWeightsFile.Data(),"READ");
              if(!weights2D_file) { printf("E-AddTaskUniFlow: Input file with weights not found!\n"); return NULL; }

              TList* weights2D_list = (TList*) weights2D_file->Get("weights");
              if(!weights2D_list) { printf("E-AddTaskUniFlow: Input list with weights not found!\n"); weights2D_file->ls(); return NULL; }

              AliAnalysisDataContainer* cInputWeights2D = mgr->CreateContainer("inputWeights",TList::Class(), AliAnalysisManager::kInputContainer);
              cInputWeights2D->SetData(weights2D_list);
              mgr->ConnectInput(task,1,cInputWeights2D);
          }
          else 
          {
              // connect existing container
              mgr->ConnectInput(task,1,weights2D);
          }
      }
      else
      {
          TObjArray* taskContainers3D = mgr->GetContainers();
          if(!taskContainers3D) { printf("E-AddTaskUniFlow: Task containers does not exists!\n"); return NULL; }

          // check if the input weights are already loaded (e.g. in different subwagon)
          AliAnalysisDataContainer* weightsOwn3D = (AliAnalysisDataContainer*) taskContainers3D->FindObject("inputWeights");
          if(!weightsOwn3D) 
          {  
              // if it does not exists create it

              // in case of non-local run, establish connection to ALiEn for loading the weights
              if(s3DWeightsFile.Contains("alien://")) { gGrid->Connect("alien://"); }

              TFile* weightsOwn3D_file = TFile::Open(s3DWeightsFile.Data(),"READ");
              if(!weightsOwn3D_file) { printf("E-AddTaskUniFlow: Input file with weights not found at all!\n"); return NULL; }

              TList* weightsOwn3D_list = (TList*) weightsOwn3D_file->Get("weights");
              if(!weightsOwn3D_list) { printf("E-AddTaskUniFlow: Input list with weights not found!\n"); weightsOwn3D_file->ls(); return NULL; }

              AliAnalysisDataContainer* cInputWeightsOwn3D = mgr->CreateContainer("inputWeights",TList::Class(), AliAnalysisManager::kInputContainer);
              cInputWeightsOwn3D->SetData(weightsOwn3D_list);
              mgr->ConnectInput(task,1,cInputWeightsOwn3D);
          }
          else 
          {
              // connect existing container
              mgr->ConnectInput(task,1,weightsOwn3D);
          }       
      } 
    }
    else 
    {
      TObjArray* taskContainersVy = mgr->GetContainers();
      if(!taskContainersVy) { printf("E-AddTaskUniFlow: Task containers does not exists!\n"); return NULL; }

      // check if the input weights are already loaded (e.g. in different subwagon)
      AliAnalysisDataContainer* weights3D = (AliAnalysisDataContainer*) taskContainersVy->FindObject("inputWeights");
      if(!weights3D) 
      {  
        // in case of non-local run, establish connection to ALiEn for loading the weights
        if(s3DWeightsFile.Contains("alien://")) { gGrid->Connect("alien://"); }

        TFile* weights3D_file = TFile::Open(s3DWeightsFile.Data(),"READ");
        if(!weights3D_file) { printf("E-AddTaskUniFlow: Input file with weights not found!\n"); return NULL; }

        TList* weights3D_list = (TList*) weights3D_file->Get("WeightList");
        if(!weights3D_list) { printf("E-AddTaskUniFlow: Input list with weights not found!\n"); weights3D_file->ls(); return NULL; }

        AliAnalysisDataContainer* cInputWeights3D = mgr->CreateContainer("inputWeights",TList::Class(), AliAnalysisManager::kInputContainer);
        cInputWeights3D->SetData(weights3D_list);
        mgr->ConnectInput(task,1,cInputWeights3D);
      }
      else 
      {
        // connect existing container
        mgr->ConnectInput(task,1,weights3D);
      }
    }
  return task;
}




