AliAnalysisDecorrTask* AddDecorrTask(TString name = "name", Bool_t IsMC = kFALSE, bool IsOnTheFly = kFALSE, TString sWeightsFile = "", const char* suffix = "")
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
    // now we create an instance of your task
    AliAnalysisDecorrTask* task = new AliAnalysisDecorrTask(name.Data(),IsMC,IsOnTheFly);   
    if(!task) return 0x0;
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("%s",suffix), TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root"));
    mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("QA_%s",suffix),TList::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root"));
    TObjArray* taskContainers = mgr->GetContainers();
    if(!taskContainers) { printf("E-AddDecorrTask: Task containers does not exists!\n"); return NULL; }
    // check if the input weights are already loaded (e.g. in different subwagon)
    if(!IsMC&&!IsOnTheFly)
    {
      AliAnalysisDataContainer* weights = (AliAnalysisDataContainer*) taskContainers->FindObject("inputWeights");
      if(!weights) 
      {  
        // in case of non-local run, establish connection to ALiEn for loading the weights
        if(sWeightsFile.Contains("alien://")) { gGrid->Connect("alien://"); }

        TFile* weights_file = TFile::Open(sWeightsFile.Data(),"READ");
        if(!weights_file) { printf("E-AddDecorrTask: Input file with weights not found!\n"); return NULL; }

        TList* weights_list = (TList*) weights_file->Get("WeightList");
        if(!weights_list) { printf("E-AddDecorrTask: Input list with weights not found!\n"); weights_file->ls(); return NULL; }

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
 
    
  return task;
}




