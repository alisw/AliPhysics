AliAnalysisTaskFlowExample* AddFlowExampleTask(TString name = "FlowTask", TString weightsFile = "", const char* suffix = "")
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":FlowExampleTask";
    TString taskName = Form("%s_%s",name.Data(),suffix);

    Bool_t useWeights = kFALSE;
    if(!weightsFile.IsNull()) useWeights = kTRUE;

    AliAnalysisTaskFlowExample* task = new AliAnalysisTaskFlowExample(taskName.Data(), useWeights);
    if(!task) return 0x0;

    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer("Output", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    if(useWeights) {
      TObjArray* taskContainers = mgr->GetContainers();
      if(!taskContainers) { printf("Task containers does not exists!\n"); return NULL; }

      AliAnalysisDataContainer* weights = (AliAnalysisDataContainer*) taskContainers->FindObject("inputWeights");
      if(!weights) {
        if(weightsFile.Contains("alien://")) { gGrid->Connect("alien://"); }

        printf("input file name: %s \n", weightsFile.Data());

        TFile* weights_file = TFile::Open(weightsFile.Data(),"READ");
        if(!weights_file) { printf("Input file with weights not found!\n"); return NULL; }

        TList* weights_list = (TList*) weights_file->Get("weightsList");
        if(!weights_list) { printf("E-AddTaskUniFlow: Input list with weights not found!\n"); weights_file->ls(); return NULL; }

        AliAnalysisDataContainer* cInputWeights = mgr->CreateContainer("inputWeights",TList::Class(), AliAnalysisManager::kInputContainer);
        cInputWeights->SetData(weights_list);
        mgr->ConnectInput(task,1,cInputWeights);
      }
      else {
        mgr->ConnectInput(task,1,weights);
      }
    }

    return task;
}
