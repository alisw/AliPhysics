AliAnalysisTaskCorrForFlow* AddCorrForFlowTask(TString name = "name", TString efficiencyFile = "",const char* suffix = "")
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":CorrForFlow";
    TString taskName = Form("%s_%s",name.Data(),suffix);

    Bool_t useEfficiency = kFALSE;
    if(!efficiencyFile.IsNull()) useEfficiency = kTRUE;

    AliAnalysisTaskCorrForFlow* task = new AliAnalysisTaskCorrForFlow(taskName.Data(), useEfficiency);
    if(!task) return 0x0;
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("Charged_%s",taskName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    if(useEfficiency) {
      TObjArray* taskContainers = mgr->GetContainers();
      if(!taskContainers) { printf("Task containers does not exists!\n"); return NULL; }

      AliAnalysisDataContainer* efficiency = (AliAnalysisDataContainer*) taskContainers->FindObject("inputEfficiency");
      if(!efficiency) {
        if(efficiencyFile.Contains("alien://")) { gGrid->Connect("alien://"); }

        printf("input file name: %s \n", efficiencyFile.Data());

        TFile* efficiency_file = TFile::Open(efficiencyFile.Data(),"READ");
        if(!efficiency_file) { printf("Input file with efficiency not found!\n"); return NULL; }

        TList* efficiency_list = (TList*) efficiency_file->Get("EffAndFD");
        if(!efficiency_list) { printf("E-AddTaskUniFlow: Input list with efficiency not found!\n"); efficiency_file->ls(); return NULL; }

        AliAnalysisDataContainer* cInputEfficiency = mgr->CreateContainer("inputEfficiency",TList::Class(), AliAnalysisManager::kInputContainer);
        cInputEfficiency->SetData(efficiency_list);
        mgr->ConnectInput(task,1,cInputEfficiency);
      }
      else {
        mgr->ConnectInput(task,1,efficiency);
      }
    }

    return task;
}
