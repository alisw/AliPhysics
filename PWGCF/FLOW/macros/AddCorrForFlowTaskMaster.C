AliAnalysisTaskCorrForFlowMaster* AddCorrForFlowTaskMaster(TString name = "name", TString efficiencyFile = "",TString calibrationFile = "",const char* suffix = "")
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

    Bool_t useCalibration = kFALSE;
    if(!calibrationFile.IsNull()) useCalibration = kTRUE;
    
    AliAnalysisTaskCorrForFlowMaster* task = new AliAnalysisTaskCorrForFlowMaster(taskName.Data(), useEfficiency, useCalibration);
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

        TList* efficiency_list = (TList*) efficiency_file->Get("Efficiency2D_wFD");
        if(!efficiency_list) { printf("E-AddTask: Input list with efficiency not found!\n"); efficiency_file->ls(); return NULL; }

        AliAnalysisDataContainer* cInputEfficiency = mgr->CreateContainer("inputEfficiency",TList::Class(), AliAnalysisManager::kInputContainer);
        cInputEfficiency->SetData(efficiency_list);
        mgr->ConnectInput(task,1,cInputEfficiency);
      }
      else {
        mgr->ConnectInput(task,1,efficiency);
      }
    }

    if(useCalibration) {
      TObjArray* taskContainers = mgr->GetContainers();
      if(!taskContainers) { printf("Task containers does not exists!\n"); return NULL; }

      AliAnalysisDataContainer* calib = (AliAnalysisDataContainer*) taskContainers->FindObject("calib");
      if(!calib) {
        if(calibrationFile.Contains("alien://")) { gGrid->Connect("alien://"); }

        printf("input file name: %s \n", calibrationFile.Data());

        TFile* calib_file = TFile::Open(calibrationFile.Data(),"READ");
        if(!calib_file) { printf("Input file with AMPT centrality calibration not found!\n"); return NULL; }

        TH1D* calib_histo = (TH1D*)calib_file->Get("hcent");
        if(!calib_histo) { printf("E-AddTask: Input list with AMPT centrality calibration not found!\n"); calib_file->ls(); return NULL; }

        AliAnalysisDataContainer* cInputCalib = mgr->CreateContainer("calib",TH1D::Class(), AliAnalysisManager::kInputContainer);
        cInputCalib->SetData(calib_histo);
        if(useEfficiency) mgr->ConnectInput(task,2,cInputCalib);
        else mgr->ConnectInput(task,1,cInputCalib);
      }
      else {
        if(useEfficiency) mgr->ConnectInput(task,2,calib);
        else mgr->ConnectInput(task,1,calib);
      }
    }

    return task;
}



