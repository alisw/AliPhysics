AliAnalysisTaskSEDmesonPIDSysProp *AddTaskDmesonPIDSysProp(int ch = AliAnalysisTaskSEDmesonPIDSysProp::kDstoKKpi,
                                                           int pid = AliAnalysisTaskSEDmesonPIDSysProp::kStrongPID,
                                                           TString cutfilename = "",
                                                           TString cutobjname = "AnalysisCuts",
                                                           TString PIDsystfilename = "",
                                                           TString postname = "")
{
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTask", "No analysis manager to connect to.");
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskImpParDistrib", "This task requires an input event handler");
        return NULL;
    }
    
    
    Bool_t stdcuts=kFALSE;
    TFile* filecuts = NULL;
    if( cutfilename.EqualTo("") ) {
        stdcuts=kTRUE;
    } 
    else {
        filecuts=TFile::Open(cutfilename.Data());
        if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
            ::Fatal("AddTaskSingleTrackPIDSysPropagation", "Cut object not found: analysis will not start!\n");
        }
        else printf("Cut object correctly found\n");
    }
    
    //Analysis Task
    AliRDHFCuts* analysiscuts = (AliRDHFCuts*)filecuts->Get(cutobjname.Data());    
    if(analysiscuts) Printf("Cut file found!");    
    
    AliAnalysisTaskSEDmesonPIDSysProp *Task = new AliAnalysisTaskSEDmesonPIDSysProp(ch, analysiscuts);
    Task->SetDebugLevel(1);
    Task->SetPIDStrategy(pid);
    bool loadfile = Task->LoadEffSystFile(PIDsystfilename);
    if(!loadfile) {
        ::Fatal("AddTaskSingleTrackPIDSysPropagation", "Histos with single-track systematic not found: analysis will not start!\n");
    }
    mgr->AddTask(Task);

    // Create containers for input/output
    TString mesonname = "";
    if(ch == AliAnalysisTaskSEDmesonPIDSysProp::kD0toKpi)               mesonname = "D0";
    else if(ch == AliAnalysisTaskSEDmesonPIDSysProp::kDplustoKpipi)     mesonname = "Dplus";
    else if(ch == AliAnalysisTaskSEDmesonPIDSysProp::kDstartoKpipi)     mesonname = "Dstar";
    else if(ch == AliAnalysisTaskSEDmesonPIDSysProp::kDstoKKpi)         mesonname = "Ds";
    TString PIDname = "";
    if(pid == AliAnalysisTaskSEDmesonPIDSysProp::kStrongPID)            PIDname = "strongPID";
    else if(pid == AliAnalysisTaskSEDmesonPIDSysProp::kConservativePID) PIDname = "conservativePID";
    else if(pid == AliAnalysisTaskSEDmesonPIDSysProp::knSigmaPID)       PIDname = "nSigmaPID";

    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    outputFileName += Form(":PWGHF_D2H_PIDeffsyst_%s_%s_%s",mesonname.Data(),PIDname.Data(),postname.Data());
    
    AliAnalysisDataContainer *coutput =0x0;
    coutput = mgr->CreateContainer(Form("systUnc_%s_%s_%s",mesonname.Data(),PIDname.Data(),postname.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName );
    
    mgr->ConnectInput(Task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(Task,1,coutput);
    
    return Task;
}
