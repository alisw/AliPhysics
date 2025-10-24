void AddTask_EtaPrimeMCStudies( TString specialSetting  = "EtaPrimeBiased",
                                Double_t etaPrimeSup    = 0.12){

    // Get AnalysisManager
    AliAnalysisManager  *mgr    = AliAnalysisManager::GetAnalysisManager();
    if(!mgr){
        Error("AddTask_EtaPrimeMCStudies","No analysis manager found.");
        return;
    }
    std::cout << "AliAnalysisManager found " << mgr << std::endl;

    // Get InputEventHandler and input container
    AliVEventHandler            *inputHandler   = mgr->GetInputEventHandler();
    AliAnalysisDataContainer    *cinput         = mgr->GetCommonInputContainer();

    // Add task to the Analysis Manager
    AliAnalysisTaskEtaPrimeMCStudies    *task   = NULL;
    task    = new AliAnalysisTaskEtaPrimeMCStudies("EtaPrimeMCStudies");
    if( !task ){
        Error("AddTask_EtaPrimeMCStudies","Failed to create the task.");
        return;
    }

    std::cout << "Task created? " << task << endl;

    // Connect containers
    AliAnalysisDataContainer    *coutput    = mgr->CreateContainer( "EtaPrimeMCStudies", TList::Class(), AliAnalysisManager::kOutputContainer, 
                                                                    specialSetting.Contains("NotSuppressed") ? Form("EtaPrimeMCStudies_%s_%.2f.root", specialSetting.Data(), etaPrimeSup) : Form("EtaPrimeMCStudies_%s.root", specialSetting.Data()) );
    AliAnalysisDataContainer    *couttree   = mgr->CreateContainer( "EtaPrimeMCStudiesTree", TTree::Class(), AliAnalysisManager::kOutputContainer, 
                                                                    specialSetting.Contains("NotSuppressed") ? Form("EtaPrimeMCStudiesTree_%s_%.2f.root", specialSetting.Data(), etaPrimeSup) : Form("EtaPrimeMCStudiesTree_%s.root", specialSetting.Data()) );

    mgr->AddTask(task);
    mgr->ConnectInput(task,0,cinput);
    mgr->ConnectOutput(task,1,coutput);
    mgr->ConnectOutput(task,2,couttree);

    return;    
}