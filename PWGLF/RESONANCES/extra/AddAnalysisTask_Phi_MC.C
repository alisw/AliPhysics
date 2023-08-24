// TODO LIST
// TODO: Set fName as last par per Train set-up indications
// TODO: Set Trigger Mask to be set from the AddAnalysisTask argument list


AliAnalysisTask_Phi_MC* AddAnalysisTask_Phi_MC( TString fName = "PhiCount_STD", Long_t kParticlePDG = 333 )
{
    // Analysis Manager
    AliAnalysisManager         *fAliAnlManager      =   AliAnalysisManager::GetAnalysisManager();
    
    // Filename
    TString fileName                                =   AliAnalysisManager::GetCommonFileName();
    
    // Task initialisation
    AliAnalysisTask_Phi_MC    *fAliAnlTask          =   new AliAnalysisTask_Phi_MC(fName.Data());
    
    // Checks
    if (!fAliAnlManager)                            return 0x0;
    if (!fAliAnlManager->GetInputEventHandler())    return 0x0;
    if (!fAliAnlTask)                               return 0x0;

    //  Task settings
    fAliAnlTask ->  SetParticlePDG          ( kParticlePDG );
    fAliAnlTask ->  SetSPCompute            ( kFALSE );   // Not implemented
    fAliAnlTask ->  SetSPWighted            ( kFALSE );   // Not implemented
    fAliAnlTask ->  SetRTCompute            ( kFALSE );   // Not implemented
    
    // Input
    fAliAnlManager->ConnectInput(fAliAnlTask,0,fAliAnlManager->GetCommonInputContainer());
    
    // Output
    // - // TLists
    fAliAnlManager->ConnectOutput(fAliAnlTask,1,fAliAnlManager->CreateContainer(Form("fAnalysisOutputList_%s",fName.Data()),                    TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    fAliAnlManager->ConnectOutput(fAliAnlTask,2,fAliAnlManager->CreateContainer(Form("fQCOutputList_%s",fName.Data()),                          TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
    // - // TTrees
    fAliAnlManager->ConnectOutput(fAliAnlTask,3,fAliAnlManager->CreateContainer(Form("PhiCandidate_%s",fName.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    
    // Add task
    fAliAnlManager->AddTask(fAliAnlTask);
    
    return fAliAnlTask;
}
