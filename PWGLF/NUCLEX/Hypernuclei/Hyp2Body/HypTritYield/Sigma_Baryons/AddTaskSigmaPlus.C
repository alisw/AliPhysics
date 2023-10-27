////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Sigma Plus Class                                                      //
//                                                                        //
//  This class is used for the reconstruction of Sigma+ baryons in the    //
//  Sigma+ -> p+ + pi0 decay channel.                                     //
//                                                                        //
//  Author: B.Heybeck (b.heybeck@cern.ch)                                 //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

AliAnalysisTaskSigmaPlus* AddTaskSigmaPlus(TString name = "SigmaPlus")
{
    // Get the analysis manager via the static access member. Since it is static there is no need
    // to create an instance of the class here to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    // Get the input event handler, via a static method. 
    // This handler is part of the managing system and feeds events to the task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    // By default, a file is opened for writing. Get the filename:
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":SigmaBaryons";  // create a subfolder in the file
    // Create an instance of the task
    AliAnalysisTaskSigmaPlus* task = new AliAnalysisTaskSigmaPlus(name.Data());   
    if(!task) return 0x0;
    // Trigger settings
        // Select Trigger in Wagon Configuration!
    // Add the task to the manager
    mgr->AddTask(task);
    // Connect the input data to the task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // Connect the output containers to the task
    mgr->ConnectOutput(task,1 ,mgr->CreateContainer("Histogram_List", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2 ,mgr->CreateContainer("Sigma_Cand_Tree", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,3 ,mgr->CreateContainer("Sigma_Pair_Tree_SE", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,4 ,mgr->CreateContainer("Sigma_Pair_Tree_ME", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,5 ,mgr->CreateContainer("Sigma_ME_Background_Tree", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,6 ,mgr->CreateContainer("Sigma_PHOS_Cand_Tree", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,7 ,mgr->CreateContainer("Sigma_Pair_Tree_PHOS_SE", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,8 ,mgr->CreateContainer("Sigma_Pair_Tree_PHOS_ME", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,9 ,mgr->CreateContainer("Sigma_PHOS_ME_Bkg_Tree", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,10,mgr->CreateContainer("Sigma_Calc_Cand_Tree", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,11,mgr->CreateContainer("Proton_Tree", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // Finally, a pointer to the created task is returned
    return task;
}
