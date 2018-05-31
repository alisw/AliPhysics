///////////////////////////////////////////////////////////////////
//                                                               //            
// AddMyTask                                                     //
// Author: Erin Gauger                                           //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTask* AddTaskTPCCalBeauty(Double_t centMin=0, Double_t centMax=10, Bool_t applySSCut = kTRUE, TString ContNameExt = " ")
{
    // get the manager via the static access member
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    
    // get the input event handler, again via a static method.
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    
    ////////////
    //  INT7  //
    ////////////

    //FOR EMCAL CLUSTERS.......................................................................
    
    AliAnalysisTaskTPCCalBeauty* taskBFEemc = new AliAnalysisTaskTPCCalBeauty("bfeemc");
    if(!taskBFEemc) return 0x0;
    // add your task to the manager
    mgr->AddTask(taskBFEemc);
    taskBFEemc->SetSSCut(applySSCut);
    taskBFEemc->SetClusterTypeEMC(kTRUE);
    taskBFEemc->SetClusterTypeDCAL(kFALSE);
    taskBFEemc->SetCentralitySelection(centMin,centMax);
    taskBFEemc->SelectCollisionCandidates(AliVEvent::kINT7);
    
    // Get the filename and make subfolders
    TString fileNameemc = mgr->AliAnalysisManager::GetCommonFileName();
    TString subContainerNameemc = ContNameExt;
    subContainerNameemc += "_INT7_EMCAL";
    AliAnalysisDataContainer *coutput3emc = mgr->CreateContainer(subContainerNameemc,TList::Class(),AliAnalysisManager::kOutputContainer,fileNameemc.Data());
        
    // connect the manager to task
    mgr->ConnectInput(taskBFEemc,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskBFEemc,1,coutput3emc);
    
    //FOR DCAL CLUSTERS....................................................................
    
    AliAnalysisTaskTPCCalBeauty* taskBFEdc = new AliAnalysisTaskTPCCalBeauty("bfedc");
    if(!taskBFEdc) return 0x0;
    // add your task to the manager
    mgr->AddTask(taskBFEdc);
    taskBFEdc->SetSSCut(applySSCut);
    taskBFEdc->SetClusterTypeEMC(kFALSE);
    taskBFEdc->SetClusterTypeDCAL(kTRUE);
    taskBFEdc->SetCentralitySelection(centMin,centMax);
    taskBFEdc->SelectCollisionCandidates(AliVEvent::kINT7);
    
    // Get the filename and make subfolders
    TString fileNamedc = mgr->AliAnalysisManager::GetCommonFileName();
    TString subContainerNamedc = ContNameExt;
    subContainerNamedc += "_INT7_DCAL";
    AliAnalysisDataContainer *coutput3dc = mgr->CreateContainer(subContainerNamedc,TList::Class(),AliAnalysisManager::kOutputContainer,fileNamedc.Data());
        
    // connect the manager to task
    mgr->ConnectInput(taskBFEdc,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskBFEdc,1,coutput3dc);
   
    /*
    //////////
    //  MB  //
    //////////
    AliAnalysisTaskTPCCalBeauty* taskBFEmb = new AliAnalysisTaskTPCCalBeauty("bfemb");
    if(!taskBFEmb) return 0x0;
    // add your task to the manager
    mgr->AddTask(taskBFEmb);
    taskBFEmb->SetClusterTypeEMC(ClsTypeEMC);
    taskBFEmb->SetClusterTypeDCAL(ClsTypeDCAL);
    taskBFEmb->SetCentralitySelection(centMin,centMax);
    taskBFEmb->SelectCollisionCandidates(AliVEvent::kMB);
    
    // Get the filename and make subfolders
    TString fileName02 = mgr->AliAnalysisManager::GetCommonFileName();
    TString subContainerName02 = ContNameExt;
    subContainerName02 += "BFE_PbPb_MB";
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(subContainerName02,TList::Class(),AliAnalysisManager::kOutputContainer,fileName02.Data());
    
    // connect the manager to task
    mgr->ConnectInput(taskBFEmb,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskBFEmb,1,coutput2);
    */
    
    ///////////////
    //  EGA EG1  //
    ///////////////
    
    //EMCAL TRIGGERED, EMCAL CLUSTERS.........................................................

    
    AliAnalysisTaskTPCCalBeauty* taskBFEeg01emc = new AliAnalysisTaskTPCCalBeauty("bfeeg01emc");
    if(!taskBFEeg01emc) return 0x0;
    // add your task to the manager
    mgr->AddTask(taskBFEeg01emc);
    taskBFEeg01emc->SetSSCut(applySSCut);
    taskBFEeg01emc->SetClusterTypeEMC(kTRUE);
    taskBFEeg01emc->SetClusterTypeDCAL(kFALSE);
    taskBFEeg01emc->SetCentralitySelection(centMin,centMax);
    taskBFEeg01emc->SetEMCalTriggerEG1(kTRUE);
    taskBFEeg01emc->SetEMCalTriggerDG1(kFALSE);
    
    // Get the filename and make subfolders
    TString fileNameEG01emc = mgr->AliAnalysisManager::GetCommonFileName();
    TString subContainerNameEG01emc = ContNameExt;
    subContainerNameEG01emc += "_TrigEG1_EMCAL";
    AliAnalysisDataContainer *coutputEG01emc = mgr->CreateContainer(subContainerNameEG01emc,TList::Class(),AliAnalysisManager::kOutputContainer,fileNameEG01emc.Data());
    
    // connect the manager to task
    mgr->ConnectInput(taskBFEeg01emc,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskBFEeg01emc,1,coutputEG01emc);
    
    //EMCAL TRIGGERED, DCAL CLUSTERS..........................................................
    
    /*AliAnalysisTaskTPCCalBeauty* taskBFEeg01dc = new AliAnalysisTaskTPCCalBeauty("bfeeg01dc");
    if(!taskBFEeg01dc) return 0x0;
    // add your task to the manager
    mgr->AddTask(taskBFEeg01dc);
    taskBFEeg01dc->SetClusterTypeEMC(kFALSE);
    taskBFEeg01dc->SetClusterTypeDCAL(kTRUE);
    taskBFEeg01dc->SetCentralitySelection(centMin,centMax);
    taskBFEeg01dc->SetEMCalTriggerEG1(kTRUE);
    taskBFEeg01emc->SetEMCalTriggerDG1(kFALSE);
    
    // Get the filename and make subfolders
    TString fileNameEG01dc = mgr->AliAnalysisManager::GetCommonFileName();
    TString subContainerNameEG01dc = ContNameExt;
    subContainerNameEG01dc += "BFE_PbPb_TrigEG1_DCClus";
    AliAnalysisDataContainer *coutputEG01dc = mgr->CreateContainer(subContainerNameEG01dc,TList::Class(),AliAnalysisManager::kOutputContainer,fileNameEG01dc.Data());
    
    // connect the manager to task
    mgr->ConnectInput(taskBFEeg01dc,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskBFEeg01dc,1,coutputEG01dc);

    //DCAL TRIGGERED, EMCAL CLUSTERS.........................................................
    
    AliAnalysisTaskTPCCalBeauty* taskBFEdg01emc = new AliAnalysisTaskTPCCalBeauty("bfedg01emc");
    if(!taskBFEdg01emc) return 0x0;
    // add your task to the manager
    mgr->AddTask(taskBFEdg01emc);
    taskBFEdg01emc->SetClusterTypeEMC(kTRUE);
    taskBFEdg01emc->SetClusterTypeDCAL(kFALSE);
    taskBFEdg01emc->SetCentralitySelection(centMin,centMax);
    taskBFEdg01emc->SetEMCalTriggerEG1(kFALSE);
    taskBFEdg01emc->SetEMCalTriggerDG1(kTRUE);
    
    // Get the filename and make subfolders
    TString fileNameDG01emc = mgr->AliAnalysisManager::GetCommonFileName();
    TString subContainerNameDG01emc = ContNameExt;
    subContainerNameDG01emc += "BFE_PbPb_TrigDG1_EMCClus";
    AliAnalysisDataContainer *coutputDG01emc = mgr->CreateContainer(subContainerNameDG01emc,TList::Class(),AliAnalysisManager::kOutputContainer,fileNameDG01emc.Data());
    
    // connect the manager to task
    mgr->ConnectInput(taskBFEdg01emc,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskBFEdg01emc,1,coutputDG01emc);
*/
    
    //DCAL TRIGGERED, DCAL CLUSTERS............................................................
    
    AliAnalysisTaskTPCCalBeauty* taskBFEdg01dc = new AliAnalysisTaskTPCCalBeauty("bfedg01dc");
    if(!taskBFEdg01dc) return 0x0;
    // add your task to the manager
    mgr->AddTask(taskBFEdg01dc);
    taskBFEdg01dc->SetSSCut(applySSCut);
    taskBFEdg01dc->SetClusterTypeEMC(kFALSE);
    taskBFEdg01dc->SetClusterTypeDCAL(kTRUE);
    taskBFEdg01dc->SetCentralitySelection(centMin,centMax);
    taskBFEdg01dc->SetEMCalTriggerEG1(kFALSE);
    taskBFEdg01dc->SetEMCalTriggerDG1(kTRUE);
        
    // Get the filename and make subfolders
    TString fileNameDG01dc = mgr->AliAnalysisManager::GetCommonFileName();
    TString subContainerNameDG01dc = ContNameExt;
    subContainerNameDG01dc += "_TrigDG1_DCAL";
    AliAnalysisDataContainer *coutputDG01dc = mgr->CreateContainer(subContainerNameDG01dc,TList::Class(),AliAnalysisManager::kOutputContainer,fileNameDG01dc.Data());
        
    // connect the manager to task
    mgr->ConnectInput(taskBFEdg01dc,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskBFEdg01dc,1,coutputDG01dc);
        
    
    return NULL;
}
