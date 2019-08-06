///////////////////////////////////////////////////////////////////
//                                                               //            
// AddFlowTPCEMCalRun2                                           //
// Author: K. Tadokoro, Univ. of Tsukuba, 2019                   //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
//class AliRDHFCuts;
AliAnalysisTask *AddTaskFlowTPCEMCalRun2(
    TString OADBfilename="",
    TString Splinefilename="", 
    Bool_t iMC = kFALSE,
    TString ContNameExt= "semicentral")
{
    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function
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
    AliAnalysisTaskFlowTPCEMCalRun2* task = new AliAnalysisTaskFlowTPCEMCalRun2("task");   
    if(!task) return 0x0;
    // add your task to the manager
    task->SetOADBFileName(OADBfilename);
    task->SetqnPercentileSelection(Splinefilename);
    if(iMC)
      {
       task->SelectCollisionCandidates(AliVEvent::kMB);
      }
    else
      {
       task->SelectCollisionCandidates(AliVEvent::kSemiCentral);
      }

    mgr->AddTask(task);

    TString containerName = mgr->GetCommonFileName();
    containerName += ":PWGHF_HFEflowTPCEMCalRun2";
    containerName += ContNameExt;
    TString SubcontainerName = Form("HFEflowTPCEMCalRun2");
    SubcontainerName += ContNameExt;
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName, TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput1);
    


    return NULL;
}
