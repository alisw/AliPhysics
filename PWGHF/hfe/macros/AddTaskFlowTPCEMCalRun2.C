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
    Double_t cmim = 30.0,
    Double_t cmax = 50.0,
    TString ContNameExt= "semicentral",
    Double_t tpcnsig = -1.0,
    Double_t emceop = 0.9,
    Double_t emcss_mim = 0.01,
    Double_t emcss_max = 0.35,
    Double_t invmass = 0.1,
    Double_t invmass_pt = 0.2,
    Bool_t cent = kFALSE,
    Bool_t semi = kTRUE,
    Bool_t TreeOn = kFALSE
 )
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
      Int_t centID = 0;
      if(cmim==0.0)centID = 0;
      if(cmim==30.0)centID = 1;

       if(centID==0)
         {
          cout << "centrality selection : " << cmim << " ; " << cmax << " ; kCentral " << endl; 
          task->SelectCollisionCandidates(AliVEvent::kCentral);
          }
       else if(centID==1)
         {
          cout << "centrality selection : " << cmim << " ; " << cmax << " ; kSemiCentral " << endl; 
          task->SelectCollisionCandidates(AliVEvent::kSemiCentral);
          }
       else
          {
           task->SelectCollisionCandidates(AliVEvent::kINT7);
          }
     }

    task->SetMinCentrality(cmim);
    task->SetMaxCentrality(cmax);
    task->SetPIDcuts(tpcnsig, emceop, emcss_mim, emcss_max);
    task->SetMasscuts(invmass,invmass_pt);
    task->SetMCCentral(cent);
    task->SetMCSemiCentral(semi);
    task->SetTree(TreeOn);

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
