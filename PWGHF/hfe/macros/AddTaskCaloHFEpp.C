///////////////////////////////////////////////////////////////////
//                                                               //            
// AddCaloHFEpp                                                  //
// Author: T. Suzuki Univ. of Tsukuba                            //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskCaloHFEpp* AddTaskCaloHFEpp(TString name = "name",
		                 TString dataname = "dataname",
		                 Bool_t flagEG1 = kFALSE,
		                 Bool_t flagEG2 = kFALSE,
		                 Bool_t flagDG1 = kFALSE,
		                 Bool_t flagDG2 = kFALSE,
		                 Bool_t SetFlagClsTypeEMC = kFALSE,
		                 Bool_t SetFlagClsTypeDCAL = kFALSE)
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
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":CaloHFEpp";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskCaloHFEpp* task = new AliAnalysisTaskCaloHFEpp(name.Data());   
    task -> SetEG1(flagEG1);
    task -> SetEG2(flagEG2);
    task -> SetDG1(flagDG1);
    task -> SetDG2(flagDG2);
    task -> SetfFlagClsTypeEMC(SetFlagClsTypeEMC);
    task -> SetfFlagClsTypeDCAL(SetFlagClsTypeDCAL);
    if(!task) return 0x0;

    // add your task to the manager
    mgr->AddTask(task);

    TString containerName = mgr->GetCommonFileName();
    containerName += ":PWGHF_hfeCalpp";
    TString SubcontainerName = Form("hfeCalpp");
    SubcontainerName += name;
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName, TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput1); 

    return task;
}
