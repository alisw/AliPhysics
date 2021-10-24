///////////////////////////////////////////////////////////////////
//                                                               //            
// AddMyTask                                                     //
// Author: Redmer A. Bertens, Utrecht University, 2012           //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskDiHadCorrelHighPt* AddTaskDiHadCorrelHighPt(TString taskName = "name", Bool_t analysisMC = kFALSE, TString container_name_extension = "",TString fileName_extension = "",Bool_t useEff = kFALSE, TString EffFileNameWithPath = "")
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
    fileName += ":AliAnalysisTaskDiHadCorrelHighPt";      // create a subfolder in the file
    fileName += fileName_extension.Data();
    // now we create an instance of your task
    AliAnalysisTaskDiHadCorrelHighPt* task = new AliAnalysisTaskDiHadCorrelHighPt(taskName.Data(),analysisMC,useEff);
    if(!task) return 0x0;
    task->SetPtTrigMin(3);
    task->SetPtAsocMin(1);
    task->SetOStatus(1);
    task->SetCutsCrosscheck(kFALSE);
    // add your task to the manager
    mgr->AddTask(task);
     AliAnalysisDataContainer *cinput1 = 0x0;
     TList * effList = 0x0;

     if(useEff){
        TString eff_container_name = "Efficiency";
        eff_container_name+=container_name_extension.Data();
        cinput1 =  mgr->CreateContainer(Form("%s",eff_container_name.Data()),
                                    TList::Class(),
                                    AliAnalysisManager::kInputContainer);
        TFile * file = TFile::Open(Form("alien:///alice/cern.ch/user/%s.root",EffFileNameWithPath.Data()));
        if(!cinput1) printf("ERROR: Input container not created!\n");
        if(!file) {
            printf("ERROR: efficiency file %s.root is not available!\n",EffFileNameWithPath.Data());
        }
        effList = (TList*)file->Get(Form("fListEffHistos"));
        if(!effList){
            printf("ERROR: no efficiency list in %s available\n", EffFileNameWithPath.Data());
        }
    }
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    if(useEff){
        cinput1->SetData(effList);
        mgr->ConnectInput(task, 1, cinput1);  
    }
    
    // same for the output
    TString container_name = "MyOutputContainer";
    container_name += container_name_extension.Data();
    mgr->ConnectOutput(task,1,mgr->CreateContainer(container_name.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
