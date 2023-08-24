AliAnalysisMCMeanPt* AddMCMeanPt(float zvtxcut1     = -10, float zvtxcut2     = 10, int trackBit    = 768, int MaxTPCclus  = 70., TString suffixName = "")
{
    // an instance of the class to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }

 	TString TaskMCMeanPt;
  	TaskMCMeanPt.Form("Taskfor_%d_%s", trackBit, suffixName.Data());
  	AliAnalysisMCMeanPt *task = new AliAnalysisMCMeanPt(TaskMCMeanPt);
    	if(!task) return 0x0;

    task->Setzvtxcut(zvtxcut1, zvtxcut2);
    task->SettrackBit(trackBit);
    task->SetMaxTPCCluster(MaxTPCclus);
    //task->SetNTPCCluster(TPCNclus);
    mgr->AddTask(task);
    
    TString finDirname      =  suffixName.Data();

    TString fileName = mgr->GetCommonFileName(); 

    printf("container name: %s\n", fileName.Data());
    printf("file name: %s\n", finDirname.Data());

    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("%s_fb%d_vz%3.2f_tpc%d_fOutputList",finDirname.Data(), trackBit, zvtxcut2, MaxTPCclus), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    return task;

}
