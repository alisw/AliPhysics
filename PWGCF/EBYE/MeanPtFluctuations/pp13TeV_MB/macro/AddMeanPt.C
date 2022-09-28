AliAnalysisMeanPt* AddMeanPt(
                             TString suffix= "fb128",
							 int trackBit = 128,
							 float vrtxz = 10.,
							 float p_low = 0.15,
							 float p_up = 2.0
							 )
{

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        return 0x0;
    }
    
    AliAnalysisMeanPt* task = new AliAnalysisMeanPt("TaskEbyE_buali");
    if(!task) return 0x0;
    
	task->SelectCollisionCandidates(AliVEvent::kINT7);
	//task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
   
    mgr->AddTask(task);
        
    task->SettrackBit(trackBit); // setting FBs
    task->Seteventvrtx(vrtxz); // setting vrtx    
    task->Settrackpt(p_low, p_up); // setting pt    


    TString fileName = "meanpt_pp_MB_";
    fileName += suffix.Data();

    //printf("container name: %s\n", fileName.Data());
  
    //connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("List%d_fOutputList_%s",trackBit, suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(), fileName.Data())));    
    mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("fTreept%d_%s",trackBit, suffix.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(), fileName.Data())));    
	
    return task;

}


