AliAnalysisHMMeanPt* AddHMMeanPt(
					int trackBit = 128,
					float vrtxz = 10.,
					float p_low = 0.15,
					float p_up = 2.0,
                    			int crossedrows = 70,
                    			TString suffix= ""
					)
{

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
        return 0x0;
    }
    
 	TString TaskMeanPt;
  	TaskMeanPt.Form("Taskfor_%d_%s", trackBit, suffix.Data());
  	AliAnalysisHMMeanPt *task = new AliAnalysisHMMeanPt(TaskMeanPt);
    	if(!task) return 0x0;
    
	//task->SelectCollisionCandidates(AliVEvent::kINT7);
	//task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
   
    mgr->AddTask(task);
        
    task->SettrackBit(trackBit); // setting FBs
    task->Seteventvrtx(vrtxz); // setting vrtx    
    task->Settrackpt(p_low, p_up); // setting pt    
    task->Settpcrows(crossedrows);	//setting TPCCrossedRows

    TString fileName = AliAnalysisManager::GetCommonFileName();
   	fileName += ":MeanpT_HM"; 

	//TString taskName = Form("%s_%s",name.Data(),suffix.Data());
	
    //connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("fOutputList_FB%d_%s",trackBit, suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));    
    //mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("fTreept%d_%s",trackBit, suffix.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(), fileName.Data())));    
	
    return task;

}


