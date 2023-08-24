AliAnalysisMeanPt* AddMeanPt(
    					int trackBit = 128,
					float vrtxz1 = -10.,
                    		float vrtxz2 = 10.,
					float p_low = 0.15,
					float p_up = 2.0,
                    		int crossedrows = 70,
                    		TString period = "LHC18",
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
  	AliAnalysisMeanPt *task = new AliAnalysisMeanPt(TaskMeanPt);
    	if(!task) return 0x0;
    
	//task->SelectCollisionCandidates(AliVEvent::kINT7);
	//task->SelectCollisionCandidates(AliVEvent::kHighMultV0);
   
    TString treename = Form("fTreept_%s_%d_%s",period.Data(), trackBit, suffix.Data());

    mgr->AddTask(task);

    task->SettrackBit(trackBit); // setting FBs
    task->Seteventvrtx(vrtxz1, vrtxz2); // setting vrtx    
    task->Settrackpt(p_low, p_up); // setting pt    
    task->Settpcrows(crossedrows);	//setting TPCCrossedRows
    task->SetTreeName(treename); //name of tree

    TString fileName = AliAnalysisManager::GetCommonFileName();
   	fileName += ":MeanpT_MB"; 

    //connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("fList_%s_MB_FB%d_%s",period.Data(), trackBit, suffix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));    
    //mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("fTreept_%s_MB_FB%d_%s",period.Data(), trackBit, suffix.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));    
	//mgr->ConnectOutput(task,2,mgr->CreateContainer("fTreept", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));    

    return task;

}


