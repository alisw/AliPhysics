AliAnalysisTaskdEdxSSDQA* AddTaskdEdxSSDQA (Float_t pcut=1.2)
{
  	 AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  	 if (!mgr) 
	 {
      		Printf("AliAnalysisTaskdEdxSSDQA No analysis manager to connect to.");
     		 return NULL;
   	}   
  	 if (!mgr->GetInputEventHandler()) 
	 {
      		Printf("AliAnalysisTaskdEdxSSDQA  no input event handler");
     		 return NULL;
   	} 
	TString type = mgr->GetInputEventHandler()->GetDataType();
	if(type!="ESD")
	{
		Printf("AliAnalysisTaskdEdxSSDQA  no ESD input event handler");
     		 return NULL;
	}
	AliESDInputHandler  * esdH =(AliESDInputHandler *) mgr->GetInputEventHandler();
	esdH ->SetReadFriends(1);
	AliAnalysisTaskdEdxSSDQA* taskdEdxSSDQA=new AliAnalysisTaskdEdxSSDQA();
	taskdEdxSSDQA->SetPcut(pcut);
	    
	mgr->AddTask(taskdEdxSSDQA);

	TString outputFileName = AliAnalysisManager::GetCommonFileName();
	outputFileName+=":PWG1dEdxSSDQA";
	AliAnalysisDataContainer *cinput1  = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
 	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("output", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
  
  //________________________________________________//
	mgr->ConnectInput(taskdEdxSSDQA,0,mgr->GetCommonInputContainer());
  	mgr->ConnectOutput(taskdEdxSSDQA,1,coutput1);	 
	return taskdEdxSSDQA;
}
