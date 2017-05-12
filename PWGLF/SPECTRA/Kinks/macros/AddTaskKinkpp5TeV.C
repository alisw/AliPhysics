AliAnalysisTaskKinkpp5TeV* AddTaskKinkpp5TeV(TString lCustomName="",TString dirName="",Float_t yRange=.5, Float_t lRadiusKLow= 130., Float_t lRadiusKUp=200., Int_t lNCluster=30., Float_t lLowQtValue=0.12)
   {
     //pp settings         
      	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      	if (!mgr) 
        {
          ::Error("AddKinkTaskpp5tev", "No analysis manager to connect to.");
          return NULL;
        }   
     // Check the analysis type using the event handlers connected to the analysis manager.
     //==============================================================================
     	if (!mgr->GetInputEventHandler()) 
       	{
         ::Error("AddKinkTaskpp5tev", "This task requires an input event handler");
        	 return NULL;
      	 }	   
     
     	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
     	if(type.Contains("AOD"))
       	{
         ::Error("AddKinkTaskpp5tev", "This task requires to run on ESD");
      	   return NULL;
       	}
     
     //TString outputFileName = AliAnalysisManager::GetCommonFileName();
     //outputFileName += ":PWG2SpectraTOF";
   
    	AliAnalysisTaskKinkpp5TeV  *task = new AliAnalysisTaskKinkpp5TeV(Form("TaskKink_%s",dirName.Data()));
   
    //task->SetMC("kFALSE"); // 26/11/12
   
   	task->SetYRange(yRange);
   	task->SetKinkRadius(lRadiusKLow, lRadiusKUp);
   	task->SetNCluster(lNCluster);
   	task->SetLowQtValue(lLowQtValue);
	   
	mgr->AddTask(task);
   
     //Attach input
     	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer(); 
   //  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());     
      	mgr->ConnectInput(task,0,cinput);
     
     	TString lContainerName="PWGLFKinks5TeV";
     	lContainerName.Append(lCustomName);
     	AliAnalysisDataContainer *coutput1= mgr->CreateContainer(lContainerName.Data(),TList::Class(), AliAnalysisManager::kOutputContainer,"MCkinkspp5TeV.root");
     	mgr->ConnectOutput(task, 1, coutput1);
    
     
     	return task;
     
   }
