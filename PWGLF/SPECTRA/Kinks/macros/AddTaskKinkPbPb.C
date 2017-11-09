AliAnalysisTaskKinkPbPb* AddTaskKinkPbPb(TString lCustomName="", Float_t lRadiusKUp=200.0,  Float_t lRadiusKLow= 130.0, Int_t lNCluster=30, Float_t lLowQtValue=0.12, Float_t yRange=0.5)
   {
     //pbpb settings         
      	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      	if (!mgr) 
        {
          ::Error("AddTaskKinkPbPb", "No analysis manager to connect to.");
          return NULL;
        }   
     // Check the analysis type using the event handlers connected to the analysis manager.
     //==============================================================================
     	if (!mgr->GetInputEventHandler()) 
       	{
         ::Error("AddTaskKinkPbPb", "This task requires an input event handler");
        	 return NULL;
      	 }	   
     
     	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
     	if(type.Contains("AOD"))
       	{
         ::Error("AddTaskKinkPbPb", "This task requires to run on ESD");
      	   return NULL;
       	}
     
   
    	AliAnalysisTaskKinkPbPb  *task = new AliAnalysisTaskKinkPbPb("AliAnalysisTaskKinkPbPb", lRadiusKUp, lRadiusKLow, lNCluster, lLowQtValue, yRange);
   
        task->SetKinkRadius(lRadiusKLow, lRadiusKUp);
   	task->SetNCluster(lNCluster);
	task->SetLowQtValue(lLowQtValue);
	task->SetYRange(yRange);
	   
	mgr->AddTask(task);

        TString outputFileName = Form("%s:PWGLFSpectra.kinkPbPb", AliAnalysisManager::GetCommonFileName());
        TString outputname0 = Form("fListDefault_RadiusUp%.1f_RadiusLow%.1f_NCluster%i_Lowqt%2f_rapidity%1f",lRadiusKUp, lRadiusKLow, lNCluster, lLowQtValue, yRange);

        outputname0.Append(Form("%s",lCustomName.Data()));

        AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(outputname0, TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);

        mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
        mgr->ConnectOutput(task, 1, coutput1);

        return task; 
   }
