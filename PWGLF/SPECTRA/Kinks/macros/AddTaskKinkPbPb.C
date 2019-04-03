AliAnalysisTaskKinkPbPb* AddTaskKinkPbPb(TString lCustomName="", Float_t lRadiusKUp=200.0,  Float_t lRadiusKLow= 130.0, Int_t lNCluster=30, Float_t lLowQtValue=0.12, Float_t yRange=0.5, Float_t lnsigma=3.5, Float_t maxChi=4.0, Float_t Zpos=2.0, Bool_t strongAntiPile = kTRUE, Bool_t sevenSigma = kTRUE, Bool_t sixSigma = kFALSE, Bool_t eightSigma = kFALSE)
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
     
   
    	AliAnalysisTaskKinkPbPb  *task = new AliAnalysisTaskKinkPbPb("AliAnalysisTaskKinkPbPb", lRadiusKUp, lRadiusKLow, lNCluster, lLowQtValue, yRange, lnsigma, maxChi, Zpos, strongAntiPile, sevenSigma, sixSigma, eightSigma );
   
        task->SetKinkRadius(lRadiusKLow, lRadiusKUp);
   	task->SetNCluster(lNCluster);
	task->SetLowQtValue(lLowQtValue);
	task->SetYRange(yRange);
	task->SetnSigma(lnsigma);
	task->SetMaximumChiSquare(maxChi);	   
	task->SetDCAToVertexZpos(Zpos);
	task->SetStrongAntiPileup(strongAntiPile);
        task->SetDCAxySeven(sevenSigma);
        task->SetDCAxySix(sixSigma);
        task->SetDCAxyEight(eightSigma);
	mgr->AddTask(task);

        TString outputFileName = Form("%s:PWGLFSpectra.kinkPbPb", AliAnalysisManager::GetCommonFileName());
        TString outputname0 = Form("fListDefault_RadiusUp%.1f_RadiusLow%.1f_NCluster%i_Lowqt%2f_rapidity%1f_nsigma%1f_maxChi%1f_Zpos%1f_strongAntiPile%d_sevenSigma%d_sixSigma%d_eightSigma%d",lRadiusKUp, lRadiusKLow, lNCluster, lLowQtValue, yRange, lnsigma, maxChi, Zpos, strongAntiPile, sevenSigma, sixSigma, eightSigma);

        outputname0.Append(Form("%s",lCustomName.Data()));

        AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(outputname0, TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);

        mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
        mgr->ConnectOutput(task, 1, coutput1);

        return task; 
   }
