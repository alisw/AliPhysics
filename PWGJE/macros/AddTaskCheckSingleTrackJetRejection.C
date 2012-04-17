AliAnalysisTaskCheckSingleTrackJetRejection *AddTaskCheckSingleTrackJetRejection(char *jf="ANTIKT",Float_t radius=0.4,UInt_t filter=256,Int_t backM=0,Float_t tPtcut=0.15,Int_t skipCone=0,Bool_t IsMC=true)
{

   // Creates a JetQA task, configures it and adds it to the analysis manager.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskCheckSingleTrackJetRejection", "No analysis manager to connect to.");
      return NULL;
   }

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCheckSingleTrackJetRejection", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================

   AliAnalysisTaskCheckSingleTrackJetRejection *jetqamcana = new AliAnalysisTaskCheckSingleTrackJetRejection("TaskCheckSingleTrackJetRejection");
   jetqamcana->SetDebugLevel(0);
	 jetqamcana->SetAlgorithm(jf);
	 jetqamcana->SetRadius(radius);
	 jetqamcana->SetFilterMask(filter);
	 jetqamcana->SetBackSubMode(backM);
	 jetqamcana->SetTrackPtCut(tPtcut);
	 jetqamcana->SetSkipCone(skipCone);
	 jetqamcana->SetMC(IsMC);
	 mgr->AddTask(jetqamcana); 

	 TString cAdd = "";
	 cAdd += Form("%02d_",(int)((radius+0.01)*10.));
	 cAdd += Form("B%d",(int)backM);
	 cAdd += Form("_Filter%05d",filter);
	 cAdd += Form("_Cut%05d",(int)(1000.*tPtcut));
	 cAdd += Form("_Skip%02d",skipCone);
	 TString Branch;
	 if(IsMC)Branch = Form("MC_clustersAOD_%s%s",jf,cAdd.Data());
	 else    Branch = Form("Data_clustersAOD_%s%s",jf,cAdd.Data());

   AliAnalysisDataContainer *cout_jetsqamc = mgr->CreateContainer("histlist", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_CheckSingleTrackJetRejection_%s",AliAnalysisManager::GetCommonFileName(),Branch.Data()));

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   mgr->ConnectInput (jetqamcana,0, mgr->GetCommonInputContainer());  
   mgr->ConnectOutput(jetqamcana,0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput(jetqamcana,1, cout_jetsqamc);


   return jetqamcana;
}
