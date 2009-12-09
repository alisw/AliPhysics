AliAnalysisTaskThreeJets * AddTaskThreeJets(char *bRec = "jets",char * bGen = "jetsAODMC_UA104")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskThreeJets", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskThreeJets", "This task requires an input event handler");
    return NULL;
   }
  
  // Create the task and configure it.
  //===========================================================================
  
  AliAnalysisTaskThreeJets * threeJets = new  AliAnalysisTaskThreeJets("Three Jet Analysis");
  
  threeJets->SetBranchRec(bRec);
  threeJets->SetBranchGen(bGen); 
  //  threeJets->SetDebugLevel(10);
  threeJets->SetR(.5); 
  //  threeJets->SetUseMC(kFALSE); // explicitly switch of use of MC/search for MC Jets

  
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if(type == "AOD"){
    threeJets->SetAODInput(kTRUE);
  }
  


  mgr->AddTask(threeJets);
   
  


      
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *coutput1_Corr = mgr->CreateContainer(Form("threeJets_%s_%s",bRec,bGen), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_threeJets_%s_%s",AliAnalysisManager::GetCommonFileName(),bRec,bGen));

  mgr->ConnectInput  (threeJets, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (threeJets, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (threeJets,  1, coutput1_Corr );
  
  return threeJets;
}


AliAnalysisTaskThreeJets * AddTaskJetCorrections(AliAnalysisManager* mgr,AliAnalysisDataContainer *cinput)
{
  if (!mgr) {
    ::Error("AddTaskJetSpectrum", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskJetSpectrum", "This task requires an input event handler");
    return NULL;
   }
  
  AliAnalysisTaskThreeJets * threeJets = new  AliAnalysisTaskThreeJets("ThreeJetAnalysis");
  
  threeJets->SetBranchGen("jetsMC"); 
  threeJets->SetBranchRec("jets");
  threeJets->SetR(.5); 
  mgr->AddTask(threeJets);

  AliAnalysisDataContainer * coutpu0 = mgr->CreateContainer("coutpu0", TTree::Class(),
				  AliAnalysisManager::kExchangeContainer);
  AliAnalysisDataContainer *coutput1_threeJets = mgr->CreateContainer("threeJets", TList::Class(),AliAnalysisManager::kOutputContainer,"threeJets.root");

   mgr->ConnectInput  (threeJets, 0, cinput);
   mgr->ConnectOutput (threeJets, 0, coutpu0);
   mgr->ConnectOutput (threeJets,  1, coutput1_Corr );
   
   return threeJets;
}

