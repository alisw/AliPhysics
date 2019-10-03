AliAnalysisTaskCountLcEta *AddTaskCountLcEta(TString type="AOD",Float_t eta=0.9, TString suffix="") {
  //

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
   Error("AddTaskCountLcEta", "No analysis manager found.");
   return 0;
  }
  TString filename =AliAnalysisManager::GetCommonFileName();
  filename += ":WP1ITSUp_Lc";
  const Int_t ncuts=3;
  //Double_t cuts[ncuts]={0.8,0.5,0.8}; // default cuts
  Double_t cuts[ncuts]={0.8,0.8,0.8};   // a bit tighter cuts
  //Double_t cuts[ncuts]={2.,2.,2.};    // for testing
  Printf("CUTS on pt = %f, \%f, %f",cuts[0],cuts[1],cuts[2]);
  //TString cutnames[ncuts]={"ptpi","ptK","ptp"};
  // Create the task, add it to the manager and configure it.
  //===========================================================================
  AliAnalysisTaskCountLcEta *hfTask = new AliAnalysisTaskCountLcEta("TaskCountLcEta",ncuts,cuts);
  hfTask->SetDataType(type);
  hfTask->SetEtaAbs(eta);
  //hfTask->SetCutNames(ncuts,cutnames);
  
  mgr->AddTask(hfTask);
  TString nameout="CountLcEta";
  nameout+=suffix;
 //AliAnalysisDataContainer *cinput= (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("cAUTO_INPUT");
 
  AliAnalysisDataContainer *cinputLambdac = mgr->CreateContainer(Form("cinputLc%s",suffix.Data()),TChain::Class(),
								 AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(hfTask,0,mgr->GetCommonInputContainer());
  
  AliAnalysisDataContainer *coutput1= mgr->CreateContainer(nameout.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,filename.Data());
 // AliAnalysisDataContainer *coutput1= mgr->GetCommonInputContainer();
  ////
  // Create containers for input/output
  mgr->ConnectOutput(hfTask,1,coutput1);

  return hfTask;
}
