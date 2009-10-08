AliAnalysisTaskSED0 *AddTaskD0(Int_t flag=0/*0 = D0,1 = LS*/)
{
  //
  // AddTask for the AliAnalysisTaskSE for D0 candidates
  // invariant mass histogram and association with MC truth 
  // (using MC info in AOD) and cut variables distributions
  // C.Bianchin  chiara.bianchin@pd.infn.it
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskD0Distr", "No analysis manager to connect to.");
    return NULL;
  }   

  TString filename="";
  if(flag==0){
    filename="D0InvMass.root"; 
  } else filename="LSD0.root";

  // Aanalysis task    
  AliAnalysisTaskSED0 *massD0Task = new AliAnalysisTaskSED0("D0MassAndDistrAnalysis");
  massD0Task->SetDebugLevel(2);
  massD0Task->SetArray(flag);
  mgr->AddTask(massD0Task);
  
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputmassD0 = mgr->CreateContainer("cinputmassD0",TChain::Class(), 
							  AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutputmassD01 = mgr->CreateContainer("coutputmassD01",TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   filename.Data());
  AliAnalysisDataContainer *coutputmassD02 = mgr->CreateContainer("coutputmassD02",TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   filename.Data());
  AliAnalysisDataContainer *coutputmassD03 = mgr->CreateContainer("nEntriesD0",TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   filename.Data());
  AliAnalysisDataContainer *coutputmassD04 = mgr->CreateContainer("coutputmassD0distr",TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   filename.Data());

  mgr->ConnectInput(massD0Task,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(massD0Task,1,coutputmassD01);
  mgr->ConnectOutput(massD0Task,2,coutputmassD02);
  mgr->ConnectOutput(massD0Task,3,coutputmassD03);
  mgr->ConnectOutput(massD0Task,4,coutputmassD04);

  return massD0Task;
}
