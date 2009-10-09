AliAnalysisTaskSED0Mass *AddTaskD0Mass(Int_t flag=0/*0 = D0,1 = LS*/)
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

  TString filename="",out1name="",out2name="",out3name="",out4name="";
  if(flag==0){
    filename="D0InvMass.root"; 
    out1name="coutputmassD01";
    out2name="coutputmassD02";
    out3name="nEntriesD0";
    out4name="coutputmassD0distr";
  } else {
    filename="LSD0.root";
    out1name="coutputmassLS1";
    out2name="coutputmassLS2";
    out3name="nEntriesLS";
    out4name="coutputmassLSdistr";

  }

  // Aanalysis task    
  AliAnalysisTaskSED0Mass *massD0Task = new AliAnalysisTaskSED0Mass("D0MassAndDistrAnalysis");
  massD0Task->SetDebugLevel(0);
  massD0Task->SetArray(flag);
  mgr->AddTask(massD0Task);
  
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputmassD0 = mgr->CreateContainer("cinputmassD0",TChain::Class(), 
							  AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutputmassD01 = mgr->CreateContainer(out1name,TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   filename.Data());
  AliAnalysisDataContainer *coutputmassD02 = mgr->CreateContainer(out2name,TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   filename.Data());
  AliAnalysisDataContainer *coutputmassD03 = mgr->CreateContainer(out3name,TH1F::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   filename.Data());
  AliAnalysisDataContainer *coutputmassD04 = mgr->CreateContainer(out4name,TList::Class(),
							   AliAnalysisManager::kOutputContainer, 
							   filename.Data());

  mgr->ConnectInput(massD0Task,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(massD0Task,1,coutputmassD01);
  mgr->ConnectOutput(massD0Task,2,coutputmassD02);
  mgr->ConnectOutput(massD0Task,3,coutputmassD03);
  mgr->ConnectOutput(massD0Task,4,coutputmassD04);


  return massD0Task;
}
