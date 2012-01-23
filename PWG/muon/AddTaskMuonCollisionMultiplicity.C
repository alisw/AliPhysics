AliAnalysisTaskMuonCollisionMultiplicity *AddTaskMuonCollisionMultiplicity() 
{
  //
  // Task for the determination of the MUON trigger chamber efficiency
  //
  // lenhardt@subatech.in2p3.fr
  //

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }   

  const Char_t* fileName = "MUON.Multiplicity.THnSparse.root";
  
  // Create the task
  AliAnalysisTaskMuonCollisionMultiplicity* taskMuonMultiplicity = new AliAnalysisTaskMuonCollisionMultiplicity("MuonMultiplicity");

  taskMuonMultiplicity->SetZCut(10.0);
  taskMuonMultiplicity->SetEtaCut(1.6);

  // Add to the manager
  mgr->AddTask(taskMuonMultiplicity);
  /*
  AliAODInputHandler* aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);
  */
  
  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  /*
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  mcHandler->SetReadTR(kTRUE); 
  //mcHandler->SetInputPath("alien:///alice/sim/LHC10d3/117222/107/");
  */

  //
  // Create containers for input/output
    AliAnalysisDataContainer* cinput0  =	mgr->GetCommonInputContainer();

    AliAnalysisDataContainer *coutput0 =
	mgr->CreateContainer("Test", TList::Class(), AliAnalysisManager::kOutputContainer,  fileName);
    AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("Trigger", TList::Class(), AliAnalysisManager::kOutputContainer,  fileName);
    AliAnalysisDataContainer *coutput2 =
	mgr->CreateContainer("SingleMuon", TList::Class(), AliAnalysisManager::kOutputContainer,  fileName);
    AliAnalysisDataContainer *coutput3 =
	mgr->CreateContainer("Dimuon", TList::Class(), AliAnalysisManager::kOutputContainer,  fileName);
    AliAnalysisDataContainer *coutput4 =
	mgr->CreateContainer("MonteCarlo", TList::Class(), AliAnalysisManager::kOutputContainer,  fileName);

  // Attach input
    mgr->ConnectInput (taskMuonMultiplicity, 0, cinput0 );
  // Attach output
    mgr->ConnectOutput(taskMuonMultiplicity, 0, coutput0);
    mgr->ConnectOutput(taskMuonMultiplicity, 1, coutput1);
    mgr->ConnectOutput(taskMuonMultiplicity, 2, coutput2);
    mgr->ConnectOutput(taskMuonMultiplicity, 3, coutput3);
    mgr->ConnectOutput(taskMuonMultiplicity, 4, coutput4);

    return taskMuonMultiplicity;
}
