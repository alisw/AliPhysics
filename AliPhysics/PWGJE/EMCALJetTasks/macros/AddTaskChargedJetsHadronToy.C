// AddTaskChargedJetsHadronToy.C

AliAnalysisTaskChargedJetsHadronToy* AddTaskChargedJetsHadronToy(
  const char *trackOutput        = "tracks_toy"
)
{
  cout << " ############ MACRO EXECUTION STARTED: AddTaskChargedJetsHadronToy.C ############\n";
  //==============================================================================
  // Prepare analysis manager, containers, etc.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr)
  {
    ::Error("AddTaskChargedJetsHadronToy", "No analysis manager to connect to.");
    return NULL;
  }  
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskChargedJetsHadronToy", "This task requires an input event handler");
    return NULL;
  }
  
  //==============================================================================
  // Adding and configuring tasks

  AliAnalysisDataContainer* contHistos = mgr->CreateContainer(Form("Toy_histos_%s", trackOutput), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));

  AliAnalysisTaskChargedJetsHadronToy* toyModel = new AliAnalysisTaskChargedJetsHadronToy();
  toyModel->SetOutputArrayName(trackOutput);
  mgr->AddTask(toyModel);

  //==============================================================================
  // Finalization

  mgr->ConnectInput  (toyModel, 0, mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (toyModel, 1, contHistos );

  cout << " ############ MACRO EXECUTION DONE: AddTaskChargedJetsHadronToy.C ############\n";
 
  return toyModel;
}
