// AddTaskChargedJetsHadronToy.C

AliAnalysisTaskChargedJetsHadronToy* AddTaskChargedJetsHadronToy(
  const char *trackInput         = "tracks",
  const char *trackOutput        = "Toymodel_tracks",
  const char *jetOutput          = "Toymodel_jets_generated",
  Bool_t      createUE           = kTRUE,
  Bool_t      createJets         = kTRUE
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

  AliAnalysisTaskChargedJetsHadronToy* toyModel = new AliAnalysisTaskChargedJetsHadronToy();
  toyModel->SetInputTracksName(trackInput);
  toyModel->SetOutputTracksName(trackOutput);
  toyModel->SetGeneratedJetsName(jetOutput);
  toyModel->SetCreateUE(createUE);
  toyModel->SetCreateJets(createJets);

  mgr->AddTask(toyModel);

  //==============================================================================
  // Finalization

  mgr->ConnectInput(toyModel, 0,  mgr->GetCommonInputContainer() );
 
  cout << " ############ MACRO EXECUTION DONE: AddTaskChargedJetsHadronToy.C ############\n";
 
  return toyModel;
}
