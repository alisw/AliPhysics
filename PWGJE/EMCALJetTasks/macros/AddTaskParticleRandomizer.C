// AddTaskParticleRandomizer.C

AliAnalysisTaskParticleRandomizer* AddTaskParticleRandomizer(
  const char *inputParticles     = "tracks",
  const char *outputParticles    = "tracks_randomized",
  Bool_t      randomizeInPhi     = kTRUE,
  Bool_t      randomizeInEta     = kFALSE
)
{  
  cout << " ############ MACRO EXECUTION STARTED: AddTaskParticleRandomizer.C ############\n";
  //==============================================================================
  // Prepare analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskParticleRandomizer", "No analysis manager to connect to.");
    return NULL;
  }
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskParticleRandomizer", "This task requires an input event handler");
    return NULL;
  }

  //==============================================================================
  // Adding and configuring tasks

  AliAnalysisTaskParticleRandomizer* randomizer = new AliAnalysisTaskParticleRandomizer();
  randomizer->SetInputArrayName(inputParticles);
  randomizer->SetOutputArrayName(outputParticles);
  randomizer->SetRandomizeInPhi(randomizeInPhi);
  randomizer->SetRandomizeInEta(randomizeInEta);

  mgr->AddTask(randomizer);

  //==============================================================================
  // Finalization

  mgr->ConnectInput(randomizer, 0, mgr->GetCommonInputContainer() );

  cout << " ############ MACRO EXECUTION DONE: AddTaskParticleRandomizer.C ############\n";

  return randomizer;
}
