/**
 * This task only tests the functionality of the
 * iterators of the containers. Therefore cuts applied
 * are of no physical meaning but are chosen in a way
 * that the amount of selected objets should differ the
 * amount of all objects in a good fraction of cases
 * in order to test the functionality of the accept_iterators.
 * @param nameClusterContainer Name of the cluster container used for the test
 * @param nameMCParticleContainer Name of the track container used for the test
 * @param nameTrackContainer Name of the track container
 * @return Pointer to the test task
 */
AliAnalysisTaskEmcalIteratorTest *AddTaskEmcalTestIterators(
    TString nameClusterContainer = "",
    TString nameMCParticleContainer = "",
    TString nameTrackContainer = "",
    TString period = ""
    )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliAnalysisTaskEmcalIteratorTest *testtask = new AliAnalysisTaskEmcalIteratorTest("emcalIteratorTest");
  mgr->AddTask(testtask);

  testtask->SelectCollisionCandidates(AliVEvent::kINT7);
  testtask->SetUseAliAnaUtils(true, 1);

  AliClusterContainer *clustercont = testtask->AddClusterContainer(nameClusterContainer.Data());
  testtask->SetClusterContainerName(nameClusterContainer);
  clustercont->SetClusECut(1.);
  clustercont->SetEtaLimits(-0.5, 0.5);

  AliMCParticleContainer *mccont = testtask->AddMCParticleContainer(nameMCParticleContainer.Data());
  testtask->SetMCParticleContainerName(nameMCParticleContainer);
  mccont->SetEtaLimits(-0.5, 0.5);
  mccont->SetMinPt(0.9);

  AliTrackContainer *trackcont = testtask->AddTrackContainer(nameTrackContainer.Data());
  testtask->SetTrackContainerName(nameTrackContainer);
  trackcont->SetEtaLimits(-0.5, 0.5);
  trackcont->SetMinPt(0.8);
  trackcont->SetTrackFilterType(AliEmcalTrackSelection::kHybridTracks);
  trackcont->SetTrackCutsPeriod(period.Data());

  TString filename = mgr->GetCommonFileName();
  filename += ":iterator_test";
  mgr->ConnectInput(testtask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(testtask, 1, mgr->CreateContainer("testresults", TList::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));

  return testtask;

}
