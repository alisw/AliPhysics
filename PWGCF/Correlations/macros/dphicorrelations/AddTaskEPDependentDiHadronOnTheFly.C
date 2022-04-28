AliEPDependentDiHadronOnTheFlyMCTask* AddTaskEPDependentDiHadronOnTheFly()
{
// Retrieving analysis manager.
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    ::Error("AddTaskEPDependentDiHadronOnTheFly", "No analysis manager to connect to.");
    return NULL;
  }

// Creating task.
  AliEPDependentDiHadronOnTheFlyMCTask* evplanetask = new AliEPDependentDiHadronOnTheFlyMCTask("AliEPDependentDiHadronOnTheFlyMCTask");
  mgr->AddTask(evplanetask);

// Attaching in- and output.
  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer("AliEPDependentDiHadronOnTheFlyMCTask", TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  mgr->ConnectInput(evplanetask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(evplanetask, 1, coutput1);

// 'Default' settings used in the initial analysis.
  evplanetask->SetTrackEtaCut(2.0);
  evplanetask->SetCustomBinning("p_t_leading: 3., 4., 6., 8., 15., 25.\np_t_assoc: 0.5, 1., 2., 3., 4., 6., 8., 15., 25.\np_t_leading_course: 3., 4., 6., 8., 15., 25.\nvertex_eff: -9., -7., -5., -3., -1., 1., 3., 5., 7., 9.\nvertex: -7., -5., -3., -1., 1., 3., 5., 7.\ndelta_phi: -1.570796, -1.472622, -1.374447, -1.276272, -1.178097, -1.079923, -0.9817477, -0.8835729, -0.7853982, -0.6872234, -0.5890486, -0.4908739, -0.3926991, -0.2945243, -0.1963495, -0.09817477, 0, 0.09817477, 0.1963495, 0.2945243, 0.3926991, 0.4908739, 0.5890486, 0.6872234, 0.7853982, 0.8835729, 0.9817477, 1.079923, 1.178097, 1.276272, 1.374447, 1.472622, 1.570796, 1.668971, 1.767146, 1.865321, 1.963495, 2.061670, 2.159845, 2.258020, 2.356195, 2.454369, 2.552544, 2.650719, 2.748894, 2.847068, 2.945243, 3.043418, 3.141593, 3.239767, 3.337942, 3.436117, 3.534292, 3.632467, 3.730641, 3.828816, 3.926991, 4.025166, 4.123340, 4.221515, 4.319690, 4.417865, 4.516039, 4.614214, 4.712389\nmultiplicity: 0., 5, 20., 50., 100.1\ndelta_eta: -4.0, -3.9, -3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1, -3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0");

  return evplanetask;
}
