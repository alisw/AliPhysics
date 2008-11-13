AliAnalysisManager *myRsnAnalysis() {
  // Create analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("My Rsn Analysis");
  if (!mgr) AliError("no Analysis Mgr");

  // Create Resonance analysis task
  AliRsnAnalysisSE *task = new AliRsnAnalysisSE("AliRsnAnalysisSE");
  // sets input type
//   task->SetInputType(AliRsnAnalysisTaskSEBase::kESD,mgr,kTRUE);
  task->SetInputType(AliRsnAnalysisTaskSEBase::kESDMC,mgr,kTRUE);
//   task->GetESDHandler()->SetInactiveBranches("*Calo*");
  
  // add our config file
  task->AddPairMgrFromConfig("myRsnConfig.C");

  // sets input output
  AliAnalysisDataContainer *input  = mgr->CreateContainer("in", TChain::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *dummy  = mgr->CreateContainer("dummy", TTree::Class(), AliAnalysisManager::kOutputContainer, "default");
  AliAnalysisDataContainer *output = mgr->CreateContainer("Histograms", TList::Class(), AliAnalysisManager::kOutputContainer, "rsnHistOut.root");

  // Reader settings
  AliRsnReader *reader = task->GetReader();

  // PID settings
  AliRsnPID *pid = task->GetPID();
  pid->SetPriorProbability(AliRsnPID::kElectron, 0.02);
  pid->SetPriorProbability(AliRsnPID::kMuon,     0.02);
  pid->SetPriorProbability(AliRsnPID::kPion,     0.83);
  pid->SetPriorProbability(AliRsnPID::kKaon,     0.07);
  pid->SetPriorProbability(AliRsnPID::kProton,   0.06);
  pid->SetMaxPt(10.0);
  pid->SetMinProb(0.5);

  // connect containers to AnalysisManager
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, input);
  mgr->ConnectOutput(task, 0, dummy);
  mgr->ConnectOutput(task, 1, output);

  // return analysis manager
  return mgr;
}
