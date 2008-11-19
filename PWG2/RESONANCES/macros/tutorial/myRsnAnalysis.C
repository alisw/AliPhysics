//
// *** myRsnAnalysis.C ***
//
// This macro initializes the classes needed to run the analysis:
//  - AliAnalysisManager
//  - AliAnalysisTask's for resonances
//  - input and output containers
//
// For the tutorials, it is recommended not to modify this macro,
// but for real analysis, one could want to change the following lines:
//  - "task->SetInputType(...)" to play with ESD with/without MC, AOD, RSN, MC only
//  - adding options to the "AliRsnReader"
//  - modifying properties of the "AliRsnPID"
//
AliAnalysisManager *myRsnAnalysis(const char *configMacro = 0)
{
  // check argument
  if (!configMacro) {
    cerr << "== Error in myRsnAnalysis.C ==" << endl;
    cerr << "This macro requires string argument containing the NAME of the CONFIG MACRO" << endl;
    cerr << "Since this macro is called by another, you should edit it and add the required value as default." << endl;
    return 0x0;
  }

  // Create analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("My Rsn Analysis");
  if (!mgr) AliError("no Analysis Mgr");

  // Create Resonance analysis task
  AliRsnAnalysisSE *task = new AliRsnAnalysisSE("AliRsnAnalysisSE");
  // sets input type
  // task->SetInputType(AliRsnAnalysisTaskSEBase::kESD,mgr,kTRUE);
  task->SetInputType(AliRsnAnalysisTaskSEBase::kESDMC,mgr,kTRUE);
  // task->GetESDHandler()->SetInactiveBranches("*Calo*");

  // add our config file
  //task->AddPairMgrFromConfig("myRsnConfig.C");
  task->AddPairMgrFromConfig(configMacro);

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
