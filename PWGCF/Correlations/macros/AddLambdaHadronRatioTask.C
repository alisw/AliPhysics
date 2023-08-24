AliAnalysisTaskLambdaHadronRatio* AddLambdaHadronRatioTask(
  TString name = "lambdaHadronRatio",
  float multLow = 0, 
  float multHigh = 80,
  float trigBit = 1048576,
  float assocBit = 1024,
  float nSigmaTPC_proton = 2,
  float nSigmaTOF_proton = 2,
  float nSigmaTPC_pion = 3,
  float nSigmaTOF_pion = 3,
  bool tofVeto = true,
  TString effFilePath = "eff_out.root",
  TString centEstimator = "V0A"
  ) {

  // NOTE: The default arguments are placeholders ONLY, everything should be set within the run macro before function is called
  // 1024 is primary tracks (tight DCA cut, BIT(10))
  // 1048576 is kIsHybridGCG (BIT(20))

  AliAnalysisManager *manage = AliAnalysisManager::GetAnalysisManager();
  if (!manage) return 0x0;

  if(!manage->GetInputEventHandler()) return 0x0;



  TString file_name = AliAnalysisManager::GetCommonFileName();

  AliAnalysisTaskLambdaHadronRatio* task = new AliAnalysisTaskLambdaHadronRatio(name.Data());

  if(!task) return 0x0;

  task->SetMultBounds(multLow, multHigh);
  task->SetTriggerBit(trigBit);
  task->SetAssociatedBit(assocBit);
  task->SetPIDCuts(nSigmaTPC_proton, nSigmaTOF_proton, nSigmaTPC_pion, nSigmaTOF_pion, tofVeto);
  task->LoadEfficiencies(effFilePath);
  task->SetCentEstimator(centEstimator);

  manage->AddTask(task);

  manage->ConnectInput(task, 0, manage->GetCommonInputContainer());
  manage->ConnectOutput(task, 1, manage->CreateContainer("h-lambda", TList::Class(), AliAnalysisManager::kOutputContainer, file_name.Data()));

  return task;

}
