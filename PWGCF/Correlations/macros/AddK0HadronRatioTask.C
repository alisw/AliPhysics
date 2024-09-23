AliAnalysisTaskK0HadronRatio *AddK0HadronRatioTask(
    TString name = "k0HadronRatio", 
    float multLow = 0, 
    float multHigh = 80,
    float trigBit = 1048576, 
    float assocBit = 1024, 
    float nSigmaTPC_pion = 3,
    float nSigmaTOF_pion = 3,
    bool tofVeto = true, 
    TString effFilePath = "eff_out_mult_eta_dep.root", 
    TString centEstimator = "V0A") {

  // NOTE: The default arguments are placeholders ONLY, everything should be set
  // within the run macro before function is called 1024 is primary tracks
  // (tight DCA cut, BIT(10)) 1048576 is kIsHybridGCG (BIT(20))


  AliAnalysisManager *manage = AliAnalysisManager::GetAnalysisManager();
  if (!manage)
    return 0x0;

  if (!manage->GetInputEventHandler())
    return 0x0;

  AliAnalysisTaskK0HadronRatio *task =
      new AliAnalysisTaskK0HadronRatio(name.Data());

  if (!task)
    return 0x0;

  task->SetMultBounds(multLow, multHigh);
  task->SetTriggerBit(trigBit);
  task->SetAssociatedBit(assocBit);
  task->SetPIDCuts(nSigmaTPC_pion, nSigmaTOF_pion, tofVeto);
  task->LoadEfficiencies(effFilePath);
  task->SetCentEstimator(centEstimator);

  manage->AddTask(task);
    
  TString file_name = AliAnalysisManager::GetCommonFileName();


  manage->ConnectInput(task, 0, manage->GetCommonInputContainer());
  manage->ConnectOutput(
      task, 1,
      manage->CreateContainer("h-k0", TList::Class(),
                              AliAnalysisManager::kOutputContainer,
                              file_name.Data()));

  return task;
}
