AliAnalysisTaskEPCalibForJet* AddTaskEPCalibForJet(
  const char *ntracks = "usedefault", 
  const char *nclusters = "usedefault",
  const char* ncells = "usedefault", 
  const char *suffix = ""
)
{
  return AliAnalysisTaskEPCalibForJet::AddTaskEPCalibForJet(
      ntracks, 
      nclusters,
      ncells, 
      suffix);
}

