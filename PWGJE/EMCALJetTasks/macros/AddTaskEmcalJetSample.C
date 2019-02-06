AliAnalysisTaskEmcalJetSample* AddTaskEmcalJetSample(
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char* ncells             = "usedefault",
  const char *suffix             = ""
)
{
  return AliAnalysisTaskEmcalJetSample::AddTaskEmcalJetSample(ntracks,
      nclusters,
      ncells,
      suffix);
}
