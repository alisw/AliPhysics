AliAnalysisTaskRawJetSpec* AddTaskRawJetSpec(
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char* ncells             = "usedefault",
  const char *suffix             = ""
)
{
  return AliAnalysisTaskRawJetSpec::AddTaskRawJetSpec(ntracks,
      nclusters,
      ncells,
      suffix);
}
