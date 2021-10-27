AliAnalysisTaskRawJetWithEP* AddTaskRawJetWithEP(
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char* ncells             = "usedefault",
  const char *suffix             = "new"
)
{
  return AliAnalysisTaskRawJetWithEP::AddTaskRawJetWithEP(ntracks,
      nclusters,
      ncells,
      suffix);
}

