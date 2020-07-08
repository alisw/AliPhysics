AliAnalysisTaskJetVn* AddTaskJetVn(
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char* ncells             = "usedefault",
  const char *suffix             = "new"
)
{
  return AliAnalysisTaskJetVn::AddTaskEmcalJetSample(
      ntracks,
      nclusters,
      ncells,
      suffix);
}
