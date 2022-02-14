AliAnalysisTaskConvJet* AddTask_GammaConvJet(
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char* ncells             = "usedefault",
  const char *suffix             = ""
)
{
  return AliAnalysisTaskConvJet::AddTask_GammaConvJet(ntracks,
      nclusters,
      ncells,
      suffix);
}
