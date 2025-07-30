AliAnalysisTaskConvJet* AddTask_GammaConvJet(
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char* ncells             = "usedefault",
  const char *suffix             = "",
  double distToEMCBorder         = 0.,
  double distToSMEdges           = 0.,
  bool addV0sToJet               = false
)
{
  return AliAnalysisTaskConvJet::AddTask_GammaConvJet(ntracks,
      nclusters,
      ncells,
      suffix,
      distToEMCBorder,
      distToSMEdges,
      addV0sToJet);
}
