AliAnalysisTaskConvJet* AddTask_GammaConvJet(
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char* ncells             = "usedefault",
  const char *suffix             = "",
  Int_t IsMC			 = 0,
  Int_t NContainers              = 1
)
{
  return AliAnalysisTaskConvJet::AddTask_GammaConvJet(ntracks,
      nclusters,
      ncells,
      suffix,
      IsMC,
      NContainers);
}
