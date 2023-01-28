
AliAnalysisTaskEmbeddingJetWithEP* AddTaskEmbeddingJetWithEP(
  const char *ntracks = "usedefault", 
  const char *nclusters = "usedefault",
  const char* ncells = "usedefault", 
  const char *suffix = ""
)
{
  return AliAnalysisTaskEmbeddingJetWithEP::AddTaskEmbeddingJetWithEP(
      ntracks, 
      nclusters,
      ncells, 
      suffix);
}
