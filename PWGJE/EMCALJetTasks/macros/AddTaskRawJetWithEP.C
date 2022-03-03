
AliAnalysisTaskRawJetWithEP1* AddTaskRawJetWithEP1(
  const char *ntracks = "usedefault", 
  const char *nclusters = "usedefault",
  const char* ncells = "usedefault", 
  const char *suffix = ""
)
{
  return AliAnalysisTaskRawJetWithEP1::AddTaskRawJetWithEP1(
      ntracks, 
      nclusters,
      ncells, 
      suffix);
}