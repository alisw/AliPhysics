AliAnalysisTaskIPCalib* AddTaskIPCalib(
  const char *ntracks            = "usedefault",
  const char *suffix             = ""
)
{
  return AliAnalysisTaskIPCalib::AddTaskIPCalib(ntracks, suffix);
}
