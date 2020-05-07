AliAnalysisTaskIPCalib* AddTaskIPCalib(
  const char *ntracks            = "usedefault",
  TString pathToCorrFunc = "",
  const char* suffix = ""
)
{
  return AliAnalysisTaskIPCalib::AddTaskIPCalib(ntracks, pathToCorrFunc, suffix);
}
