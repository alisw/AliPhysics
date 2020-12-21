AliAnalysisTaskIPCalib* AddTaskIPCalib(
  const char *ntracks            = "usedefault",
  TString pathToCorrFuncPscat = "",
  TString pathToCorrFuncNvtxContrib = "",
  const char* suffix = ""
)
{
  return AliAnalysisTaskIPCalib::AddTaskIPCalib(ntracks, pathToCorrFuncPscat, pathToCorrFuncNvtxContrib, suffix);
}
