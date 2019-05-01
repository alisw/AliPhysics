// AddTaskEventExtractor.C
AliAnalysisTaskEventExtractor* AddTaskEventExtractor (TString trackArray, TString particleArray, const char* taskNameSuffix = 0)
{
  return AliAnalysisTaskEventExtractor::AddTaskEventExtractor(trackArray, particleArray,taskNameSuffix);
}
