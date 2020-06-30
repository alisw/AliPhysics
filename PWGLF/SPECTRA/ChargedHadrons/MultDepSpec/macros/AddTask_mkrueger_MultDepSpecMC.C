#if defined(__CLING__)
#include "AliMultDepSpecAnalysisTask.h"
#endif

AliMultDepSpecAnalysisTask* AddTask_mkrueger_MultDepSpecMC(const string& dataSet, int cutModeLow = 100, int cutModeHigh = 119, TString options = "")
{
  return AliMultDepSpecAnalysisTask::AddTaskMultDepSpec(dataSet, cutModeLow, cutModeHigh, options, true);
}
