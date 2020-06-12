#if defined(__CLING__)
#include "AliMultDepSpecAnalysisTask.h"
#endif

AliMultDepSpecAnalysisTask* AddTask_mkrueger_MultDepSpec(string dataSet, TString options, int cutModeLow = 100, int cutModeHigh = 119)
{
  return AliMultDepSpecAnalysisTask::AddTaskMultDepSpec(dataSet, options, cutModeLow, cutModeHigh);
}
