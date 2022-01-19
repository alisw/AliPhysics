#if defined(__CLING__)
#include "AliMultDepSpecAnalysisTask.h"
#endif

AliMultDepSpecAnalysisTask* AddTask_mkrueger_MultDepSpecMC(const string& dataSet,
                                                           TString options = "",
                                                           int cutModeLow = 100,
                                                           int cutModeHigh = 121)
{
  return AliMultDepSpecAnalysisTask::AddTaskMultDepSpec(dataSet, options, cutModeLow, cutModeHigh, true);
}
