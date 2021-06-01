#if defined(__CLING__)
#include "AliMultDepSpecAnalysisTask.h"
#endif

AliMultDepSpecAnalysisTask* AddTask_mkrueger_MultDepSpecMC(const string& dataSet,
                                                           int cutModeLow = 99,
                                                           int cutModeHigh = 123,
                                                           TString options = "")
{
  return AliMultDepSpecAnalysisTask::AddTaskMultDepSpec(dataSet, cutModeLow, cutModeHigh, options, true);
}
