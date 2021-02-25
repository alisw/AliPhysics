#if defined(__CLING__)
#include "AliMultDepSpecAnalysisTaskOLD.h"
#endif

AliMultDepSpecAnalysisTaskOLD* AddTask_mkrueger_MultDepSpec(const string& dataSet,
                                                         int cutModeLow = 100,
                                                         int cutModeHigh = 119,
                                                         TString options = "")
{
  return AliMultDepSpecAnalysisTaskOLD::AddTaskMultDepSpec(dataSet, cutModeLow, cutModeHigh, options);
}
