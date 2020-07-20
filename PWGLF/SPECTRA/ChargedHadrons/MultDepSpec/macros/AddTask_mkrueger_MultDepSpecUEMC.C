#if defined(__CLING__)
#include "AliMultDepSpecAnalysisTaskUE.h"
#endif

AliMultDepSpecAnalysisTaskUE* AddTask_mkrueger_MultDepSpecUEMC(const string& dataSet, int cutModeLow = 100, int cutModeHigh = 119, TString options = "")
{
  return AliMultDepSpecAnalysisTaskUE::AddTaskMultDepSpecUE(dataSet, cutModeLow, cutModeHigh, options, true);
}
