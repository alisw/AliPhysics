#if defined(__CLING__)
#include "AliMultDepSpecAnalysisTaskUE.h"
#endif

AliMultDepSpecAnalysisTaskUE* AddTask_mkrueger_MultDepSpecUE(const string& dataSet,
                                                             int cutModeLow = 100,
                                                             int cutModeHigh = 119,
                                                             TString options = "", bool isUE = true,
                                                             float ptLeadMIN = 0.01)
{
  return AliMultDepSpecAnalysisTaskUE::AddTaskMultDepSpecUE(dataSet, cutModeLow, cutModeHigh,
                                                            options, false, isUE, ptLeadMIN);
}
