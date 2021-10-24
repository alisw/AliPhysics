#if defined(__CLING__)
#include "AliMultDepSpecAnalysisTaskUE.h"
#endif

AliMultDepSpecAnalysisTaskUE* AddTask_mkrueger_MultDepSpecUE(const string& dataSet,
                                                             int cutModeLow = 100,
                                                             int cutModeHigh = 119,
                                                             TString options = "",
                                                             AliMultDepSpecAnalysisTaskUE::PhiRange phiRange = AliMultDepSpecAnalysisTaskUE::PhiRange::full,
                                                             float ptLeadMIN = 0.01,
                                                             const char* suffix = "")
{
  return AliMultDepSpecAnalysisTaskUE::AddTaskMultDepSpecUE(dataSet, cutModeLow, cutModeHigh,
                                                            options, false, phiRange, ptLeadMIN,  suffix);
}
