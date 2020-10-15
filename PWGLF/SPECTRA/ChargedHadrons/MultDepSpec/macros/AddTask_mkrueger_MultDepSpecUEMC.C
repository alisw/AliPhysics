#if defined(__CLING__)
#include "AliMultDepSpecAnalysisTaskUE.h"
#endif

AliMultDepSpecAnalysisTaskUE* AddTask_mkrueger_MultDepSpecUEMC(
  const string& dataSet, int cutModeLow = 100, int cutModeHigh = 119, TString options = "",
  AliMultDepSpecAnalysisTaskUE::PhiRange phiRange = AliMultDepSpecAnalysisTaskUE::PhiRange::full,
  double ptLeadMIN = 0.01, const char* suffix = "")
{
  return AliMultDepSpecAnalysisTaskUE::AddTaskMultDepSpecUE(dataSet, cutModeLow, cutModeHigh,
                                                            options, true, phiRange, ptLeadMIN, suffix);
}
