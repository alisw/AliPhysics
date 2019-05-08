#if defined(__CLING__)
#include "AliMultDepSpecAnalysisTask.h"
#endif

AliMultDepSpecAnalysisTask* AddTask_mkrueger_MultDepSpec(TString controlstring, Int_t cutModeLow = 100, Int_t cutModeHigh = 121)
{
  return AliMultDepSpecAnalysisTask::AddTaskMultDepSpec(controlstring, cutModeLow, cutModeHigh);
}
