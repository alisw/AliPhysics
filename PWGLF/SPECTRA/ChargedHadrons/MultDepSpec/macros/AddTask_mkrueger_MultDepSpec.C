#if defined(__CLING__)
#include "AliMultDepSpecAnalysisTask.h"
#endif

AliMultDepSpecAnalysisTask* AddTask_mkrueger_MultDepSpec(TString controlstring, Int_t cutModeLow = 100, Int_t cutModeHigh = 119, Bool_t useDataDrivenCorrections = kFALSE, string pccTrainOutputPath = "",  Int_t pccSysFlag = 0,  Int_t secSysFlag = 0)
{
  return AliMultDepSpecAnalysisTask::AddTaskMultDepSpec(controlstring, cutModeLow, cutModeHigh, useDataDrivenCorrections,  pccTrainOutputPath, pccSysFlag, secSysFlag);
}
