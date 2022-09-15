#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TString.h>
#include "AliAnalysisTaskHe3piKF.h"
#endif

AliAnalysisTaskHe3piKF *AddTaskHe3piKF(bool isMC, TString tskname = "HyperTask", TString suffix = "")
{
  return AliAnalysisTaskHe3piKF::AddTask(isMC, tskname, suffix);
}
