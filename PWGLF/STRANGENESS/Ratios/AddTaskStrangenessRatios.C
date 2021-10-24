#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TString.h>
#include "AliAnalysisTaskStrangenessRatios.h"
#endif

AliAnalysisTaskStrangenessRatios *AddTaskStrangenessRatios(bool isMC, TString tskname = "StrangenessRatios", TString suffix = "")
{
  return AliAnalysisTaskStrangenessRatios::AddTask(isMC, tskname, suffix);
}
