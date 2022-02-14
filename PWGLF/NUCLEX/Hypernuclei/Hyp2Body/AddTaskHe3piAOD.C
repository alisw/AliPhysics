#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TString.h>
#include "AliAnalysisTaskHe3piAOD.h"
#endif

AliAnalysisTaskHe3piAOD *AddTaskHe3piAOD(bool isMC, TString tskname = "Hypertriton", TString suffix = "")
{
  return AliAnalysisTaskHe3piAOD::AddTask(isMC, tskname, suffix);
}
