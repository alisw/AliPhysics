#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TString.h>
#include "AliAnalysisTaskKaonXiCorrelation.h"
#endif

AliAnalysisTaskKaonXiCorrelation *AddTaskKaonXiCorrelation(bool isMC, TString tskname = "KaonXiCorrelation", TString suffix = "")
{
  return AliAnalysisTaskKaonXiCorrelation::AddTask(isMC, tskname, suffix);
}
