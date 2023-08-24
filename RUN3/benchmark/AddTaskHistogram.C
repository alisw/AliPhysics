#include "AliAnalysisTaskHistogram.h"

AliAnalysisTaskHistogram* AddTaskHistogram(TString suffix = "")
{
  return AliAnalysisTaskHistogram::AddTask(suffix);
}
