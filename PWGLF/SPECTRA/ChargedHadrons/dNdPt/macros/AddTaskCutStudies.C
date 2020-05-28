#if defined(__CLING__)
#include "AliAnalysisTaskCutStudies.h"
#endif

AliAnalysisTaskCutStudies* AddTaskCutStudies(const char* name = "TaskCutStudies") {
  return AliAnalysisTaskCutStudies::AddTaskCutStudies(name, outfile);
}
