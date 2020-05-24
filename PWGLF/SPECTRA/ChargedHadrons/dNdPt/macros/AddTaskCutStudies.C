#if defined(__CLING__)
#include "AliAnalysisTaskCutStudies.h"
#endif

AliAnalysisTaskCutStudies* AddTaskCutStudies(const char* name = "TaskCutStudies", const char* outfile = 0) {
  return AliAnalysisTaskCutStudies::AddTaskCutStudies(name, outfile);
}
