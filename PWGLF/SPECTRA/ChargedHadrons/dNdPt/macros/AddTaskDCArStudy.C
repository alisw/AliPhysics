#if defined(__CLING__)
#include "AliAnalysisTaskDCArStudy.h"
#endif

AliAnalysisTaskDCArStudy* AddTaskDCArStudy(const char* name = "TaskDCArStudy", const char* outfile = 0) {
  return AliAnalysisTaskDCArStudy::AddTaskDCArStudy(name, outfile);
}
