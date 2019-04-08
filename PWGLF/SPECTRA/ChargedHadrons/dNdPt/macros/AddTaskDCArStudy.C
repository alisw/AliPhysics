#if defined(__CLING__)
#include "AliAnalysisTaskDCArStudy.h"
#endif

AliAnalysisTaskDCArStudy* AddTaskDCArStudy() {
  return AliAnalysisTaskDCArStudy::AddTaskDCArStudy("TaskDCArStudy");
}
