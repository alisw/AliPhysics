#if defined(__CLING__)
#include "AliAnalysisTaskMKBase.h"
#endif

AliAnalysisTaskMKBase* AliAnalysisTaskMKBase(const char* name = "TaskMKBase", const char* outfile = 0) {
  return AliAnalysisTaskMKBase::AddTaskMKBase(name, outfile);
}
