#if defined(__CLING__)
#include "AliAnalysisTaskUEStudy.h"
#endif

AliAnalysisTaskUEStudy* AddTaskUEStudy(const char* name = "TaskUEStudy", const char* outfile = 0) {
  return AliAnalysisTaskUEStudy::AddTaskUEStudy(name, outfile);
}
