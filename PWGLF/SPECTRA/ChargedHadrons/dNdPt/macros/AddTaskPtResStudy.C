#if defined(__CLING__)
#include "AliAnalysisTaskPtResStudy.h"
#endif

AliAnalysisTaskPtResStudy* AddTaskPtResStudy(const char* name = "TaskPtResStudy", const char* outfile = 0) {
  return AliAnalysisTaskPtResStudy::AddTaskPtResStudy(name, outfile);
}
