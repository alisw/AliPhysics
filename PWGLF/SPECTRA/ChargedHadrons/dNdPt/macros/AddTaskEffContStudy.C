#if defined(__CLING__)
#include "AliAnalysisTaskEffContStudy.h"
#endif

AliAnalysisTaskEffContStudy* AddTaskEffContStudy(const char* name = "TaskEffContStudy", const char* outfile = 0) {
  return AliAnalysisTaskEffContStudy::AddTaskEffContStudy(name, outfile);
}
