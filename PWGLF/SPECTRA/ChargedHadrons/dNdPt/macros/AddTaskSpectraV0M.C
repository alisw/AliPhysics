#if defined(__CLING__)
#include "AliAnalysisTaskSpectraV0M.h"
#endif

AliAnalysisTaskSpectraV0M* AddTaskSpectraV0M(const char* name = "TaskSpectraV0M", const char* outfile = 0) {
  return AliAnalysisTaskSpectraV0M::AddTaskSpectraV0M(name, outfile);
}
