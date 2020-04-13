#if defined(__CLING__)
#include "AliAnalysisTaskSpectra.h"
#endif

AliAnalysisTaskSpectra* AddTaskSpectra(const char* name = "TaskSpectra", const char* outfile = 0) {
  return AliAnalysisTaskSpectra::AddTaskSpectra(name, outfile);
}
