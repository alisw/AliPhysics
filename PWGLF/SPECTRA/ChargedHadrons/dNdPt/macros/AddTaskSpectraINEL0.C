#if defined(__CLING__)
#include "AliAnalysisTaskSpectraINEL0.h"
#endif

AliAnalysisTaskSpectraINEL0* AddTaskSpectraINEL0(const char* name = "TaskSpectraINEL0", const char* outfile = 0) {
  return AliAnalysisTaskSpectraINEL0::AddTaskSpectraINEL0(name, outfile);
}
