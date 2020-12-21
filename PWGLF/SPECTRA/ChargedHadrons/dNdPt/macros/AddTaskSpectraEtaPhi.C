#if defined(__CLING__)
#include "AliAnalysisTaskSpectraEtaPhi.h"
#endif

AliAnalysisTaskSpectraEtaPhi* AddTaskSpectraEtaPhi(const char* name = "TaskSpectraEtaPhi", const char* outfile = 0) {
  return AliAnalysisTaskSpectraEtaPhi::AddTaskSpectra(name, outfile);
}
