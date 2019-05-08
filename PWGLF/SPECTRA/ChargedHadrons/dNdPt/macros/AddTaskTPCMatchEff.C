#if defined(__CLING__)
#include "AliAnalysisTaskTPCMatchEff.h"
#endif

AliAnalysisTaskTPCMatchEff* AddTaskTPCMatchEff(const char* name = "TaskTPCMatchEff", const char* outfile = 0) {
  return AliAnalysisTaskTPCMatchEff::AddTaskTPCMatchEff(name, outfile);
}
