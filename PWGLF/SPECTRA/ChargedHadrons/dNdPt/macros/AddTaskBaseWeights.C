#if defined(__CLING__)
#include "AliAnalysisTaskBaseWeights.h"
#endif

AliAnalysisTaskBaseWeights* AddTaskBaseWeights(const char* name = "TaskBaseWeights", const char* outfile = 0) {
  return AliAnalysisTaskBaseWeights::AddTaskBaseWeights(name, outfile);
}
