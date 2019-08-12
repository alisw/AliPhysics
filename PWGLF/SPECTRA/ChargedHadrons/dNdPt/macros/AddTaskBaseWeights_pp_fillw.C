#if defined(__CLING__)
#include "AliAnalysisTaskBaseWeights.h"
#endif

AliAnalysisTaskBaseWeights* AddTaskBaseWeights_pp_fillw(const char* name = "TaskBaseWeights", const char* outfile = 0, const char* collisionSystem = "pp", Int_t sysFlag = 0, const char* prevTrainOutputPath = 0) {
  return AliAnalysisTaskBaseWeights::AddTaskBaseWeights(name, outfile, collisionSystem, sysFlag, prevTrainOutputPath);
}
