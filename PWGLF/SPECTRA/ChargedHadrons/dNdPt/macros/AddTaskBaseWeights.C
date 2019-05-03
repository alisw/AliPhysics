#if defined(__CLING__)
#include "AliAnalysisTaskBaseWeights.h"
#endif

AliAnalysisTaskBaseWeights* AddTaskBaseWeights(const char* name = "TaskBaseWeights", const char* outfile = 0, const char* collisionSystem = 0, Int_t sysFlag = 0, const char* prevTrainOutputPath = 0) {
  return AliAnalysisTaskBaseWeights::AddTaskBaseWeights(name, outfile, collisionSystem, sysFlag, prevTrainOutputPath);
}
