#if defined(__CLING__)
#include "AliMCSpectraWeightsAnalysisTask.h"
#endif

AliMCSpectraWeightsAnalysisTask* AddTaskWeights(const char* collisionSystem = "pp",
const char* previousTrain = 0,
const char* name = "MCWeightsAnalysisTask",
const char* outfile = 0) {
  return AliMCSpectraWeightsAnalysisTask::AddTaskWeights(collisionSystem, previousTrain, name, outfile);
}
