#if defined(__CLING__)
#include "AliAnalysisTaskBaseWeights.h"
#endif

AliAnalysisTaskBaseWeights* AddTaskBaseWeights_pp(const char* name = "TaskBaseWeights", const char* outfile = 0, const char* collisionSystem = "pp", Int_t sysFlag = 0, const char* prevTrainOutputPath = "alien:///alice/sim/2017/LHC17l3b_cent/282247/PWGLF/LF_pp_MC/1110_20190427-1332_child_1/AnalysisResults.root") {
  return AliAnalysisTaskBaseWeights::AddTaskBaseWeights(name, outfile, collisionSystem, sysFlag, prevTrainOutputPath);
}
    