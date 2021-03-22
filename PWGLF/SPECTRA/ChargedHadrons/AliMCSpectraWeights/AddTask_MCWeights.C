#if defined(__CLING__)
#include "AliMCWeightsTask.h"
#else
enum MCGeneratorType {
    NONE=-1,
    PP_PYTHIA=0,
    PPB_EPOS,
    PBPB_HIJING,
};
#endif

AliMCWeightsTask* AddTask_MCWeights(MCGeneratorType gen = MCGeneratorType::NONE, const char* collType = "pp") {
    return AliMCWeightsTask::AddTaskAliMCWeightsTask(gen, collType);
}
