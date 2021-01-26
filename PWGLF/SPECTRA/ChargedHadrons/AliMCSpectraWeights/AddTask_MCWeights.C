#if defined(__CLING__)
#include "AliMCWeightsTask.h"
#endif

AliMCWeightsTask* AddTask_MCWeights() {
    return AliMCWeightsTask::AddTaskAliMCWeightsTask();
}
