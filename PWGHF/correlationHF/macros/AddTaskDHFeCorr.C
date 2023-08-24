#include "AliAnalysisTaskDHFeCorr.h"
#include <string>

AliAnalysisTaskDHFeCorr *AddTaskDHFeCorr(std::string name, std::string config_file, Int_t trigger = AliVEvent::kINT7) {
    return AliAnalysisTaskDHFeCorr::AddTask(name, config_file, trigger);
}
