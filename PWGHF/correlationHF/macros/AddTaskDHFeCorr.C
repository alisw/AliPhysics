#include "AliAnalysisTaskDHFeCorr.h"
#include <string>

AliAnalysisTaskDHFeCorr* AddTaskDHFeCorr(std::string name, std::string config_file){
    return AliAnalysisTaskDHFeCorr::AddTask(name,config_file);
}
