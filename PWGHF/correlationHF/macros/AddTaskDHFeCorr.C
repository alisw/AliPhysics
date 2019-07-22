#include "AliAnalysisTaskDHFeCorr.h"
#include <string>

AliAnalysisTaskDHFeCorr* AddTaskDHFeCorr(std::string name, std::string config_file){
    return AliAnalysisTaskMyTask::AddTask(name,config_file);
}
