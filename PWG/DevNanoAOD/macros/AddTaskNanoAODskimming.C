#ifndef __CINT__
#include "AliAnalysisTaskNanoAODskimming.h"
#include <string>
#endif

AliAnalysisTaskNanoAODskimming* AddTaskNanoAODskimming(std::string name = "NanoAODSkimming") {
  return AliAnalysisTaskNanoAODskimming::AddTask(name); 
}
