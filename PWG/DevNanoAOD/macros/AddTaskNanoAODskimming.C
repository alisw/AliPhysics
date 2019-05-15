#ifndef __CINT__
#include "AliAnalysisTaskNanoAODskimming.h"
#include <string>
#endif

AliAnalysisTaskNanoAODskimming* AddTaskNanoAODskimming(std::string name = "NanoAODskimming") {
  return AliAnalysisTaskNanoAODskimming::AddTask(name); 
}
