#ifndef __CINT__
#include "AliAnalysisTaskNanoAODnormalisation.h"
#include <string>
#endif

AliAnalysisTaskNanoAODnormalisation* AddTaskNanoAODnormalisation(std::string name = "NanoAODnormalisation") {
  return AliAnalysisTaskNanoAODnormalisation::AddTask(name); 
}
