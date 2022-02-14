#include "AliAnalysisTaskAO2Dconverter.h"

AliAnalysisTaskAO2Dconverter* AddTaskAO2Dconverter(TString suffix = "")
{
  return AliAnalysisTaskAO2Dconverter::AddTask(suffix);
}
