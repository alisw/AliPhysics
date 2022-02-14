AliAnalysisTaskFlowVectorCorrections* AddTaskFlowQnVectorCorrections(const char* configFilename)
{
  AliAnalysisTaskFlowVectorCorrections* task =
   PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHUtils::AddTaskFlowQnVectorCorrections(configFilename);
  return task;
}
