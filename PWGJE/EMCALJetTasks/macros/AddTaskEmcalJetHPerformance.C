PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHPerformance* AddTaskEmcalJetHPerformance(
    const char * suffix = ""
)
{  
  PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHPerformance * task = PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHPerformance::AddTaskEmcalJetHPerformance(suffix);
  return task;
}
