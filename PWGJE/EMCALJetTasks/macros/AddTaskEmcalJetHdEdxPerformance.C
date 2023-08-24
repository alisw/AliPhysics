PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHdEdxPerformance* AddTaskEmcalJetHdEdxPerformance(
    const char * suffix = ""
)
{  
  PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHdEdxPerformance * task = PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetHdEdxPerformance::AddTaskEmcalJetHdEdxPerformance(suffix);
  return task;
}
