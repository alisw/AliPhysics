PWGJE::EMCALJetTasks::AliAnalysisTaskTrackSkim* AddTaskTrackSkim(const std::string& suffix = "")
{
  PWGJE::EMCALJetTasks::AliAnalysisTaskTrackSkim* task =
   PWGJE::EMCALJetTasks::AliAnalysisTaskTrackSkim::AddTaskTrackSkim(suffix);
  return task;
}
