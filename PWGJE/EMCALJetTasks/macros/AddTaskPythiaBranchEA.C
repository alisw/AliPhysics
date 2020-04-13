PWGJE::EMCALJetTasks::AliAnalysisTaskPythiaBranchEA* AddTaskPythiaBranchEA(const char *outputName = "pyparticles",
                                                                           const char *pyfilepath ="",
                                                                           const char *pyfilemask ="",
                                                                           const char *suffix = ""){

   return PWGJE::EMCALJetTasks::AliAnalysisTaskPythiaBranchEA::AddTaskPythiaBranchEA(
       outputName,
       pyfilepath,
       pyfilemask,
       suffix);
}
