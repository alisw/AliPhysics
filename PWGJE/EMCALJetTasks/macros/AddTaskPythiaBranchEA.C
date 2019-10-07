PWGJE::EMCALJetTasks::AliAnalysisTaskPythiaBranchEA* AddTaskPythiaBranchEA(const char *outputName = "pyparticles",
                                                     const char *pyfile ="",
                                                     const char *suffix = ""){

   return PWGJE::EMCALJetTasks::AliAnalysisTaskPythiaBranchEA::AddTaskPythiaBranchEA(
       outputName,
       pyfile,
       suffix);
}
