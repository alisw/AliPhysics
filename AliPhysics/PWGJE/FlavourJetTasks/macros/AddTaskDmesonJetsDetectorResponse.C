/// \file AddTaskDmesonJetsDetectorResponse.C
/// \brief AddTask macro for the AliAnalysisTaskDmesonJetsDetectorResponse class.
///
/// AddTask macro for the AliAnalysisTaskDmesonJetsDetectorResponse class.
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Jun 10, 2016

AliAnalysisTaskDmesonJetsDetectorResponse* AddTaskDmesonJetsDetectorResponse(
    const char *ntracks    = "usedefault",
    const char *nclusters  = "usedefault",
    const char *nMCpart    = "",
    Int_t       nMaxTrees  = 2,
    const char *suffix     = ""
)
{  
  return AliAnalysisTaskDmesonJetsDetectorResponse::AddTaskDmesonJetsDetectorResponse(ntracks, nclusters, nMCpart, nMaxTrees, suffix);
}
