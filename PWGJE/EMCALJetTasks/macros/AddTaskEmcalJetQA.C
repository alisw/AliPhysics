// AddTaskEmcalJetQA.C

AliAnalysisTaskEmcalJetQA* AddTaskEmcalJetQA(
  const char* ntracks            = "usedefault",
  const char* nclusters          = "usedefault",
  const char* ncells             = "usedefault",
  const char* suffix             = ""
)
{  
  return AliAnalysisTaskEmcalJetQA::AddTaskEmcalJetQA(ntracks, nclusters, ncells, "", suffix);
}
