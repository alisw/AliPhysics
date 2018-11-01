// AddTaskJetExtractor.C
AliAnalysisTaskJetExtractor* AddTaskJetExtractor(const char* configFile,
                                                 const char* taskNameSuffix = 0)
{
  AliRDHFJetsCutsVertex* cuts = new AliRDHFJetsCutsVertex("jetCuts");
  return AliAnalysisTaskJetExtractor::AddTaskJetExtractor(configFile, cuts, taskNameSuffix);
}
