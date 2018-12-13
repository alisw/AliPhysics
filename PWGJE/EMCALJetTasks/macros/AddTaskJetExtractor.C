// AddTaskJetExtractor.C
AliAnalysisTaskJetExtractor* AddTaskJetExtractor (TString trackArray, TString jetArray, TString rhoObject, Double_t jetRadius, TString configFile, const char* taskNameSuffix)
{
  AliRDHFJetsCutsVertex* cuts = new AliRDHFJetsCutsVertex("jetCuts");
  return AliAnalysisTaskJetExtractor::AddTaskJetExtractor(trackArray, jetArray, rhoObject, jetRadius, configFile, cuts, taskNameSuffix);
}
