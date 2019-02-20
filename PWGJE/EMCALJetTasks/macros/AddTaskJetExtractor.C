// AddTaskJetExtractor.C
AliAnalysisTaskJetExtractor* AddTaskJetExtractor (TString trackArray, TString clusterArray, TString jetArray, TString rhoObject, Double_t jetRadius, const char* taskNameSuffix)
{
  AliRDHFJetsCutsVertex* cuts = new AliRDHFJetsCutsVertex("jetCuts");
  return AliAnalysisTaskJetExtractor::AddTaskJetExtractor(trackArray, clusterArray, jetArray, rhoObject, jetRadius, cuts, taskNameSuffix);
}
