// AddTaskJetExtractorPP.C
AliAnalysisTaskJetExtractorPP* AddTaskJetExtractorPP (TString trackArray, TString clusterArray, TString jetArray, TString rhoObject, Double_t jetRadius, const char* taskNameSuffix)
{
  AliRDHFJetsCutsVertex* cuts = new AliRDHFJetsCutsVertex("jetCuts");
  return AliAnalysisTaskJetExtractorPP::AddTaskJetExtractorPP(trackArray, clusterArray, jetArray, rhoObject, jetRadius, cuts, taskNameSuffix);
}
