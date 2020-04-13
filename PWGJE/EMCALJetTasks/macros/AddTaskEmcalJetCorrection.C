// AddTaskEmcalJetCorrection.C
AliAnalysisTaskEmcalJetCorrection* AddTaskEmcalJetCorrection(TString modelName, TString trackArray, TString jetArray, TString rhoObject, Double_t jetRadius, const char* taskNameSuffix)
{
  return AliAnalysisTaskEmcalJetCorrection::AddTaskEmcalJetCorrection(modelName, trackArray, jetArray, rhoObject, jetRadius, taskNameSuffix);
}
