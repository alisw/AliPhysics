//__________________________________________________//
AliBalance *GetBalanceFunctionObject(const char* analysisLevel = "ESD") {
  //Function to setup the AliProtonAnalysis object and return it
  AliBalance *gBalance = new AliBalance();
  gBalance->SetAnalysisLevel(analysisLevel);
  gBalance->SetAnalysisType(1);
  gBalance->SetNumberOfBins(18);
  gBalance->SetInterval(-0.9,0.9);

  return gBalance;
}
