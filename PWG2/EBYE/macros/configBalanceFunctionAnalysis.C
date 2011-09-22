//__________________________________________________//
AliBalance *GetBalanceFunctionObject(const char* analysisLevel = "ESD", Bool_t bShuffle = kFALSE) {
  //Function to setup the AliProtonAnalysis object and return it
  AliBalance *gBalance = new AliBalance();
  gBalance->SetAnalysisLevel(analysisLevel);
  gBalance->SetShuffle(bShuffle);

  //Set all the same
  //gBalance->SetNumberOfBins(-1,9);
  //gBalance->SetInterval(-1,0.0,0.9);

  //Set all analyses separately
  //Rapidity
  //gBalance->SetNumberOfBins(0,64);
  gBalance->SetInterval(AliBalance::kRapidity,-0.8,0.8,64,0.0,1.6);
  
  //Eta
  //gBalance->SetNumberOfBins(1,64);
  gBalance->SetInterval(AliBalance::kEta,-0.8,0.8,64,0.0,1.6);
  
  //Qlong
  //gBalance->SetNumberOfBins(2,64);
  gBalance->SetInterval(AliBalance::kQlong,-1,1,100,0.0,2.0);
  
  //Qout
  //gBalance->SetNumberOfBins(3,64);
  gBalance->SetInterval(AliBalance::kQout,-1,1,100,0.0,2.0);
  
  //Qside
  //gBalance->SetNumberOfBins(4,64);
  gBalance->SetInterval(AliBalance::kQside,-1,1,100,0.0,2.0);
  
  //Qinv
  //gBalance->SetNumberOfBins(5,64);
  gBalance->SetInterval(AliBalance::kQinv,-1,1,100,0.0,2.0);
  
  //Phi
  //gBalance->SetNumberOfBins(6,180);
  gBalance->SetInterval(0,360,90,0,180.0);

  //Init the histograms
  gBalance->InitHistograms();
  
  return gBalance;
}

//__________________________________________________//
AliESDtrackCuts *GetTrackCutsObject() {
  // only used for ESDs
  // Function to setup the AliESDtrackCuts object and return it
  AliESDtrackCuts *cuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  cuts->SetMinNClustersTPC(70);
  
  cuts->SetPtRange(0.3,1.5);
  cuts->SetEtaRange(-0.8,0.8);
  cuts->DefineHistograms(1);
  //cuts->SaveHistograms("trackCuts");

  return cuts;
}

