void ConfigLegoTrainPWGJE(int iFlag = 0){

  // Adds the Global Variables depending on the flag
  // 

  if(iFlag==1013){ // 10h
    AliAnalysisManager::SetGlobalStr("kJetRunPeriod","LHC10h");
  }

}
